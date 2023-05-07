### to run the AR method (traditional and ML), IV ranges from weak to strong
library(drf)
library(ivmodel)

n <- 1000  # {1000, 3000, 5000}
a <- 0  # {0 ~ 1}
grid <- 0.01  # {>0 ~ 0.1}

### simulation ----------------------------------------------------------------

sim_func <- function(n, a, grid=0.01) {
    ### generate data -------------------------------------------------------------
    
    # dimension
    p <- 10
    # sample size
    n <- n
    # interaction value
    inter_val <- 0
    # the IV strength
    a <- a
    # violation strength
    tau <- 0
    f_1 <- function(z, a) {
        0.5*z + a * ( z^2 + 0.125 * z^4) - 25 / 12
    }
    rho <- 0.5
    Cov <- stats::toeplitz(rho^c(0 : p))
    mu <- rep(0, p + 1)
    # true effect
    beta <- 1
    alpha <- as.matrix(rep(-0.3, p))
    gamma <- as.matrix(rep(0.2, p))
    inter <- as.matrix(c(rep(inter_val, 5),rep(0, p - 5)))
    
    # generate the data
    mu_error <- rep(0,2)
    Cov_error <- matrix(c(1, 0.5, 0.5,1), 2, 2)
    Error <- MASS::mvrnorm(n, mu_error, Cov_error)
    delta <- Error[, 1]
    eps <- Error[, 2]
    W_original <- MASS::mvrnorm(n, mu, Cov)
    W <- pnorm(W_original)
    # instrument variable
    Z <- W[, 1]
    # baseline covariates
    X <- W[, -1]
    # generate the treatment variable D
    D <- f_1(Z, a) + X %*% alpha + Z * X %*% inter + delta
    f <- D - Error[,1]
    # generate the outcome variable Y
    Y <- D * beta + tau * Z + X %*% gamma + eps
    
    ### Stage 1: DRF --------------------------------------------------------------
    
    n1 <- round(n * 2 / 3)
    A1_ind <- 1:n1
    inputX_train <- cbind(Z, X)[-A1_ind, ]
    inputY_train <- cbind(Y, D)[-A1_ind, ]
    inputX_test <- cbind(Z, X)[A1_ind, ]
    inputY_test <- cbind(Y, D)[A1_ind, ]
    forest <- drf(inputX_train, inputY_train)
    weight <- as.matrix(get_sample_weights(forest, newdata=inputX_test))
    
    ### Stage 2: AR ---------------------------------------------------------------
    
    beta_cand_test <- function(beta_cand, Y, D, X, A1_ind, weight) {
        Y_A1 <- as.matrix(Y[A1_ind, ])
        D_A1 <- as.matrix(D[A1_ind, ])
        X_A1 <- as.matrix(X[A1_ind, ])
        Y_A2 <- as.matrix(Y[-A1_ind, ])
        D_A2 <- as.matrix(D[-A1_ind, ])
        X_A2 <- as.matrix(X[-A1_ind, ])
        Y_rep <- as.matrix(weight %*% Y_A2)
        D_rep <- as.matrix(weight %*% D_A2)
        X_rep <- as.matrix(weight %*% X_A2)
        
        Tstat_nume <- sum(resid(lm(Y_rep-beta_cand*D_rep ~ X_rep))^2)
        weight_resid <- resid(lm(weight ~ X_rep))
        trace <- sum(diag(tcrossprod(weight_resid)))
        e_hat1 <- Y_A1 - Y_rep
        delta_hat <- D_A1 - D_rep
        SigmaSq_e1 <- mean(e_hat1^2)
        SigmaSq_delta <- mean(delta_hat^2)
        Sigma_e1delta <- mean(e_hat1 * delta_hat)
        Tstat_deno1 <- trace * (SigmaSq_e1 - 2*beta_cand*Sigma_e1delta + beta_cand^2*SigmaSq_delta)
        Tstat1 <- Tstat_nume / Tstat_deno1
        thol <- qchisq(0.95, df=1)
        result1 <- as.numeric(Tstat1 <= thol)
        
        # the alternative estimation for e
        e_hat2 <- resid(lm(Y_A1-Y_rep ~ X_A1))
        SigmaSq_e2 <- mean(e_hat2^2)
        Sigma_e2delta <- mean(e_hat2 * delta_hat)
        Tstat_deno2 <- trace * (SigmaSq_e2 - 2*beta_cand*Sigma_e2delta + beta_cand^2*SigmaSq_delta)
        Tstat2 <- Tstat_nume / Tstat_deno2
        result2 <- as.numeric(Tstat2 <= thol)
        
        # another test statistics Tstat_prime (numerators)
        Tstat_prime_nume1 <- Tstat_nume - Tstat_deno1
        
        # another test statistics Tstat_prime with alternative e (numerators)
        Tstat_prime_nume2 <- Tstat_nume - Tstat_deno2
        
        
        returnResults <- c(beta = beta_cand
                           , result1 = result1
                           , result2 = result2
                           , Tstat1 = Tstat1
                           , Tstat2 = Tstat2
                           , Tstat_prime_nume1 = Tstat_prime_nume1
                           , Tstat_prime_nume2 = Tstat_prime_nume2
                           # for calculating the variance -- denominator
                           , Tstat_nume = Tstat_nume)
        returnResults
    }
    CI_set <- NULL
    niter <- 1
    beta_cand_left <- beta_cand_right <- beta
    results_left <- results_right <- c(beta, 999, 999)
    while ((sum(results_left[2:3]) > 0 || sum(results_right[2:3]) > 0) && niter <= 150) {
        # print(c(niter=niter, left=beta_cand_left, right=beta_cand_right))
        if (sum(results_left[2:3]) > 0) {
            results_left <- beta_cand_test(beta_cand_left, Y, D, X, A1_ind, weight)
        }
        if (sum(results_right[2:3]) > 0 && beta_cand_left != beta_cand_right) {
            results_right <- beta_cand_test(beta_cand_right, Y, D, X, A1_ind, weight)
        }
        if (length(which(CI_set[, 1] == results_left[1])) == 0) {
            CI_set <- rbind(CI_set, results_left)
        }
        if (length(which(CI_set[, 1] == results_right[1])) == 0) {
            CI_set <- rbind(CI_set, results_right)
        }
        niter <- niter + 1
        beta_cand_left <- beta_cand_left - grid
        beta_cand_right <- beta_cand_right + grid
    }
    # continue calculating Tstat_prime
    Tstat_prime_deno <- var(CI_set[, 8])
    Tstat_prime1 <- CI_set[, 6] / Tstat_prime_deno
    Tstat_prime2 <- CI_set[, 7] / Tstat_prime_deno
    result3 <- qnorm(0.025) <= Tstat_prime1 & Tstat_prime1 <= qnorm(0.975)
    result4 <- qnorm(0.025) <= Tstat_prime2 & Tstat_prime2 <= qnorm(0.975)
    CI_set <- cbind(CI_set, result3, result4)
    
    CI_set <- CI_set[order(CI_set[, 1]), ]
    length1_set <- CI_set[max(which(CI_set[, 2]==1)), 1] - 
        CI_set[min(which(CI_set[, 2]==1)), 1]
    cover1_set <- CI_set[which(CI_set[, 1]==beta), 2]
    length2_set <- CI_set[max(which(CI_set[, 3]==1)), 1] -
        CI_set[min(which(CI_set[, 3]==1)), 1]
    cover2_set <- CI_set[which(CI_set[, 1]==beta), 3]
    length3_set <- CI_set[max(which(CI_set[, 9]==1)), 1] -
        CI_set[min(which(CI_set[, 9]==1)), 1]
    cover3_set <- CI_set[which(CI_set[, 1]==beta), 9]
    length4_set <- CI_set[max(which(CI_set[, 10]==1)), 1] -
        CI_set[min(which(CI_set[, 10]==1)), 1]
    cover4_set <- CI_set[which(CI_set[, 1]==beta), 10]
    Tstat1_quantile <- quantile(CI_set[, 4])
    Tstat2_quantile <- quantile(CI_set[, 5])
    Tstat_prime1_quantile <- quantile(Tstat_prime1)
    Tstat_prime2_quantile <- quantile(Tstat_prime2)
    
    ### AR test in library --------------------------------------------------------
    
    Y_A1 <- as.matrix(Y[A1_ind, ])
    D_A1 <- as.matrix(D[A1_ind, ])
    X_A1 <- as.matrix(X[A1_ind, ])
    Y_A2 <- as.matrix(Y[-A1_ind, ])
    D_A2 <- as.matrix(D[-A1_ind, ])
    X_A2 <- as.matrix(X[-A1_ind, ])
    Y_rep <- as.matrix(weight %*% Y_A2)
    D_rep <- as.matrix(weight %*% D_A2)
    X_rep <- as.matrix(weight %*% X_A2)
    
    reg <- ivmodel(Y_A1, D_A1, D_rep, X_A1, k=1)
    beta_hat <- coef(reg)['TSLS', 'Estimate']
    sd_hat <- coef(reg)['TSLS', 'Std. Error']
    CI_ar <- reg$AR$ci
    length <- sum(CI_ar[, 2] - CI_ar[, 1])
    cover <- any(CI_ar[, 1] <= beta & beta <= CI_ar[, 2])
    
    ### return --------------------------------------------------------------------
    
    returnValue <- c(
        beta_true = beta
        , length1_set = as.numeric(length1_set)
        , length2_set = as.numeric(length2_set)
        , length3_set = as.numeric(length3_set)
        , length4_set = as.numeric(length4_set)
        , length = length
        , cover1_set = as.numeric(cover1_set)
        , cover2_set = as.numeric(cover2_set)
        , cover3_set = as.numeric(cover3_set)
        , cover4_set = as.numeric(cover4_set)
        , cover = cover
        , Tstat1_quantile = Tstat1_quantile
        , Tstat2_quantile = Tstat2_quantile
        , Tstat_prime1_quantile = Tstat_prime1_quantile
        , Tstat_prime2_quantile = Tstat_prime2_quantile
    )
    returnValue
    
}

nsim <- 300
total <- NULL
for(i in 1:nsim) {
    print(c(n, a, i))
    result <- sim_func(n, a)
    total <- rbind(total, result)
}
total

### generate data -------------------------------------------------------------

# dimension
p <- 10
# sample size
n <- n
# interaction value
inter_val <- 0
# the IV strength
a <- a
# violation strength
tau <- 0
f_1 <- function(z, a) {
  0.5*z + a * ( z^2 + 0.125 * z^4) - 25 / 12
}
rho <- 0.5
Cov <- stats::toeplitz(rho^c(0 : p))
mu <- rep(0, p + 1)
# true effect
beta <- 1
alpha <- as.matrix(rep(-0.3, p))
gamma <- as.matrix(rep(0.2, p))
inter <- as.matrix(c(rep(inter_val, 5),rep(0, p - 5)))

# generate the data
mu_error <- rep(0,2)
Cov_error <- matrix(c(1, 0.5, 0.5,1), 2, 2)
Error <- MASS::mvrnorm(n, mu_error, Cov_error)
delta <- Error[, 1]
eps <- Error[, 2]
W_original <- MASS::mvrnorm(n, mu, Cov)
W <- pnorm(W_original)
# instrument variable
Z <- W[, 1]
# baseline covariates
X <- W[, -1]
# generate the treatment variable D
D <- f_1(Z, a) + X %*% alpha + Z * X %*% inter + delta
f <- D - Error[,1]
# generate the outcome variable Y
Y <- D * beta + tau * Z + X %*% gamma + eps

### Stage 1: DRF --------------------------------------------------------------

n1 <- round(n * 2 / 3)
A1_ind <- 1:n1
inputX_train <- cbind(Z, X)[-A1_ind, ]
inputY_train <- cbind(Y, D)[-A1_ind, ]
inputX_test <- cbind(Z, X)[A1_ind, ]
inputY_test <- cbind(Y, D)[A1_ind, ]
forest <- drf(inputX_train, inputY_train)
weight <- as.matrix(get_sample_weights(forest, newdata=inputX_test))

### Stage 2: AR ---------------------------------------------------------------

beta_cand_test <- function(beta_cand, Y, D, X, A1_ind, weight) {
    Y_A1 <- as.matrix(Y[A1_ind, ])
    D_A1 <- as.matrix(D[A1_ind, ])
    X_A1 <- as.matrix(X[A1_ind, ])
    Y_A2 <- as.matrix(Y[-A1_ind, ])
    D_A2 <- as.matrix(D[-A1_ind, ])
    X_A2 <- as.matrix(X[-A1_ind, ])
    Y_rep <- as.matrix(weight %*% Y_A2)
    D_rep <- as.matrix(weight %*% D_A2)
    X_rep <- as.matrix(weight %*% X_A2)
    
    Tstat_nume <- sum(resid(lm(Y_rep-beta_cand*D_rep ~ X_rep))^2)
    weight_resid <- resid(lm(weight ~ X_rep))
    trace <- sum(diag(tcrossprod(weight_resid)))
    e_hat1 <- Y_A1 - Y_rep
    delta_hat <- D_A1 - D_rep
    SigmaSq_e1 <- mean(e_hat1^2)
    SigmaSq_delta <- mean(delta_hat^2)
    Sigma_e1delta <- mean(e_hat1 * delta_hat)
    Tstat_deno1 <- trace * (SigmaSq_e1 - 2*beta_cand*Sigma_e1delta + beta_cand^2*SigmaSq_delta)
    Tstat1 <- Tstat_nume / Tstat_deno1
    thol <- qchisq(0.95, df=1)
    result1 <- as.numeric(Tstat1 <= thol)
    
    # the alternative estimation for e
    e_hat2 <- resid(lm(Y_A1-Y_rep ~ X_A1))
    SigmaSq_e2 <- mean(e_hat2^2)
    Sigma_e2delta <- mean(e_hat2 * delta_hat)
    Tstat_deno2 <- trace * (SigmaSq_e2 - 2*beta_cand*Sigma_e2delta + beta_cand^2*SigmaSq_delta)
    Tstat2 <- Tstat_nume / Tstat_deno2
    result2 <- as.numeric(Tstat2 <= thol)
    
    # another test statistics Tstat_prime (numerators)
    Tstat_prime_nume1 <- Tstat_nume - Tstat_deno1

    # another test statistics Tstat_prime with alternative e (numerators)
    Tstat_prime_nume2 <- Tstat_nume - Tstat_deno2

    
    returnResults <- c(beta = beta_cand
                       , result1 = result1
                       , result2 = result2
                       , Tstat1 = Tstat1
                       , Tstat2 = Tstat2
                       , Tstat_prime_nume1 = Tstat_prime_nume1
                       , Tstat_prime_nume2 = Tstat_prime_nume2
                       # for calculating the variance -- denominator
                       , Tstat_nume = Tstat_nume)
    returnResults
}
CI_set <- NULL
niter <- 1
beta_cand_left <- beta_cand_right <- beta
results_left <- results_right <- c(beta, 999, 999)
while ((sum(results_left[2:3]) > 0 || sum(results_right[2:3]) > 0) && niter <= 150) {
    # print(c(niter=niter, left=beta_cand_left, right=beta_cand_right))
    if (sum(results_left[2:3]) > 0) {
        results_left <- beta_cand_test(beta_cand_left, Y, D, X, A1_ind, weight)
    }
    if (sum(results_right[2:3]) > 0 && beta_cand_left != beta_cand_right) {
        results_right <- beta_cand_test(beta_cand_right, Y, D, X, A1_ind, weight)
    }
    if (length(which(CI_set[, 1] == results_left[1])) == 0) {
        CI_set <- rbind(CI_set, results_left)
    }
    if (length(which(CI_set[, 1] == results_right[1])) == 0) {
        CI_set <- rbind(CI_set, results_right)
    }
    niter <- niter + 1
    beta_cand_left <- beta_cand_left - grid
    beta_cand_right <- beta_cand_right + grid
}
# continue calculating Tstat_prime
Tstat_prime_deno <- var(CI_set[, 8])
Tstat_prime1 <- CI_set[, 6] / Tstat_prime_deno
Tstat_prime2 <- CI_set[, 7] / Tstat_prime_deno
result3 <- qnorm(0.025) <= Tstat_prime1 & Tstat_prime1 <= qnorm(0.975)
result4 <- qnorm(0.025) <= Tstat_prime2 & Tstat_prime2 <= qnorm(0.975)
CI_set <- cbind(CI_set, result3, result4)

CI_set <- CI_set[order(CI_set[, 1]), ]
length1_set <- CI_set[max(which(CI_set[, 2]==1)), 1] - 
                 CI_set[min(which(CI_set[, 2]==1)), 1]
cover1_set <- CI_set[which(CI_set[, 1]==beta), 2]
length2_set <- CI_set[max(which(CI_set[, 3]==1)), 1] -
                 CI_set[min(which(CI_set[, 3]==1)), 1]
cover2_set <- CI_set[which(CI_set[, 1]==beta), 3]
length3_set <- CI_set[max(which(CI_set[, 9]==1)), 1] -
                 CI_set[min(which(CI_set[, 9]==1)), 1]
cover3_set <- CI_set[which(CI_set[, 1]==beta), 9]
length4_set <- CI_set[max(which(CI_set[, 10]==1)), 1] -
                 CI_set[min(which(CI_set[, 10]==1)), 1]
cover4_set <- CI_set[which(CI_set[, 1]==beta), 10]
Tstat1_quantile <- quantile(CI_set[, 4])
Tstat2_quantile <- quantile(CI_set[, 5])
Tstat_prime1_quantile <- quantile(Tstat_prime1)
Tstat_prime2_quantile <- quantile(Tstat_prime2)

### MLIV: AR test in library --------------------------------------------------

Y_A1 <- as.matrix(Y[A1_ind, ])
D_A1 <- as.matrix(D[A1_ind, ])
X_A1 <- as.matrix(X[A1_ind, ])
Y_A2 <- as.matrix(Y[-A1_ind, ])
D_A2 <- as.matrix(D[-A1_ind, ])
X_A2 <- as.matrix(X[-A1_ind, ])
Y_rep <- as.matrix(weight %*% Y_A2)
D_rep <- as.matrix(weight %*% D_A2)
X_rep <- as.matrix(weight %*% X_A2)

reg <- ivmodel(Y_A1, D_A1, D_rep, X_A1, k=1)
beta_hat <- coef(reg)['TSLS', 'Estimate']
sd_hat <- coef(reg)['TSLS', 'Std. Error']
CI_ar <- reg$AR$ci
length <- sum(CI_ar[, 2] - CI_ar[, 1])
cover <- any(CI_ar[, 1] <= beta & beta <= CI_ar[, 2])

### return --------------------------------------------------------------------

returnValue <- c(
    beta_true = beta
    , length1_set = as.numeric(length1_set)
    , length2_set = as.numeric(length2_set)
    , length3_set = as.numeric(length3_set)
    , length4_set = as.numeric(length4_set)
    , length = length
    , cover1_set = as.numeric(cover1_set)
    , cover2_set = as.numeric(cover2_set)
    , cover3_set = as.numeric(cover3_set)
    , cover4_set = as.numeric(cover4_set)
    , cover = cover
    , Tstat1_quantile = Tstat1_quantile
    , Tstat2_quantile = Tstat2_quantile
    , Tstat_prime1_quantile = Tstat_prime1_quantile
    , Tstat_prime2_quantile = Tstat_prime2_quantile
)
returnValue


