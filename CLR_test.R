### to run the CLR method (traditional and ML), IV ranges from weak to strong
library(drf)
library(ivmodel)

n <- 1000  # {1000, 3000, 5000}
a <- 0  # {0 ~ 1}

### simulation ----------------------------------------------------------------

sim_func <- function(n, a, B=300, grid_val=0.01) {
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
    
    ### Stage 2: CLR --------------------------------------------------------------
    
    Y_A1 <- as.matrix(Y[A1_ind, ])
    D_A1 <- as.matrix(D[A1_ind, ])
    X_A1 <- as.matrix(X[A1_ind, ])
    Y_A2 <- as.matrix(Y[-A1_ind, ])
    D_A2 <- as.matrix(D[-A1_ind, ])
    X_A2 <- as.matrix(X[-A1_ind, ])
    Y_rep <- as.matrix(weight %*% Y_A2)
    D_rep <- as.matrix(weight %*% D_A2)
    X_rep <- as.matrix(weight %*% X_A2)
    
    quantile_generation <- function(beta_cand, beta, B, Y_rep, D_rep, X_rep) {
        D_rep_resid <- resid(lm(D_rep ~ X_rep))
        Tstat <- resid(lm(Y_rep-D_rep*beta_cand ~ X_rep))
        reg <- lm(D_rep_resid ~ Tstat - 1)
        tau <- coef(reg)
        RemBeta <- resid(reg)
        beta_sample <- NULL
        for (j in 1:B) {
            T_sd <- sqrt(mean(Tstat^2))
            T_sample <- rnorm(length(Tstat), mean=Tstat, sd=T_sd)
            G_sample <- tau * T_sample + RemBeta
            beta_m <- t(G_sample) %*% (T_sample + beta_cand * G_sample) / crossprod(G_sample)
            beta_sample <- c(beta_sample, beta_m)
        }
        quantile_val <- quantile(beta_sample, c(0.025, 0.975))
        cover <- as.numeric(quantile_val[1] <= beta & beta <= quantile_val[2])
        returnValue <- c(
            beta = beta_cand
            , cover = cover
            , quantile_val_lower = as.numeric(quantile_val[1])
            , quantile_val_upper = as.numeric(quantile_val[2])
        )
        returnValue
    }
    
    ### MLIV: CLR -----------------------------------------------------------------
    
    MLIV_quantile <- function(beta_cand, beta, Y_A1, D_A1, D_rep, X_A1) {
        reg <- ivmodel(Y_A1, D_A1, D_rep, X_A1)
        test_stat <- CLR(reg, beta0=beta_cand, alpha=0.05)$test.stat
        test_pvalue <- CLR(reg, beta0=beta_cand, alpha=0.05)$p.value
        returnValue <- c(
            beta = beta_cand
            , stat = test_stat
            , pvalue = test_pvalue
        )
        returnValue
    }
    
    ### start searching beta ------------------------------------------------------
    
    beta_cand_set <- seq(-1, 3, grid_val)
    CI_set <- NULL
    for (i in 1:length(beta_cand_set)) {
        # print(i / length(beta_cand_set))
        result_test <- quantile_generation(beta_cand_set[i], beta, B, Y_rep, D_rep, X_rep)
        result_mliv <- MLIV_quantile(beta_cand_set[i], beta, Y_A1, D_A1, D_rep, X_A1)
        result <- c(beta = beta_cand_set[i]
                    , test = result_test[-1]
                    , mliv = result_mliv[-1])
        CI_set <- rbind(CI_set, result)
    }
    CI_test <- as.matrix(CI_set[CI_set[,2]==1, ], ncol=ncol(CI_set))
    if (nrow(CI_test)==0) {
        length_test <- NA 
        cover_test <- 0
    } else {
        length_test <- max(CI_test[, 1]) - min(CI_test[, 1])
        cover_test <- as.numeric(beta %in% CI_test[, 1])
    }
    CI_mliv <- CI_set[CI_set[,6]>0.05, ]
    length_mliv <- max(CI_mliv[, 1]) - min(CI_mliv[, 1])
    cover_mliv <- as.numeric(beta %in% CI_mliv[, 1])
    
    ### return --------------------------------------------------------------------
    
    returnValue <- c(
        beta_true = beta
        , length_test = length_test
        , length_mliv = length_mliv
        , cover_test = cover_test
        , cover_mliv = cover_mliv
        , quantile_lower = quantile(CI_set[, 3])
        , quantile_upper = quantile(CI_set[, 4])
    )
    returnValue
}

nsim <- 30
total <- NULL
for (i in 1:nsim) {
    print(c(a, n, i))
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

### Stage 2: CLR --------------------------------------------------------------

Y_A1 <- as.matrix(Y[A1_ind, ])
D_A1 <- as.matrix(D[A1_ind, ])
X_A1 <- as.matrix(X[A1_ind, ])
Y_A2 <- as.matrix(Y[-A1_ind, ])
D_A2 <- as.matrix(D[-A1_ind, ])
X_A2 <- as.matrix(X[-A1_ind, ])
Y_rep <- as.matrix(weight %*% Y_A2)
D_rep <- as.matrix(weight %*% D_A2)
X_rep <- as.matrix(weight %*% X_A2)

# beta_hat
init_reg <- lm(Y_rep ~ D_rep + X_rep)
beta_init <- coef(init_reg)[2]
delta_hat <- D_A1 - D_rep
eps_hat <- resid(lm(Y_A1-D_A1*beta_init ~ X_A1))

diag_val = rep(NA, n-n1)
for (j in 1:(n-n1)) {
    diag_val[j] <- sum(resid(lm(weight[, j] ~ X_rep))^2)
}
D_rep_resid <- resid(lm(D_rep ~ X_rep))
deno <- sum(D_rep_resid^2)
bias_hat <- sum(diag_val*delta_hat*eps_hat) / deno




quantile_generation <- function(beta_cand, beta, B, Y_rep, D_rep, X_rep) {
    D_rep_resid <- resid(lm(D_rep ~ X_rep))
    Tstat <- resid(lm(Y_rep-D_rep*beta_cand ~ X_rep))
    reg <- lm(D_rep_resid ~ Tstat - 1)
    tau <- coef(reg)
    RemBeta <- resid(reg)
    beta_sample <- NULL
    for (j in 1:B) {
        T_sd <- sqrt(mean(Tstat^2))
        T_sample <- rnorm(length(Tstat), mean=Tstat, sd=T_sd)
        G_sample <- tau * T_sample + RemBeta
        beta_m <- t(G_sample) %*% (T_sample + beta_cand * G_sample) / crossprod(G_sample)
        beta_sample <- c(beta_sample, beta_m)
    }
    quantile_val <- quantile(beta_sample, c(0.025, 0.975))
    cover <- as.numeric(quantile_val[1] <= beta & beta <= quantile_val[2])
    returnValue <- c(
        beta = beta_cand
        , cover = cover
        , quantile_val_lower = as.numeric(quantile_val[1])
        , quantile_val_upper = as.numeric(quantile_val[2])
    )
    returnValue
}

### MLIV: CLR -----------------------------------------------------------------

MLIV_quantile <- function(beta_cand, beta, Y_A1, D_A1, D_rep, X_A1) {
    reg <- ivmodel(Y_A1, D_A1, D_rep, X_A1)
    test_stat <- CLR(reg, beta0=beta_cand, alpha=0.05)$test.stat
    test_pvalue <- CLR(reg, beta0=beta_cand, alpha=0.05)$p.value
    returnValue <- c(
        beta = beta_cand
        , stat = test_stat
        , pvalue = test_pvalue
    )
    returnValue
}

### start searching beta ------------------------------------------------------

beta_cand_set <- seq(-1, 3, grid_val)
CI_set <- NULL
for (i in 1:length(beta_cand_set)) {
    # print(i / length(beta_cand_set))
    result_test <- quantile_generation(beta_cand_set[i], beta, B, Y_rep, D_rep, X_rep)
    result_mliv <- MLIV_quantile(beta_cand_set[i], beta, Y_A1, D_A1, D_rep, X_A1)
    result <- c(beta = beta_cand_set[i]
                , test = result_test[-1]
                , mliv = result_mliv[-1])
    CI_set <- rbind(CI_set, result)
}
CI_test <- as.matrix(CI_set[CI_set[,2]==1, ], ncol=ncol(CI_set))
if (nrow(CI_test)==0) {
    length_test <- NA 
    cover_test <- 0
} else {
    length_test <- max(CI_test[, 1]) - min(CI_test[, 1])
    cover_test <- as.numeric(beta %in% CI_test[, 1])
}
CI_mliv <- CI_set[CI_set[,6]>0.05, ]
length_mliv <- max(CI_mliv[, 1]) - min(CI_mliv[, 1])
cover_mliv <- as.numeric(beta %in% CI_mliv[, 1])

### return --------------------------------------------------------------------

returnValue <- c(
    beta_true = beta
    , length_test = length_test
    , length_mliv = length_mliv
    , cover_test = cover_test
    , cover_mliv = cover_mliv
    , quantile_lower = quantile(CI_set[, 3])
    , quantile_upper = quantile(CI_set[, 4])
)
returnValue