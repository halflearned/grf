rm(list = ls())

library(grf)




athey2017estimation1 <- function(n, p = 10) {
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  noise <- stats::rnorm(n = n)
  main_effect <- rep(0, times = n)
  e <- rep(0.5, n)
  W <- matrix(stats::rbinom(n = n, size = 1, prob = e), nrow = n, ncol = 1)
  zeta1 <- 1 + 1 / (1 + exp(-20 * (X[, 1] - (1 / 3))))
  zeta2 <- 1 + 1 / (1 + exp(-20 * (X[, 2] - (1 / 3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + W * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau, e = e)
  output
}


athey2017estimation2 <- function(n, p = 10) {
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  noise <- stats::rnorm(n = n)
  main_effect <- rep(0, times = n)
  e <- rep(0.5, n)
  W <- matrix(stats::rbinom(n = n, size = 1, prob = e), nrow = n, ncol = 1)
  zeta1 <- 2 / (1 + exp(-12 * (X[, 1] - (1 / 2))))
  zeta2 <- 2 / (1 + exp(-12 * (X[, 2] - (1 / 2))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + W * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, W = W, Y = Y, Tau = Tau, e = e)
  output
}

xlearners1 <- function(n, p = 10) {
  e <- rep(0.05, n)
  sigma <- diag(rep(1, p)) + matrix(1, p, p)
  means <- rep(0, p)
  X <- MASS::mvrnorm(n = n, mu = means, Sigma = sigma)
  beta <- runif(p, -5, 5)
  mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
  mu_1 <- mu_0 + 8 * (X[, 2] > 0.1) + stats::rnorm(n = n)
  W <- matrix(stats::rbinom(n = n, size = 1, prob = e), nrow = n, ncol = 1)
  y <- matrix(ifelse(W == 0, mu_0, mu_1), nrow = n, ncol = 1)
  tau <- matrix(8 * (X[, 2] > 0.1), nrow = n, ncol = 1)
  output <- list(X = X, Y = y, W = W, Tau = tau, e = e)
  output
}

xlearners4 <- function(n, p = 10) {
  e <- rep(0.5, n)
  sigma <- diag(rep(1, p)) + matrix(1, p, p)
  means <- rep(0, p)
  X <- MASS::mvrnorm(n = n, mu = means, Sigma = sigma)
  Beta <- runif(p)
  mu_0 <- X %*% Beta
  mu_1 <- mu_0
  W <- matrix(stats::rbinom(n = n, size = 1, prob = e), nrow = n, ncol = 1)
  y <- matrix(ifelse(W == 0, mu_0, mu_1), nrow = n, ncol = 1)
  tau <- matrix(mu_1 - mu_0, nrow = n, ncol = 1)
  output <- list(X = X, Y = y, W = W, Tau = tau, e = e)
  output
}

xlearners6 <- function(n) {
  X <- matrix(runif(n = n * 20, min = 0, max = 1), nrow = n, ncol = 20)
  e <- 0.25 * (1 + stats::dbeta(X[, 1], shape1 = 2, shape2 = 4))
  mu_0 <- 2 * X[, 1] - 1
  mu_1 <- mu_0
  W <- matrix(stats::rbinom(n = n, size = 1, prob = e), nrow = n, ncol = 1)
  y <- matrix(ifelse(W == 0, mu_0, mu_1), nrow = n, ncol = 1)
  tau <- matrix(mu_1 - mu_0, nrow = n, ncol = 1)
  output <- list(X = X, Y = y, W = W, Tau = tau, e = e)
  output
}


datasets <- list(
  athey2017estimation1 = athey2017estimation1,
  athey2017estimation2 = athey2017estimation2,
  xlearners1 = xlearners1,
  xlearners4 = xlearners4,
  xlearners6 = xlearners6
)


filename <- sprintf("result_%d.csv", sample(100000000, 1))

for (i in seq(1000)) {
  print(i)
  data_idx <- sample.int(length(datasets), 1)
  data_function <- datasets[[data_idx]]
  data_name <- names(datasets)[data_idx]
  n <- sample(c(500, 1000, 2000), 1)

  df <- data_function(n)

  X <- df$X
  Y <- as.numeric(df$Y)
  W <- as.numeric(df$W)
  var_w <- df$e * (1 - df$e)
  Tau <- as.numeric(df$Tau)

  tau.forest.kriging <- causal_forest(X, Y, W, tune.parameters = TRUE)
  tau.forest.notune <- causal_forest(X, Y, W, tune.parameters = FALSE)

  tauhat.kriging <- predict(tau.forest.kriging)$predictions
  tauhat.notune <- predict(tau.forest.notune)$predictions

  mse.kriging <- mean((tauhat.kriging - Tau)^2)
  mse.notune <- mean((tauhat.notune - Tau)^2)
  wmse.kriging <- mean(var_w * (tauhat.kriging - Tau)^2)
  wmse.notune <- mean(var_w * (tauhat.notune - Tau)^2)

  cfg <- c(data = data_name, n = n, p = dim(X)[2])
  v1 <- data.frame(rbind(c(cfg, tau.forest.kriging$tunable.params, method = "kriging", stat = "mse", value = mse.kriging)))
  v2 <- data.frame(rbind(c(cfg, tau.forest.notune$tunable.params, method = "notune", stat = "mse", value = mse.notune)))
  v3 <- data.frame(rbind(c(cfg, tau.forest.kriging$tunable.params, method = "kriging", stat = "wmse", value = wmse.kriging)))
  v4 <- data.frame(rbind(c(cfg, tau.forest.notune$tunable.params, method = "notune", stat = "wmse", value = wmse.notune)))

  result <- rbind(v1, v2, v3, v4)

  write.table(file = filename, x = result, append = TRUE, quote = FALSE, row.names = FALSE, col.names = i == 1, sep = ",")
  print(result)
}
