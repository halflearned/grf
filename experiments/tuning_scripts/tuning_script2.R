rm(list = ls())

library(grf)
library(stats)



varsigma <- function(x)  {
  1 + 1 / (1 + exp(-20 * (x - (1 / 3))))
}


athey2017estimation1 <- function(n, p) {
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  mu <- 2 * X[,1] - 1
  e <- 1/4 * dbeta(X[, 1], shape1 = 2, shape2 = 4)
  W <- rbinom(n = n, size = 1, prob = e)
  Tau <- 0
  Y <- mu + W * Tau + rnorm(n = n)
  output <- list(X = X, W = W, Y = Y, Tau = Tau, e = e, mu = mu)
  output
}


athey2017estimation2 <- function(n, p) {
  X <- matrix(runif(n * p, min = 0, max = 1), nrow = n, ncol = p)
  mu <- 0
  e <- 0.5
  W <- rbinom(n = n, size = 1, prob = e)
  Tau <- varsigma(X[,1]) * varsigma(X[,2])
  Y <- mu + W * Tau + rnorm(n = n)
  output <- list(X = X, W = W, Y = Y, Tau = Tau, e = e, mu = mu)
  output
}



datasets <- list(
  athey2017estimation1 = athey2017estimation1,
  athey2017estimation2 = athey2017estimation2
)


filename <- sprintf("result_%d.csv", sample(100000000, 1))

for (i in seq(1000)) {
  tryCatch({

    data_idx <- sample.int(length(datasets), 1)
    data_function <- datasets[[data_idx]]
    data_name <- names(datasets)[data_idx]

    n <- sample(c(500, 1000, 2000), 1)
    p <- sample(seq(5, 10), 1)
    df <- data_function(n, p)

    # Forests, with and without "oracle" centering
    tau.forest.kriging <- causal_forest(df$X, df$Y, df$W, tune.parameters = TRUE)
    tau.forest.notune <- causal_forest(df$X, df$Y, df$W, tune.parameters = FALSE)
    tau.forest.kriging.centered <- causal_forest(df$X, df$Y, df$W, Y.hat = df$mu, W.hat = df$e, tune.parameters = TRUE)
    tau.forest.notune.centered <- causal_forest(df$X, df$Y, df$W, Y.hat = df$mu, W.hat = df$e, tune.parameters = FALSE)

    # Predictions
    tauhat.kriging <- predict(tau.forest.kriging)$predictions
    tauhat.notune <- predict(tau.forest.notune)$predictions
    tauhat.kriging.centered <- predict(tau.forest.kriging.centered)$predictions
    tauhat.notune.centered <- predict(tau.forest.notune.centered)$predictions

    # 'Raw' oracle MSE
    mse.kriging <- mean((tauhat.kriging - df$Tau)^2)
    mse.notune <- mean((tauhat.notune - df$Tau)^2)
    mse.kriging.centered <- mean((tauhat.kriging.centered - df$Tau)^2)
    mse.notune.centered <- mean((tauhat.notune.centered - df$Tau)^2)

    # 'Weighted' oracle MSE
    var_w <- df$e * (1 - df$e)
    wmse.kriging <- mean(var_w * (tauhat.kriging - df$Tau)^2)
    wmse.notune <- mean(var_w * (tauhat.notune - df$Tau)^2)
    wmse.kriging.centered <- mean(var_w * (tauhat.kriging.centered - df$Tau)^2)
    wmse.notune.centered <- mean(var_w * (tauhat.notune.centered - df$Tau)^2)

    # Tally up
    cfg <- c(data = data_name, n = n, p = p)
    v1 <- data.frame(rbind(c(cfg, tau.forest.kriging$tunable.params, method = "kriging", stat = "mse", value = mse.kriging)))
    v2 <- data.frame(rbind(c(cfg, tau.forest.notune$tunable.params, method = "notune", stat = "mse", value = mse.notune)))
    v3 <- data.frame(rbind(c(cfg, tau.forest.kriging.centered$tunable.params, method = "kriging:centered", stat = "mse", value = mse.kriging.centered)))
    v4 <- data.frame(rbind(c(cfg, tau.forest.notune.centered$tunable.params, method = "notune:centered", stat = "mse", value = mse.notune.centered)))
    v5 <- data.frame(rbind(c(cfg, tau.forest.kriging$tunable.params, method = "kriging", stat = "wmse", value = wmse.kriging)))
    v6 <- data.frame(rbind(c(cfg, tau.forest.notune$tunable.params, method = "notune", stat = "wmse", value = wmse.notune)))
    v7 <- data.frame(rbind(c(cfg, tau.forest.kriging.centered$tunable.params, method = "kriging:centered", stat = "wmse", value = wmse.kriging.centered)))
    v8 <- data.frame(rbind(c(cfg, tau.forest.notune.centered$tunable.params, method = "notune:centered", stat = "wmse", value = wmse.notune.centered)))

    result <- rbind(v1, v2, v3, v4, v5, v6, v7, v8)
    write.table(file = filename, x = result, append = TRUE, quote = FALSE, row.names = FALSE, col.names = i == 1, sep = ",")
    print(result)
  }, error = function(e) print("Error"))
}




