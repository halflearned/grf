#!/usr/bin/env Rscript

rm(list=ls())
library(grf)
library(stats)

# These environment variables exist in the Sherlock cluster
filename = paste0(paste0("results",
       Sys.getenv('SLURM_JOB_ID'),
       Sys.getenv('SLURM_LOCALID'),
       Sys.getenv('SLURM_JOB_NAME'),
       as.integer(Sys.time()), collapse="_"), ".csv")

num_sims = 500
for (s in seq(num_sims)) {
    print(paste0("Simulation ", s))

    # Generate data.
    dgp = sample(c("simple", "aw1", "aw2", "aw3", "ai1", "ai2", "kunzel"), 1)
    n = 10 #sample(c(250, 1000, 5000), 1)
    p = sample(c(10, 20), 1)
    tune = sample(c(TRUE, FALSE), 1)
    nft = sample(c(200, 1000), 1)

    # Create data
    if (dgp == "simple") {
      X = matrix(rnorm(n*p), n, p)
      W = rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))
      Tau = pmax(X[,1], 0)
      Y = X[,2] + pmin(X[,3], 0) + Tau * W + rnorm(n)
    } else if (dgp == "aw1") {
      X <- matrix(runif(n * p), n, p)
      W <- rbinom(n, 1, 0.5)
      zeta1 <- 1 + 1/(1 + exp(-20 * (X[, 1] - (1/3))))
      zeta2 <- 1 + 1/(1 + exp(-20 * (X[, 2] - (1/3))))
      Tau <- zeta1 * zeta2
      Y <-  W * Tau + rnorm(n)
    } else if (dgp == "aw2") {
      X <- matrix(runif(n * p), n, p)
      W <- rbinom(n, 1, 0.5)
      zeta1 <- 2/(1 + exp(-12 * (X[, 1] - (1/2))))
      zeta2 <- 2/(1 + exp(-12 * (X[, 2] - (1/2))))
      Tau <- matrix(zeta1 * zeta2, n, 1)
      Y <- W * Tau + rnorm(n = n)
    } else if (dgp == "aw3") {
      X <- matrix(runif(n * p, min = 0, max = 1), n, p)
      Tau <- 0
      treatment_propensity <- (1/4) * (1 + dbeta(X[, 1], 2, 4))
      W <- rbinom(n = n, size = 1, prob = treatment_propensity)
      Y <- 2 * X[, 1] - 1 + rnorm(n = n)
    } else if (dgp == "ai1") {
      X <- matrix(rnorm(n, p), n, p)
      W <- rbinom(n = n, size = 1, prob = 0.5)
      nu_x <- 0.5 * X[, 1] + X[, 2]
      Tau <- 0.25 * X[, 1]
      Y <- nu_x +  Tau * W + rnorm(n = n, sd = 0.1)
    } else if (dgp == "ai2") {
      X <- matrix(rnorm(n, p), n, p)
      W <- rbinom(n = n, size = 1, prob = 0.5)
      nu_x <- 0.5 * X[, 1] + 0.5 * X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6]
      Tau <- 0.5 * ((X[, 1] > 0) * X[, 1] + (X[, 2] > 0) * X[, 2])
      Y <- nu_x +  Tau * W + rnorm(n = n, sd = 0.1)
    } else if (dgp == "kunzel") {
      X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = toeplitz(0.5^seq(0, p-1)))
      beta <- runif(p, -5, 5)
      mu_0 <- X %*% beta + 5 * (X[, 1] > 0.5) + rnorm(n = n)
      mu_1 <- mu_0 + 8 * (X[, 2] > 0.1) + rnorm(n = n)
      W <- rbinom(n = n, size = 1, prob = 0.01)
      Y <- ifelse(W == 0, mu_0, mu_1)
      Tau <- 8 * (X[, 2] > 0.1)
    } else  {
      stop(paste0("Unknown dgp: ", dgp))
    }

    # Estimate the forest
    cf = causal_forest(X, Y, W, tune.parameters=tune)

    # Estimate treatment effect on oob samples
    tau.hat.oob = predict(cf)$predictions

    # Compute mse
    mse.oob = mean((Tau - tau.hat.oob)^2)

    status = ifelse(tune, cf$tuning.output$status, "notune")

    # Save results
    res = rbind(c(dgp=dgp, n=n, p=p, tune=tune,
      status=status, num.fit.trees=nft, mse.oob=mse.oob))
    print(res)
    write.table(res, file=filename, col.names=s == 0, row.names=F, append=T)
}
