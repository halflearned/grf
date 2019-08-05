rm(list=ls())
library(grf)

num_sims = 1000
filename = paste0(paste0("results",
                   Sys.getenv('SLURM_JOB_ID'),
                   Sys.getenv('SLURM_LOCALID'),
                   Sys.getenv('SLURM_JOB_NAME'),
                   as.integer(Sys.time(), collapse="_"),
                   ".csv"))

for (s in seq(num_sims)) {

  print(paste0("Simulation ", s))

  # Generate data.
  dgp = sample(c("simple", "aw1", "aw2", "amle1", "amle2"), 1)
  n = sample(c(200, 1000, 5000), 1)
  p = sample(c(5, 10, 20), 1)

  if (dgp == "simple") {
    X = matrix(rnorm(n*p), n, p)
    W = rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))
    TAU = pmax(X[,1], 0)
    Y = X[,2] + pmin(X[,3], 0) + TAU * W + rnorm(n)
  } else if (dgp == "aw1") {
    X = matrix(rnorm(n*p), n, p)
    e = 1/4*(1 + dbeta(X[,1], 2, 4))
    W = rbinom(n, 1, e)
    m <- 2*X[,1] - 1
    TAU <- rep(0, n)
    Y <- m + rnorm(n)
  } else if (dgp == "aw2") {
    X = matrix(rnorm(n*p), n, p)
    e = 0.5
    W = rbinom(n, 1, e)
    zeta1 <- 1 + 1/(1 + exp(-20 * (X[, 1] - (1/3))))
    zeta2 <- 1 + 1/(1 + exp(-20 * (X[, 2] - (1/3))))
    m <- 2*X[,1] - 1
    TAU <- zeta1 * zeta2
    Y <- m + W*TAU + rnorm(n)
  } else if (dgp == "amle1") {
    X = matrix(rnorm(n*p), n, p)
    zeta = apply(X, 1, sum) / sqrt(p)
    eta = sign(zeta) * zeta^2
    alpha = pmax(0.05, pmin(0.95, 1/(1 + exp(-eta))))
    W = rbeta(n, alpha, 1-alpha)
    mu = eta + 0.2*(alpha - 0.5)
    TAU = -0.2
    Y = mu + W*TAU + rnorm(n)
  } else if (dgp == "amle2") {
    X = matrix(rnorm(n*p), n, p)
    eta = apply(X, 1, prod) * 2^(p-1)
    mu = sign(eta) * sqrt(abs(eta))
    lambda = 0.1*sign(mu) + mu
    W = rnorm(n, lambda, abs(lambda))
    TAU = pmax(X[,1] + X[,2], 0)
    Y = W*TAU + rnorm(n)
  } else if (dgp == "ai1") {
     X = matrix(rnorm(n*p), n, p)
     W <- rbinom(n = n, size = 1, prob = 0.5)
     nu_x <- 0.5 * X[, 1] + X[, 2]
     Tau <- 0.25 * X[, 1]
     Y <- nu_x +  Tau * W + rnorm(n = n, sd = 0.1)
   } else if (dgp == "ai2") {
     X = matrix(rnorm(n*p), n, p)
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
  }

  # Estimate forests
  cf.tuned = causal_forest(X, Y, W, tune.parameters=T)
  cf.not = causal_forest(X, Y, W, tune.parameters=F)

  # Estimate treatment effect on oob samples
  tau.hat.tuned = predict(cf.tuned)$predictions
  tau.hat.not = predict(cf.not)$predictions

  # Compute mse
  mse.tuned = mean((TAU - tau.hat.tuned)^2)
  mse.not = mean((TAU - tau.hat.not)^2)

  # Save results
  sim.params = c(dgp=dgp, n=n, p=p)
  res1 = rbind(c(sim.params, cf.tuned$tuning.output$status, mse.oob=mse.tuned))
  res2 = rbind(c(sim.params, "notune", mse.oob=mse.not))
  write.table(res1, file=filename, col.names=s == 0, row.names=F, append=T)
  write.table(res2, file=filename, col.names=s == 0, row.names=F, append=T)

}
