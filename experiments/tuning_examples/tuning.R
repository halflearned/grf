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
  X = matrix(rnorm(n*p), n, p)

  if (dgp == "simple") {
    W = rbinom(n, 1, 0.4 + 0.2 * (X[,1] > 0))
    TAU = pmax(X[,1], 0)
    Y = X[,2] + pmin(X[,3], 0) + TAU * W + rnorm(n)
  } else if (dgp == "aw1") {
    e = 1/4*(1 + dbeta(X[,1], 2, 4))
    W = rbinom(n, 1, e)
    m <- 2*X[,1] - 1
    TAU <- rep(0, n)
    Y <- m + rnorm(n)
  } else if (dgp == "aw2") {
    e = 0.5
    W = rbinom(n, 1, e)
    zeta1 <- 1 + 1/(1 + exp(-20 * (X[, 1] - (1/3))))
    zeta2 <- 1 + 1/(1 + exp(-20 * (X[, 2] - (1/3))))
    m <- 2*X[,1] - 1
    TAU <- zeta1 * zeta2
    Y <- m + W*TAU + rnorm(n)
  } else if (dgp == "amle1") {
    zeta = apply(X, 1, sum) / sqrt(p)
    eta = sign(zeta) * zeta^2
    alpha = pmax(0.05, pmin(0.95, 1/(1 + exp(-eta))))
    W = rbeta(n, alpha, 1-alpha)
    mu = eta + 0.2*(alpha - 0.5)
    TAU = -0.2
    Y = mu + W*TAU + rnorm(n)
  } else if (dgp == "amle2") {
    eta = apply(X, 1, prod) * 2^(p-1)
    mu = sign(eta) * sqrt(abs(eta))
    lambda = 0.1*sign(mu) + mu
    W = rnorm(n, lambda, abs(lambda))
    TAU = pmax(X[,1] + X[,2], 0)
    Y = W*TAU + rnorm(n)
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
