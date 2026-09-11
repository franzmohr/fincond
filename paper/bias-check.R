# Does the standard small-sample bias of the least-squares AR(1) coefficient
# account for what the simulation shows? Marriott and Pope (1954) and Kendall
# (1954) give, for an AR(1) with an intercept estimated on T observations,
#
#     E[rho_hat] - rho  ~  -(1 + 3 rho) / T
#
# This isolates that channel: theta = 0, so the target is a pure AR(1) and the
# only thing acting on rho_hat is the estimator's own bias.

set.seed(11)

sim_rho <- function(rho, tt, reps = 20000, burn = 200) {
  v <- vapply(seq_len(reps), function(i) {
    n <- burn + tt + 1
    u <- numeric(n)
    e <- rnorm(n)
    for (j in 2:n) u[j] <- rho * u[j - 1] + e[j]
    u <- u[(burn + 1):n]
    coef(lm(u[-1] ~ u[-length(u)]))[2]
  }, numeric(1))
  mean(v)
}

cat("T = 80: simulated mean rho_hat vs the Marriott-Pope approximation\n")
cat(" rho    simulated   predicted   implied h=8 attenuation\n")
for (rho in c(0.90, 0.95, 0.97, 0.99, 1.00)) {
  s <- sim_rho(rho, 80)
  pred <- rho - (1 + 3 * rho) / 80
  # forecast attenuation at h = 8: rho^8 - rho_hat^8, relative to rho^8
  att <- (rho^8 - s^8) / rho^8
  cat(sprintf("%5.2f   %9.3f   %9.3f   %20.1f%%\n", rho, s, pred, 100 * att))
}

cat("\nBias of imposing a unit root instead, for comparison:\n")
cat(" rho   (1 - rho)   vs estimation bias (1+3rho)/80\n")
for (rho in c(0.90, 0.95, 0.97, 0.99, 1.00)) {
  cat(sprintf("%5.2f   %8.3f   %25.3f\n", rho, 1 - rho, (1 + 3 * rho) / 80))
}
