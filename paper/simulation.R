# Monte Carlo for the mechanism: a near-unit-root target, an estimation sample
# the length of the paper's, and the same four forecasts.
#
#   x_t = phi x_{t-1} + eta_t                      the "index"
#   u_t = rho u_{t-1} + theta x_{t-1} + eps_t      the target
#
# theta > 0 gives the index genuine predictive content at every horizon, so any
# failure below is the benchmark's, not a missing signal.

set.seed(20260911)

fast_ols <- function(X, y) {
  qr.solve(crossprod(X), crossprod(X, y))
}

one_rep <- function(rho, theta, h, tt = 80, phi = 0.8, burn = 200) {
  n <- burn + tt + h + 1
  eta <- rnorm(n); eps <- rnorm(n)
  x <- numeric(n); u <- numeric(n)
  for (i in 2:n) {
    x[i] <- phi * x[i - 1] + eta[i]
    u[i] <- rho * u[i - 1] + theta * x[i - 1] + eps[i]
  }
  x <- x[(burn + 1):n]; u <- u[(burn + 1):n]

  # Estimation sample: origins 1..tt, path target over the next h quarters.
  lead <- vapply(seq_len(tt - h), function(t) mean(u[(t + 1):(t + h)]), numeric(1))
  own <- u[seq_len(tt - h)]
  fci <- x[seq_len(tt - h)]

  origin <- tt
  realised <- mean(u[(origin + 1):(origin + h)])

  X1 <- cbind(1, own)
  X2 <- cbind(1, own, fci)
  b1 <- fast_ols(X1, lead)
  b2 <- fast_ols(X2, lead)

  c(mean = mean(lead),
    rw = u[origin],
    ar = as.numeric(c(1, u[origin]) %*% b1),
    fci = as.numeric(c(1, u[origin], x[origin]) %*% b2),
    actual = realised,
    beta = b1[2])
}

run <- function(rho, theta, h, reps = 4000) {
  m <- t(vapply(seq_len(reps), function(i) one_rep(rho, theta, h), numeric(6)))
  rmse <- function(col) sqrt(mean((m[, "actual"] - m[, col])^2))
  c(rho = rho, h = h,
    beta_hat = mean(m[, "beta"]),
    rmse_mean = rmse("mean"), rmse_rw = rmse("rw"),
    rmse_ar = rmse("ar"), rmse_fci = rmse("fci"),
    ar_over_rw = rmse("ar") / rmse("rw"),
    fci_over_ar = rmse("fci") / rmse("ar"),
    fci_over_best = rmse("fci") / min(rmse("mean"), rmse("rw"), rmse("ar")))
}

grid <- expand.grid(rho = c(0.90, 0.95, 0.97, 0.99, 1.00), h = c(1, 2, 4, 8))
out <- as.data.frame(t(mapply(function(r, h) run(r, h, theta = 0.30),
                              grid$rho, grid$h)))
saveRDS(out, "sim.rds")

cat("theta = 0.30, phi = 0.8, T = 80, 4000 replications\n\n")
cat("=== beta_hat: the estimated own-lag coefficient (truth is rho^h for h=1) ===\n")
b <- reshape(out[, c("rho", "h", "beta_hat")], idvar = "rho", timevar = "h",
             direction = "wide")
print(round(b, 3), row.names = FALSE)

cat("\n=== AR benchmark relative to a random walk ===\n")
a <- reshape(out[, c("rho", "h", "ar_over_rw")], idvar = "rho", timevar = "h",
             direction = "wide")
print(round(a, 3), row.names = FALSE)

cat("\n=== index model relative to the AR benchmark ===\n")
f <- reshape(out[, c("rho", "h", "fci_over_ar")], idvar = "rho", timevar = "h",
             direction = "wide")
print(round(f, 3), row.names = FALSE)

cat("\n=== index model relative to the BEST benchmark ===\n")
g <- reshape(out[, c("rho", "h", "fci_over_best")], idvar = "rho", timevar = "h",
             direction = "wide")
print(round(g, 3), row.names = FALSE)
