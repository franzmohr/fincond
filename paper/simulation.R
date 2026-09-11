# Consolidated Monte Carlo for the manuscript.
#
#   x_t = phi x_{t-1} + eta_t                      the "index"
#   u_t = rho u_{t-1} + theta x_{t-1} + eps_t      the target
#
# phi = 0.4 keeps the first-order autocorrelation of u within about 0.02 of
# rho; with a more persistent predictor, x feeds its own autocorrelation into u
# and rho stops describing the target's persistence. theta = 0.3 gives the
# predictor genuine content at every horizon, so anything the benchmark does
# wrong below is the benchmark's doing.

set.seed(20260911)

PHI <- 0.4; THETA <- 0.3; TT <- 80
RHOS <- c(0.80, 0.90, 0.95, 0.97, 0.99)
HS <- c(1, 2, 4, 8)

hac <- function(z, lag) {
  n <- length(z); z <- z - mean(z); v <- sum(z^2) / n
  for (l in seq_len(lag)) v <- v + 2*(1-l/(lag+1))*sum(z[(l+1):n]*z[1:(n-l)])/n
  v
}
draw <- function(rho, n, burn = 200) {
  N <- n + burn
  eta <- rnorm(N); eps <- rnorm(N)
  x <- stats::filter(eta, PHI, method = "recursive")
  u <- stats::filter(eps + THETA * c(0, x[-N]), rho, method = "recursive")
  list(x = as.numeric(x)[(burn + 1):N], u = as.numeric(u)[(burn + 1):N])
}
# Trailing mean of the h-step lead, computed once for the whole series.
lead_of <- function(u, h) {
  cs <- cumsum(c(0, u))
  n <- length(u)
  out <- rep(NA_real_, n)
  idx <- seq_len(n - h)
  out[idx] <- (cs[idx + h + 1] - cs[idx + 1]) / h
  out
}
# OLS on an intercept plus 1 or 2 regressors, returning the fitted value at
# the supplied point. Solved directly; the systems are 2x2 and 3x3.
fit_predict <- function(y, Z, z_new) {
  X <- cbind(1, Z)
  b <- qr.solve(crossprod(X), crossprod(X, y))
  as.numeric(c(1, z_new) %*% b)
}

## 1. Population quantities --------------------------------------------------
cat("population quantities ...\n")
pop <- do.call(rbind, lapply(RHOS, function(rho) {
  d <- draw(rho, 2e5)
  ac <- as.numeric(acf(d$u, lag.max = 1, plot = FALSE)$acf[2])
  b <- vapply(HS, function(h) {
    L <- lead_of(d$u, h); k <- which(!is.na(L))
    as.numeric(coef(lm(L[k] ~ d$u[k]))[2]) }, numeric(1))
  data.frame(rho = rho, eff_rho = ac, t(setNames(round(b, 3), paste0("pop", HS))))
}))

## 2. Accuracy ---------------------------------------------------------------
cat("accuracy grid ...\n")
one_rep <- function(rho, h) {
  d <- draw(rho, TT + h + 1); u <- d$u; x <- d$x
  L <- lead_of(u, h); k <- seq_len(TT - h)
  y <- L[k]; own <- u[k]; fci <- x[k]
  c(actual = mean(u[(TT + 1):(TT + h)]), mean = mean(y), rw = u[TT],
    ar = fit_predict(y, own, u[TT]),
    fci = fit_predict(y, cbind(own, fci), c(u[TT], x[TT])),
    beta = as.numeric(qr.solve(crossprod(cbind(1, own)),
                               crossprod(cbind(1, own), y))[2]))
}
acc <- do.call(rbind, lapply(RHOS, function(rho) do.call(rbind, lapply(HS, function(h) {
  m <- t(vapply(seq_len(6000), function(i) one_rep(rho, h), numeric(6)))
  r <- function(col) sqrt(mean((m[, "actual"] - m[, col])^2))
  data.frame(rho = rho, h = h, beta_hat = mean(m[, "beta"]),
             ar_over_rw = r("ar") / r("rw"),
             fci_over_ar = r("fci") / r("ar"),
             fci_over_best = r("fci") / min(r("mean"), r("rw"), r("ar")))
}))))

## 3. Inference --------------------------------------------------------------
cat("inference grid ...\n")
one_path <- function(rho, h, n_org = 79) {
  d <- draw(rho, TT + n_org + h + 1); u <- d$u; x <- d$x
  L <- lead_of(u, h)
  out <- matrix(NA_real_, n_org, 4,
                dimnames = list(NULL, c("actual", "rw", "ar", "fci")))
  for (j in seq_len(n_org)) {
    o <- TT + j - 1
    k <- seq_len(o - h)
    y <- L[k]; own <- u[k]; fci <- x[k]
    out[j, ] <- c(mean(u[(o + 1):(o + h)]), u[o],
                  fit_predict(y, own, u[o]),
                  fit_predict(y, cbind(own, fci), c(u[o], x[o])))
  }
  out
}
cw <- function(m, bench, h) {
  eb <- m[, "actual"] - m[, bench]; ef <- m[, "actual"] - m[, "fci"]
  adj <- eb^2 - (ef^2 - (m[, bench] - m[, "fci"])^2)
  mean(adj) / sqrt(hac(adj, max(h - 1, 0)) / length(adj))
}
inf <- do.call(rbind, lapply(RHOS, function(rho) do.call(rbind, lapply(HS, function(h) {
  z <- t(vapply(seq_len(500), function(i) {
    m <- one_path(rho, h); c(cw(m, "ar", h), cw(m, "rw", h)) }, numeric(2)))
  data.frame(rho = rho, h = h,
             rej_ar = mean(z[, 1] > 1.645), rej_rw = mean(z[, 2] > 1.645))
}))))

saveRDS(list(pop = pop, acc = acc, inf = inf), "sim4.rds")

wide <- function(d, col) {
  w <- reshape(d[, c("rho", "h", col)], idvar = "rho", timevar = "h",
               direction = "wide")
  colnames(w) <- c("rho", paste0("h=", HS)); w
}
cat("\nphi =", PHI, " theta =", THETA, " T =", TT, "\n")
cat("\n=== effective persistence and population projection coefficient ===\n")
print(pop, row.names = FALSE)
for (v in c("beta_hat", "ar_over_rw", "fci_over_ar", "fci_over_best")) {
  cat("\n===", v, "===\n"); print(round(wide(acc, v), 3), row.names = FALSE)
}
cat("\n=== CW rejection rate vs own-lag regression ===\n")
print(round(wide(inf, "rej_ar"), 3), row.names = FALSE)
cat("\n=== CW rejection rate vs random walk ===\n")
print(round(wide(inf, "rej_rw"), 3), row.names = FALSE)
