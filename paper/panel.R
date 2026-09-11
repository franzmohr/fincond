suppressMessages(library(fincond))
set.seed(1234567); data("eamd")

last_of <- function(z) as.numeric(z)[length(z)]
hac <- function(z, lag) { n <- length(z); z <- z - mean(z); v <- sum(z^2)/n
  for (l in seq_len(lag)) v <- v + 2*(1-l/(lag+1))*sum(z[(l+1):n]*z[1:(n-l)])/n; v }

vintages_of <- function(d, from, n = 1, anchor = "esi") {
  o <- as.numeric(time(window(d$macro, start = from)))
  list(origins = o, macro = d$macro, vintage = lapply(o, function(oo) {
    yr <- floor(oo + 1e-8); q <- round((oo - yr) * 4) + 1
    x <- window(d$financial, end = c(yr, q))
    x <- x[, colSums(!is.na(x)) >= 2, drop = FALSE]
    m <- add_posterior_coefficients(add_priors(create_fcimodel(
      y = window(d$macro, end = c(yr, q)), x = x, p = 4, n = n)))
    fci(m, sign_on = anchor, standardize = FALSE) }))
}
forecasts <- function(target, h, set, cumulative = TRUE) {
  y <- set$macro[, target]; st <- as.numeric(time(y)); ay <- as.numeric(y)
  do.call(rbind, lapply(seq_along(set$origins), function(k) {
    path <- set$vintage[[k]]
    ah <- vapply(seq_len(h), function(j) {
      z <- which(abs(st - (set$origins[k] + j/4)) < 1e-6)
      if (!length(z)) NA_integer_ else z[1] }, integer(1))
    if (anyNA(ah)) return(NULL)
    realised <- if (cumulative) mean(ay[ah]) else ay[ah[h]]
    seen <- window(y, end = end(path))
    lead <- stats::lag(seen, h)
    if (cumulative) lead <- stats::filter(lead, rep(1/h, h), sides = 1)
    d <- as.data.frame(na.omit(stats::ts.intersect(lead = lead, own = seen, fci = path)))
    if (nrow(d) < 12) return(NULL)
    now <- data.frame(own = last_of(seen), fci = last_of(path))
    p_ <- function(f) as.numeric(predict(lm(f, data = d), newdata = now))
    data.frame(origin = set$origins[k], actual = realised,
               mean = mean(d$lead), rw = last_of(seen),
               ar = p_(lead ~ own), fci = p_(lead ~ own + fci)) }))
}

cat("building vintages ...\n")
V <- lapply(eamd, vintages_of, from = c(2006, 1), n = 1)

# Per-origin Clark-West adjusted loss differential, country by country, against
# a fixed benchmark. Averaging across countries first absorbs the common shocks
# into the time-series variance of the average, which is what makes a panel
# test valid when the cross-section is correlated.
adjusted <- function(f, bench) {
  eb <- f$actual - f[[bench]]; ef <- f$actual - f$fci
  eb^2 - (ef^2 - (f[[bench]] - f$fci)^2)
}
panel_test <- function(h, bench) {
  fs <- lapply(V, function(s) forecasts("u", h, s))
  org <- Reduce(intersect, lapply(fs, function(f) f$origin))
  D <- sapply(fs, function(f) adjusted(f[f$origin %in% org, ], bench))
  dbar <- rowMeans(D)
  t_stat <- mean(dbar) / sqrt(hac(dbar, max(h - 1, 0)) / length(dbar))
  # Average pairwise cross-sectional correlation of the loss differentials.
  cc <- cor(D); rho_bar <- mean(cc[upper.tri(cc)])
  c(h = h, n_origins = length(org), t = t_stat, rho_bar = rho_bar)
}

cat("\n=== Panel Clark-West, index vs RANDOM WALK ===\n")
pw <- t(sapply(c(1, 2, 4, 8), panel_test, bench = "rw"))
print(round(pw, 3))
cat("\n=== Panel Clark-West, index vs OWN-LAG REGRESSION ===\n")
pa <- t(sapply(c(1, 2, 4, 8), panel_test, bench = "ar"))
print(round(pa, 3))

cat("\n(one-sided 5% critical value 1.645; rho_bar is the average pairwise\n")
cat(" cross-sectional correlation of the adjusted loss differentials)\n")
saveRDS(list(rw = pw, ar = pa), "panel.rds")
