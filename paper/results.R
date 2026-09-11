# Replication for "Benchmark choice and the apparent long-horizon predictive
# content of financial conditions indices".
#
# Run from the package root with fincond installed. Writes bench3.rds, which
# paper/tables.R turns into the numbers in the manuscript. Takes about a minute.

library("fincond")
set.seed(1234567); data("eamd"); data("atlong")

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

# Four predictors of the same target, all estimated on data through the origin:
#   mean  the recursive mean of past realisations of this target
#   rw    the last observed level of the variable (a random walk)
#   ar    the own-lag regression used so far
#   fci   that regression with the index added
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
    d <- as.data.frame(na.omit(stats::ts.intersect(
      lead = lead, own = seen, fci = path)))
    if (nrow(d) < 12) return(NULL)
    now <- data.frame(own = last_of(seen), fci = last_of(path))
    p_ <- function(f) as.numeric(predict(lm(f, data = d), newdata = now))
    data.frame(origin = set$origins[k], actual = realised,
               mean = mean(d$lead),
               rw = last_of(seen),
               ar = p_(lead ~ own),
               fci = p_(lead ~ own + fci)) }))
}

# Clark-West of `fci` against a named benchmark. The larger model nests all
# three: the mean sets the own-lag and index coefficients to zero, the random
# walk fixes the own-lag coefficient at one and the intercept at zero.
cw_vs <- function(f, bench, h) {
  eb <- f$actual - f[[bench]]; ef <- f$actual - f$fci
  adj <- eb^2 - (ef^2 - (f[[bench]] - f$fci)^2)
  mean(adj) / sqrt(hac(adj, max(h - 1, 0)) / length(adj))
}
summarise <- function(target, h, set, cumulative = TRUE) {
  f <- forecasts(target, h, set, cumulative)
  r <- function(col) sqrt(mean((f$actual - f[[col]])^2))
  b <- c(mean = r("mean"), rw = r("rw"), ar = r("ar"))
  best <- names(b)[which.min(b)]
  c(n = nrow(f), sd = sd(f$actual), mean = b[["mean"]], rw = b[["rw"]],
    ar = b[["ar"]], fci = r("fci"),
    best_rmse = min(b), ratio_best = r("fci") / min(b),
    ratio_ar = r("fci") / b[["ar"]],
    cw_best = cw_vs(f, best, h), cw_ar = cw_vs(f, "ar", h),
    best_is = match(best, c("mean", "rw", "ar")))
}

cat("building vintages ...\n"); t0 <- Sys.time()
V1 <- lapply(eamd, vintages_of, from = c(2006, 1), n = 1)
VL <- vintages_of(atlong, from = c(1999, 1), n = 1)
cat("  ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

hs <- c(1, 2, 4, 8)
grid <- expand.grid(region = names(eamd), h = hs, cumulative = c(TRUE, FALSE),
                    stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  s <- summarise("u", g$h, V1[[g$region]], g$cumulative)
  data.frame(g, t(s)) }))
names(res) <- sub("^mean$", "rmse_mean", names(res))
saveRDS(res, "bench3.rds")

show <- function(cum) {
  z <- res[res$cumulative == cum, ]
  cat(sprintf("\n--- %s target ---\n", if (cum) "PATH" else "ENDPOINT"))
  cat("region  h   sd   mean    rw     ar    fci | ratio/ar ratio/best  cw/ar cw/best  best\n")
  for (r in names(eamd)) for (h in hs) {
    x <- z[z$region == r & z$h == h, ]
    cat(sprintf("%-4s %3d %5.2f %6.3f %6.3f %6.3f %6.3f | %7.3f %9.3f %6.2f %7.2f  %s\n",
        toupper(r), h, x$sd, x$rmse_mean, x$rw, x$ar, x$fci,
        x$ratio_ar, x$ratio_best, x$cw_ar, x$cw_best,
        c("mean","rw","ar")[x$best_is]))
  }
}
show(TRUE); show(FALSE)

cat("\n=== HEADLINE: index vs the BEST of three benchmarks ===\n")
cat("           ---------- endpoint ----------   ------------ path ------------\n")
cat("  h   wins  meanratio  CW>1.645   wins  meanratio  CW>1.645\n")
for (h in hs) {
  e <- res[!res$cumulative & res$h == h, ]; p <- res[res$cumulative & res$h == h, ]
  cat(sprintf("%3d %5d %10.3f %9d %6d %10.3f %9d\n", h,
      sum(e$ratio_best < 1), mean(e$ratio_best), sum(e$cw_best > 1.645),
      sum(p$ratio_best < 1), mean(p$ratio_best), sum(p$cw_best > 1.645)))
}
cat("\nwhich benchmark wins, by horizon and target (mean/rw/ar counts):\n")
for (cum in c(FALSE, TRUE)) for (h in hs) {
  z <- res[res$cumulative == cum & res$h == h, ]
  cat(sprintf("%-8s h=%d : %s\n", if (cum) "path" else "endpoint", h,
      paste(table(factor(c("mean","rw","ar")[z$best_is],
                         levels = c("mean","rw","ar"))), collapse = " / ")))
}
