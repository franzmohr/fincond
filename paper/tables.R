# Turns the output of paper/results.R into the four tables of the manuscript.
# Run from the directory holding bench3.rds.

d <- readRDS("bench3.rds")
nm <- c(ea = "Euro area", at = "Austria", be = "Belgium", de = "Germany",
        el = "Greece", es = "Spain", fr = "France", ie = "Ireland",
        it = "Italy", nl = "Netherlands", pt = "Portugal")
p <- d[d$cumulative, ]
e <- d[!d$cumulative, ]

cat("=== T1: how bad is the AR benchmark? (path target) ===\n")
cat("region   h=1              h=2              h=4              h=8   (ar/rw ratio)\n")
for (r in names(nm)) {
  v <- vapply(c(1,2,4,8), function(h) {
    x <- p[p$region == r & p$h == h, ]; x$ar / x$rw }, numeric(1))
  cat(sprintf("%-12s %6.2f %16.2f %16.2f %16.2f\n", nm[r], v[1], v[2], v[3], v[4]))
}
cat("\nbenchmark chosen, counts of (mean/rw/ar) by horizon, path target:\n")
for (h in c(1,2,4,8)) {
  z <- p[p$h == h, ]
  cat(sprintf("h=%d : %s\n", h, paste(table(factor(c("mean","rw","ar")[z$best_is],
      levels = c("mean","rw","ar"))), collapse = " / ")))
}
cat("\nmedian and max ar/rw ratio by horizon:\n")
for (h in c(1,2,4,8)) {
  z <- p[p$h == h, ]
  cat(sprintf("h=%d  median %.2f  max %.2f (%s)\n", h, median(z$ar/z$rw),
      max(z$ar/z$rw), nm[z$region[which.max(z$ar/z$rw)]]))
}

cat("\n\n=== T2: index vs best benchmark, path target (ratio | CW) ===\n")
for (r in names(nm)) {
  rr <- vapply(c(1,2,4,8), function(h) p$ratio_best[p$region==r & p$h==h], numeric(1))
  cc <- vapply(c(1,2,4,8), function(h) p$cw_best[p$region==r & p$h==h], numeric(1))
  bb <- vapply(c(1,2,4,8), function(h)
    c("mean","rw","ar")[p$best_is[p$region==r & p$h==h]], character(1))
  cat(sprintf("%-12s %5.3f %5.3f %5.3f %5.3f  & %5.2f %5.2f %5.2f %5.2f   %s\n",
      nm[r], rr[1],rr[2],rr[3],rr[4], cc[1],cc[2],cc[3],cc[4], paste(bb, collapse=" ")))
}
cat("\nsummary, path target:\n")
for (h in c(1,2,4,8)) {
  z <- p[p$h == h, ]
  cat(sprintf("h=%d  wins %2d/11  mean %.3f  median %.3f  CW>1.645 %2d/11\n",
      h, sum(z$ratio_best < 1), mean(z$ratio_best), median(z$ratio_best),
      sum(z$cw_best > 1.645)))
}

cat("\n\n=== T3: what the own-lag benchmark alone would have shown ===\n")
for (h in c(1,2,4,8)) {
  z <- p[p$h == h, ]
  cat(sprintf("h=%d  vs AR: wins %2d/11 mean %.3f CW>1.645 %2d | vs best: wins %2d/11 mean %.3f CW>1.645 %2d\n",
      h, sum(z$ratio_ar < 1), mean(z$ratio_ar), sum(z$cw_ar > 1.645),
      sum(z$ratio_best < 1), mean(z$ratio_best), sum(z$cw_best > 1.645)))
}

cat("\n\n=== T4: target definition, against the best benchmark ===\n")
for (h in c(1,2,4,8)) {
  ze <- e[e$h == h, ]; zp <- p[p$h == h, ]
  cat(sprintf("h=%d  endpoint: wins %2d mean %.3f CW %2d | path: wins %2d mean %.3f CW %2d\n",
      h, sum(ze$ratio_best < 1), mean(ze$ratio_best), sum(ze$cw_best > 1.645),
      sum(zp$ratio_best < 1), mean(zp$ratio_best), sum(zp$cw_best > 1.645)))
}
cat("\nn per region-horizon:", paste(sort(unique(p$n)), collapse = " "), "\n")
