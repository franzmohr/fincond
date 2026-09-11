# Builds `atlong`: Austria alone, over as long a sample as the data allow.
#
# `eamd` stops at 2000 Q1 because EA-MD-QD does: every one of its eleven
# workbooks begins in 2000-01, and the point of that dataset is that the eleven
# regions share one span and one definition, so it cannot be extended for one
# member without giving up what makes it comparable. This dataset gives that
# comparability up deliberately, in exchange for thirty more quarters of one
# country. `eamd` is untouched; use it for anything cross-country.
#
# Austria is the country this is possible for because AustriaMacroData
# publishes a national panel reaching back to the 1950s. Germany has one too,
# so the same script would work for it; the other nine regions have nothing
# before 2000 in either source.
#
#### What sets the start ####
#
# create_fcimodel() tolerates missing values in the financial block -- it drops
# a series over the periods it lacks rather than imputing -- but not in the
# macro block, which is the state itself. So the sample begins where the last
# of the three macro variables does, and that is the unemployment rate, which
# neither AustriaMacroData nor FRED carries before 1993 Q1 in a form that also
# runs to the present. Real GDP reaches back to 1960 and would allow far more.
#
# The financial block is therefore ragged by construction, which is the same
# shape koop has -- its VIX starts in 1990 against a 1970 sample -- and is what
# the filter is built to carry. Eight series are present in 1993 and the rest
# arrive as they begin, the sectoral balance sheets last, in 2000.
#
#### Two departures from eamd ####
#
#  1. Inflation is the consumer price index, where eamd uses the GDP deflator.
#     The national panel publishes no deflator and its CPI starts only in 1996,
#     so the twelve quarters before that are taken from FRED's OECD series
#     AUTCPIALLQINMEI and the national index is used from 1996 Q2 onward.
#     Growth rates are spliced rather than levels, so no chain-linking is
#     involved and no level break can be introduced; the two sources correlate
#     0.92 in quarterly inflation over their 116-quarter overlap.
#
#  2. Where a concept exists in both sources, the national one is used over the
#     whole sample rather than switching to EA-MD-QD in 2000, so that a column
#     means one thing throughout. That is why the short rate is Austria's own
#     money-market rate before 1999 and the euro rate after: the correct series
#     for a national index, and not the same construction eamd uses.

rm(list = ls())

library(dplyr)
library(readxl)

md_qd <- "data-raw/ea-md-qd/extracted/EA-MD-QD-2026-08"
amd_dir <- "data-raw/austria-macro-data"
fred_dir <- "data-raw/fred"

start_q <- c(1993, 1)

#### Helpers ####

amd_panel <- function(file) {
  d <- read.csv(file.path(amd_dir, file), stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d
}

amd_series <- function(panel, id) {
  yr <- as.integer(format(panel$date[1], "%Y"))
  qtr <- (as.integer(format(panel$date[1], "%m")) - 1L) %/% 3L + 1L
  ts(panel[[id]], start = c(yr, qtr), frequency = 4)
}

# A FRED series that is already quarterly, unlike VIXCLS which is daily.
fred_quarterly <- function(fred_id) {
  file <- file.path(fred_dir, paste0(fred_id, ".csv"))
  if (!file.exists(file)) {
    dir.create(fred_dir, showWarnings = FALSE, recursive = TRUE)
    utils::download.file(
      paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", fred_id),
      destfile = file, quiet = TRUE)
  }
  d <- read.csv(file, stringsAsFactors = FALSE)
  names(d) <- c("date", "value")
  d$date <- as.Date(d$date)
  ts(suppressWarnings(as.numeric(d$value)),
     start = c(as.integer(format(d$date[1], "%Y")),
               (as.integer(format(d$date[1], "%m")) - 1L) %/% 3L + 1L),
     frequency = 4)
}

# VIXCLS is daily; quarters still in progress are dropped rather than averaged
# over part of themselves. Same rule as data-raw/eamd.R.
fred_daily_mean <- function(fred_id, min_obs = 55) {
  file <- file.path(fred_dir, paste0(fred_id, ".csv"))
  if (!file.exists(file)) {
    dir.create(fred_dir, showWarnings = FALSE, recursive = TRUE)
    utils::download.file(
      paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", fred_id),
      destfile = file, quiet = TRUE)
  }
  d <- read.csv(file, stringsAsFactors = FALSE)
  names(d) <- c("date", "value")
  d$date <- as.Date(d$date)
  d$value <- suppressWarnings(as.numeric(d$value))
  d <- d %>%
    mutate(yr = as.integer(format(date, "%Y")),
           qtr = (as.integer(format(date, "%m")) - 1L) %/% 3L + 1L) %>%
    filter(!is.na(value)) %>%
    group_by(yr, qtr) %>%
    filter(n() >= min_obs) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(yr, qtr)
  ts(d$value, start = c(d$yr[1], d$qtr[1]), frequency = 4)
}

read_md_qd <- function(code) {
  raw <- as.data.frame(read_excel(file.path(md_qd, paste0(code, "data.xlsx")),
                                  sheet = "data"))
  raw$Time <- as.Date(raw$Time)
  raw$.yr <- as.integer(format(raw$Time, "%Y"))
  raw$.qtr <- (as.integer(format(raw$Time, "%m")) - 1L) %/% 3L + 1L
  raw
}

read_info <- function(code) {
  as.data.frame(read_excel(file.path(md_qd, paste0(code, "data.xlsx")),
                           sheet = "info"))
}

as_quarterly_native <- function(raw, id) {
  keep <- !is.na(raw[[id]])
  ts(raw[[id]][keep], start = c(raw$.yr[keep][1], raw$.qtr[keep][1]),
     frequency = 4)
}

as_quarterly_mean <- function(raw, id) {
  d <- data.frame(yr = raw$.yr, qtr = raw$.qtr, value = raw[[id]]) %>%
    filter(!is.na(value)) %>%
    group_by(yr, qtr) %>%
    filter(n() == 3) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(yr, qtr)
  ts(d$value, start = c(d$yr[1], d$qtr[1]), frequency = 4)
}

series <- function(raw, info, id) {
  freq <- info$Frequency[match(id, info$Name)]
  if (is.na(freq)) return(NULL)
  if (freq == "Q") as_quarterly_native(raw, id) else as_quarterly_mean(raw, id)
}

dlog <- function(x) if (is.null(x)) NULL else 100 * diff(log(x))
ddiff <- function(x) if (is.null(x)) NULL else diff(x)

# Prefer `first`, fall back to `second` only where `first` has nothing. Used
# once, for the pre-1996 stretch of the CPI.
prefer <- function(first, second) {
  both <- ts.union(first, second)
  out <- both[, 1]
  gap <- is.na(out) & !is.na(both[, 2])
  out[gap] <- both[gap, 2]
  ts(out, start = start(both), frequency = 4)
}

#### Sources ####

aut <- amd_panel("aut_panel.csv")
raw_at <- read_md_qd("AT")
info_at <- read_info("AT")
raw_ea <- read_md_qd("EA")
info_ea <- read_info("EA")

#### Macro ####

# Inflation: national CPI where it exists, OECD via FRED before it. Growth
# rates are spliced, not levels; see the header.
infl <- prefer(dlog(amd_series(aut, "cpi_index")),
               dlog(fred_quarterly("AUTCPIALLQINMEI")))

macro <- na.omit(ts.intersect(
  y = dlog(amd_series(aut, "real_gdp")),
  p = infl,
  u = amd_series(aut, "unemployment_rate")
))
dimnames(macro)[[2]] <- c("y", "p", "u")
macro <- window(macro, start = start_q)

stopifnot(!anyNA(macro))

#### Financial ####

short <- amd_series(aut, "short_term_rate")
long <- amd_series(aut, "long_term_rate")

# The euro-area money market, which is a single market and so is not national.
# EA-MD-QD is the only source for it and it starts in 2000, like everything
# else there.
irt3m <- series(raw_ea, info_ea, "IRT3M_EACC")
irt6m <- series(raw_ea, info_ea, "IRT6M_EACC")

at <- function(x) series(raw_at, info_at, paste0(x, "_AT"))

lst <- list(
  # Long national series, present well before 2000.
  equity          = dlog(amd_series(aut, "share_price_index")),
  reer            = ddiff(amd_series(aut, "real_effective_exchange_rate")),
  # AustriaMacroData quotes national currency per USD, the reciprocal of the
  # USD-per-EUR rate eamd uses, so the growth rate is negated to keep the
  # column meaning what it means there.
  usd             = -dlog(amd_series(aut, "fx_rate_to_usd")),
  irate_short     = ddiff(short),
  irate_long      = ddiff(long),
  spread_term     = long - short,
  gpr             = amd_series(aut, "geopolitical_risk"),
  vix             = fred_daily_mean("VIXCLS"),
  house_prices    = dlog(amd_series(aut, "house_price_real")),
  conf_ind        = amd_series(aut, "industrial_confidence"),
  esi             = amd_series(aut, "economic_sentiment_indicator"),
  conf_cons       = amd_series(aut, "consumer_confidence"),
  conf_constr     = amd_series(aut, "construction_confidence"),
  credit_private  = dlog(amd_series(aut, "credit_to_private_nonfin_sector")),
  credit_gdp_hh   = ddiff(amd_series(aut, "household_credit_to_gdp")),
  credit_gdp_nfc  = ddiff(amd_series(aut, "corporate_credit_to_gdp")),
  spread_mortgage = amd_series(aut, "mortgage_rate") - long,

  # From EA-MD-QD, so 2000 onwards only. These are the sectoral balance sheets
  # and the money-market spread, which no national source carries.
  spread_money    = irt6m - irt3m,
  hh_debt         = dlog(at("HHLB")),
  hh_loans_lt     = dlog(at("HHLB.LLN")),
  hh_loans_st     = dlog(at("HHLB.SLN")),
  hh_assets       = dlog(at("HHASS")),
  nfc_debt        = dlog(at("NFCLB")),
  nfc_loans_lt    = dlog(at("NFCLB.LLN")),
  nfc_assets      = dlog(at("NFCASS")),
  bank_deposits   = dlog(at("TLB.SDB")),
  bank_loans_lt   = dlog(at("TASS.LLN"))
)

lst <- lst[!vapply(lst, is.null, logical(1))]
financial <- do.call(ts.union, lst)
dimnames(financial)[[2]] <- names(lst)
financial <- window(financial, start = start_q)
financial <- financial[, colSums(!is.na(financial)) > 0, drop = FALSE]

#### Sector-account breaks ####
# The same simultaneity rule data-raw/eamd.R applies, and for the same reason:
# a reclassification moves a whole sector account at once and nothing else.
# Only the EA-MD-QD stock columns can show it, so only they are tested.
stock_cols <- intersect(
  c("hh_debt", "hh_loans_lt", "hh_loans_st", "hh_assets", "nfc_debt",
    "nfc_loans_lt", "nfc_assets", "bank_deposits", "bank_loans_lt"),
  colnames(financial))

z <- scale(financial[, stock_cols, drop = FALSE])
hit <- abs(z) > 5
hit[is.na(hit)] <- FALSE
for (i in which(rowSums(hit) >= 3)) {
  when <- as.numeric(stats::time(financial))[i]
  yr <- floor(when + 1e-8)
  cat(sprintf("  break blanked: %d Q%d -- %s\n", yr, round((when - yr) * 4) + 1,
              paste(stock_cols[hit[i, ]], collapse = " ")))
  financial[i, stock_cols[hit[i, ]]] <- NA
}

#### Assemble ####

atlong <- list(financial = financial, macro = macro)

#### Report ####

cat("\nmacro: ", paste(start(macro), collapse = " Q"), " to ",
    paste(end(macro), collapse = " Q"), " (", nrow(macro), " quarters)\n",
    sep = "")
cat("financial: ", ncol(financial), " columns, ",
    paste(start(financial), collapse = " Q"), " to ",
    paste(end(financial), collapse = " Q"), "\n\n", sep = "")

firsts <- sapply(colnames(financial), function(cl) {
  i <- which(!is.na(financial[, cl]))[1]
  when <- as.numeric(stats::time(financial))[i]
  yr <- floor(when + 1e-8)
  paste0(yr, " Q", round((when - yr) * 4) + 1)
})
print(data.frame(starts = firsts))

cat("\nseries available in the first quarter:",
    sum(!is.na(financial[1, ])), "of", ncol(financial), "\n")

usethis::use_data(atlong, overwrite = TRUE)
