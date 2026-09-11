# Builds `eamd`: the euro area and its ten largest member countries, as the
# European counterpart of the Koop and Korobilis (2014) dataset shipped as
# `koop`.
#
# Three sources are combined.
#
#  1. EA-MD-QD (Barigozzi and Lissona), release 2026-08. One workbook per
#     region -- EAdata.xlsx, ATdata.xlsx, ... -- each with a "data" and an
#     "info" sheet. It supplies the macroeconomic block for all eleven regions
#     and most of the financial block.
#
#  2. AustriaMacroData (https://github.com/franzmohr/AustriaMacroData),
#     output/aut_panel.csv and output/deu_panel.csv, copied to
#     data-raw/austria-macro-data/. It supplies the geopolitical risk index
#     that every region's block uses, and four extra columns for the two
#     countries it publishes a panel for. See "Country extras" below.
#
#  3. FRED, series VIXCLS, for the CBOE volatility index -- koop's own `vix`,
#     and the one series of the paper's eighteen that has a direct counterpart
#     here rather than a substitute. Fetched from the public fredgraph.csv
#     export, which needs no API key, and cached under data-raw/fred/ so the
#     build is reproducible without the network. Delete the cached file to
#     refresh it.
#
#### Transformations ####
#
# EA-MD-QD's "info" sheet carries a transformation code per series (column
# TR1), numbered
#     1 = 100 * log(x), 2 = 100 * Dlog(x), 3 = 100 * D2log(x),
#     4 = x (none),     5 = Dx
# and the block is transformed rather than shipped in levels: a factor
# extracted from trending levels is a trend, so the block has to be made
# stationary before it reaches create_fcimodel().
#
# The code applied to a concept is however the same in every region, which the
# TR1 column is not -- it is assigned per country, so REER42 is code 5
# everywhere but Spain, HHLB.SLN is code 2 everywhere but Greece and Spain, and
# so on for NFCLB.LLN, TASS.LLN and KCONFIX. Honouring each country's own code
# would make the same column mean a growth rate in one region and an
# acceleration in another, and the point of this dataset is that the eleven
# indices can be read against each other. Uniformity wins, and where the two
# rules disagree the series is differenced once rather than twice. HHLB and
# HHLB.LLN are code 3 everywhere and are also differenced once here, so that
# every stock series in the block carries one concept -- a growth rate, as
# koop's cmdebt does.
#
#### Composition ####
#
# The financial block is curated rather than exhaustive. EA-MD-QD's financial
# class is mostly sectoral balance sheets -- 26 of them for the euro area -- and
# a factor of all of it is a factor of sectoral balance sheets. What is selected
# below mirrors the composition of Koop and Korobilis's own 18 series: asset
# prices, spreads, credit growth, surveys, an uncertainty measure.

rm(list = ls())

library(dplyr)
library(readxl)

md_qd <- "data-raw/ea-md-qd/extracted/EA-MD-QD-2026-08"
amd_dir <- "data-raw/austria-macro-data"
fred_dir <- "data-raw/fred"

# The euro area aggregate first, then its ten member countries alphabetically.
regions <- c(ea = "EA", at = "AT", be = "BE", de = "DE", el = "EL", es = "ES",
             fr = "FR", ie = "IE", it = "IT", nl = "NL", pt = "PT")

#### Helpers ####

# EA-MD-QD stores quarterly series in the last month of the quarter and monthly
# series in every month. Both have to come out as one quarterly `ts`.
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

# A quarterly-native column: drop the empty months, keep what is left.
as_quarterly_native <- function(raw, id) {
  keep <- !is.na(raw[[id]])
  if (!any(keep)) stop("no observations for ", id)
  ts(raw[[id]][keep], start = c(raw$.yr[keep][1], raw$.qtr[keep][1]),
     frequency = 4)
}

# A monthly column: quarterly average, and only over complete quarters, so a
# part-published quarter is left out rather than averaged over one month.
as_quarterly_mean <- function(raw, id) {
  d <- data.frame(yr = raw$.yr, qtr = raw$.qtr, value = raw[[id]]) %>%
    filter(!is.na(value)) %>%
    group_by(yr, qtr) %>%
    filter(n() == 3) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(yr, qtr)
  if (nrow(d) == 0) stop("no complete quarter for ", id)
  ts(d$value, start = c(d$yr[1], d$qtr[1]), frequency = 4)
}

# Either of the two, chosen by the frequency the info sheet reports. Returns
# NULL when the region does not publish the series, which is how the per-region
# gaps below are handled.
series <- function(raw, info, id) {
  freq <- info$Frequency[match(id, info$Name)]
  if (is.na(freq)) return(NULL)
  if (freq == "Q") as_quarterly_native(raw, id) else as_quarterly_mean(raw, id)
}

# AustriaMacroData panels are already quarterly, one row per quarter start.
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

# One FRED series, from the public fredgraph.csv export -- no API key -- cached
# on first use. The daily observations are averaged to quarters, over quarters
# that are complete: a full quarter of trading days is 59 to 66 here (2001 Q3
# is the low one, since the market was shut for a week after 11 September),
# while a quarter still in progress has fewer, so 55 separates the two.
fred_series <- function(fred_id, min_obs = 55) {
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
  d$value <- suppressWarnings(as.numeric(d$value))  # holidays come back blank

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

# 100 * Dlog(x), the EA-MD-QD code-2 transformation, applied to every stock and
# price series in every region. Both this and the plain difference pass NULL
# through, so that a series a region does not publish stays NULL all the way to
# where the block is assembled and can be dropped there.
dlog <- function(x) if (is.null(x)) NULL else 100 * diff(log(x))
ddiff <- function(x) if (is.null(x)) NULL else diff(x)

#### Series shared by every region ####

info_ea <- read_info("EA")
raw_ea <- read_md_qd("EA")

# The money market is a single euro-area market, so the same 3- and 6-month
# rates apply to every member country; EA-MD-QD publishes them once, suffixed
# _EACC. So is the euro's dollar rate. Only the long rate is country-specific.
irt3m <- series(raw_ea, info_ea, "IRT3M_EACC")
irt6m <- series(raw_ea, info_ea, "IRT6M_EACC")
erus <- series(raw_ea, info_ea, "ERUS_EA")

# The Caldara and Iacoviello (2022) geopolitical risk index, taken from the
# Austrian panel, where it is the authors' *global* index because they publish
# no Austrian one. Used for every region here rather than the country-specific
# index where one exists, so that the eleven blocks share one uncertainty
# variable and their indices stay comparable. It is what fci() is signed on,
# and it stands in for koop's vix and move, which have no European counterpart
# in either source.
aut <- amd_panel("aut_panel.csv")
deu <- amd_panel("deu_panel.csv")
gpr <- amd_series(aut, "geopolitical_risk")

# The CBOE volatility index, in levels, as koop's `vix` is. It is US equity
# volatility, so for these eleven regions it is a global risk-aversion proxy
# rather than a domestic one -- the same standing as gpr. The euro area's own
# counterpart would be the VSTOXX, which FRED does not carry.
vix <- fred_series("VIXCLS")

#### Country extras ####
# AustriaMacroData publishes a panel for Austria and Germany but not for the
# other nine regions, and EA-MD-QD has no mortgage rate, house price index or
# credit-to-GDP ratio for any of them. So these four columns -- the closest
# counterparts either source has to koop's loanhpi and spread_30_mort -- exist
# for AT and DE only, and are appended after the common columns rather than
# mixed in among them. Subset a block to the common columns when comparing
# regions; see ?eamd.
extras_for <- list(at = aut, de = deu)

#### Build ####

build_macro <- function(raw, suffix) {
  y <- dlog(as_quarterly_native(raw, paste0("GDP", suffix)))
  p <- dlog(as_quarterly_native(raw, paste0("DFGDP", suffix)))
  u <- as_quarterly_mean(raw, paste0("UNETOT", suffix))
  out <- na.omit(ts.intersect(y, p, u))
  dimnames(out)[[2]] <- c("y", "p", "u")
  out
}

build_financial <- function(raw, info, suffix, code) {
  id <- function(x) paste0(x, suffix)
  get <- function(x) series(raw, info, id(x))

  # The euro area's own interest rates are suffixed _EACC, not _EA, so the
  # aggregate's long rate has to be looked up under the money-market suffix.
  ltirt <- get("LTIRT")
  if (is.null(ltirt)) ltirt <- series(raw, info, "LTIRT_EACC")
  if (is.null(ltirt)) stop("no long rate for ", code)

  # The Koop and Korobilis analogue of each is given in the comment. A NULL
  # element is a series the region does not publish and is dropped below.
  lst <- list(
    equity        = dlog(get("SHIX")),                       # sp500
    reer          = ddiff(get("REER42")),                    # twexmmth
    usd           = dlog(erus),                              # twexmmth
    irate_short   = diff(irt3m),
    irate_long    = diff(ltirt),
    spread_term   = ltirt - irt3m,                           # spread_10_2
    spread_money  = irt6m - irt3m,                           # spread_2_3m
    gpr           = gpr,                                     # move
    vix           = vix,                                     # vix
    hh_debt       = dlog(get("HHLB")),                       # cmdebt
    hh_loans_lt   = dlog(get("HHLB.LLN")),                   # loanhpi
    hh_loans_st   = dlog(get("HHLB.SLN")),                   # totalsl
    hh_assets     = dlog(get("HHASS")),
    nfc_debt      = dlog(get("NFCLB")),
    nfc_loans_lt  = dlog(get("NFCLB.LLN")),                  # spread_comm
    nfc_assets    = dlog(get("NFCASS")),
    bank_deposits = dlog(get("TLB.SDB")),
    bank_loans_lt = dlog(get("TASS.LLN")),
    conf_cons     = get("CCONFIX"),                          # cci
    conf_ind      = get("ICONFIX"),                          # mich1
    conf_constr   = get("KCONFIX"),                          # mich2
    esi           = get("ESENTIX")                           # mich3
  )

  panel <- extras_for[[tolower(code)]]
  if (!is.null(panel)) {
    # EA-MD-QD's ERUS is USD per EUR; AustriaMacroData's fx_rate_to_usd is
    # national currency per USD, so the mortgage rate is the only one of these
    # that needs no reorientation.
    lst$house_prices    <- dlog(amd_series(panel, "house_price_real"))
    lst$spread_mortgage <- amd_series(panel, "mortgage_rate") - ltirt
    lst$credit_gdp_hh   <- diff(amd_series(panel, "household_credit_to_gdp"))
    lst$credit_gdp_nfc  <- diff(amd_series(panel, "corporate_credit_to_gdp"))
  }

  lst <- lst[!vapply(lst, is.null, logical(1))]
  out <- do.call(ts.union, lst)
  dimnames(out)[[2]] <- names(lst)
  out
}

# Every block is cut to the span EA-MD-QD covers, 2000 Q1 onwards. The
# AustriaMacroData columns reach back to the 1950s and would otherwise pad the
# financial block with decades the macro block cannot match.
trim <- function(x, from = c(2000, 1)) {
  x <- window(x, start = from)
  x[, colSums(!is.na(x)) > 0, drop = FALSE]
}

#### Sector-account breaks ####
# National financial accounts are occasionally rebased or reclassified, and in
# a growth rate that shows up as one quarter in which a whole sector account
# jumps at once. The Dutch accounts do it in 2009 Q4 -- non-financial corporate
# assets grow 19% in a quarter against 2-4% either side -- and left alone that
# single quarter dominates the Dutch factor, which then correlates 0.10 with
# the euro area index instead of 0.87.
#
# What separates a reclassification from a genuine shock is that it moves the
# stock series together and nothing else: the 2008 Q4 crash reaches at most 3
# standard deviations in any region's balance-sheet columns, while these breaks
# reach 6 to 9 in several columns of the same quarter. So the rule is a
# simultaneity rule, applied uniformly rather than country by country -- blank
# a balance-sheet observation when it exceeds `z_min` standard deviations and
# at least `n_min` balance-sheet series in that region do so in the same
# quarter. It flags two quarters in the whole eleven-region panel and no part
# of 2008-09; the count is stable for z_min anywhere in 4 to 5.
#
# The observations are set to NA rather than interpolated, because that is what
# create_fcimodel() is built to carry: the filter drops a series over the
# periods it is missing instead of imputing anything.
stock_cols <- c("hh_debt", "hh_loans_lt", "hh_loans_st", "hh_assets",
                "nfc_debt", "nfc_loans_lt", "nfc_assets", "bank_deposits",
                "bank_loans_lt")

blank_breaks <- function(x, region, z_min = 5, n_min = 3) {
  cols <- intersect(stock_cols, colnames(x))
  z <- scale(x[, cols, drop = FALSE])
  hit <- abs(z) > z_min
  hit[is.na(hit)] <- FALSE
  for (i in which(rowSums(hit) >= n_min)) {
    when <- as.numeric(stats::time(x))[i]
    yr <- floor(when + 1e-8)
    cat(sprintf("  break blanked: %-3s %d Q%d -- %s\n", toupper(region), yr,
                round((when - yr) * 4) + 1, paste(cols[hit[i, ]], collapse = " ")))
    x[i, cols[hit[i, ]]] <- NA
  }
  x
}

eamd <- list()
cat("\n")
for (r in names(regions)) {
  code <- regions[[r]]
  raw <- if (code == "EA") raw_ea else read_md_qd(code)
  info <- if (code == "EA") info_ea else read_info(code)
  # The euro area's own series are suffixed _EA; a country's, _AT, _BE, ...
  fin <- trim(build_financial(raw, info, paste0("_", code), code))
  eamd[[r]] <- list(
    financial = blank_breaks(fin, r),
    macro = build_macro(raw, paste0("_", code))
  )
}

#### Report ####

# The 21 columns the block is defined to have, whether or not a given region
# publishes all of them.
core <- c("equity", "reer", "usd", "irate_short", "irate_long", "spread_term",
          "spread_money", "gpr", "vix", "hh_debt", "hh_loans_lt",
          "hh_loans_st", "hh_assets", "nfc_debt", "nfc_loans_lt",
          "nfc_assets", "bank_deposits", "bank_loans_lt", "conf_cons",
          "conf_ind", "conf_constr", "esi")

shared <- Reduce(intersect, lapply(eamd, function(z) colnames(z$financial)))
cat("\ncore columns:", length(core),
    "| present in every region:", length(shared), "\n\n")

summary_tbl <- do.call(rbind, lapply(names(eamd), function(r) {
  f <- eamd[[r]]$financial
  m <- eamd[[r]]$macro
  data.frame(
    region = r,
    n_fin = ncol(f),
    missing = paste(setdiff(core, colnames(f)), collapse = " "),
    extra = paste(setdiff(colnames(f), core), collapse = " "),
    macro = paste0(paste(start(m), collapse = "Q"), "-",
                   paste(end(m), collapse = "Q")),
    n_obs = nrow(m)
  )
}))
print(summary_tbl, row.names = FALSE)

stopifnot(identical(sort(shared), sort(intersect(core, shared))))

usethis::use_data(eamd, overwrite = TRUE)
