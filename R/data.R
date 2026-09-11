#' Data from Koop and Korobilis (2014)
#'
#' The macroeconomic and financial variables behind the financial conditions
#' index of Koop and Korobilis (2014), quarterly, from 1970 Q1.
#'
#' @usage data("koop")
#'
#' @format A named list with two elements:
#' \describe{
#'   \item{financial}{A time-series object of 18 financial variables, 1970 Q1 to
#'   2013 Q3. Several of them start later than the sample does and are \code{NA}
#'   until they do -- the VIX from 1990, for instance -- which is what
#'   \code{\link{create_fcimodel}} carries through to the filter rather than
#'   imputing.}
#'   \item{macro}{A list of 96 real-time vintages of three macroeconomic
#'   variables: inflation (\code{p}), the unemployment rate (\code{u}) and
#'   output (\code{y}). Each vintage begins in 1970 Q1 and ends one quarter later
#'   than the one before it, so \code{macro[[1]]} ends in 1989 Q4 and
#'   \code{macro[[96]]} in 2013 Q3.}
#' }
#'
#' @details Use \code{macro[[96]]} for an index over the whole sample. The
#' earlier vintages are what a recursive, real-time exercise runs over, and note
#' that the short ones cannot be combined with the full set of financial
#' variables: over the periods \code{macro[[1]]} covers, the series that start in
#' the 1990s have no observations at all.
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
"koop"


#' Euro area and member country data
#'
#' Macroeconomic and financial variables for the euro area and its ten largest
#' member countries, built to the design of \code{\link{koop}} so that the index
#' of Koop and Korobilis (2014) can be estimated for each of them on comparable
#' data, quarterly, from 2000 Q1.
#'
#' @usage data("eamd")
#'
#' @format A named list of eleven regions -- \code{ea} (the euro area
#' aggregate), then \code{at}, \code{be}, \code{de}, \code{el}, \code{es},
#' \code{fr}, \code{ie}, \code{it}, \code{nl} and \code{pt}. Each is a list of
#' two time-series objects:
#' \describe{
#'   \item{financial}{The block the index is extracted from, 2000 Q1 to 2026 Q3.
#'   The first 22 columns are one common definition applied to every region and
#'   carry the same name throughout, so the eleven indices can be read against
#'   each other; see 'Details' for the four regions that do not have all of
#'   them, and for the two that have more. Missing values are carried through to
#'   the filter by \code{\link{create_fcimodel}} rather than imputed.}
#'   \item{macro}{Three macroeconomic variables, 2000 Q2 to 2025 Q4: output
#'   growth (\code{y}), inflation (\code{p}) and the unemployment rate
#'   (\code{u}), in the units \code{koop$macro} uses. Unlike \code{koop} there
#'   are no real-time vintages -- neither source publishes them -- so this is a
#'   single, fully revised series rather than a list.}
#' }
#'
#' @details The financial block is stationary as shipped. EA-MD-QD publishes a
#' transformation code per series and \code{data-raw/eamd.R} applies one:
#' growth rates for stocks and prices, first differences for interest rates and
#' the effective exchange rate, levels for survey balances, the spreads and the
#' uncertainty index. This matters, because a factor extracted from trending
#' levels is a trend.
#'
#' The code applied to a concept is the same in every region, which EA-MD-QD's
#' own TR1 column is not -- it is assigned per country, so the same series is
#' differenced once for one member and twice for another. Honouring that would
#' make a column mean a growth rate in one region and an acceleration in
#' another, which would defeat the purpose of a comparable panel, so the block
#' is differenced once throughout.
#'
#' It is also curated rather than exhaustive: EA-MD-QD's financial class is
#' mostly sectoral balance sheets, and a factor of all of it is a factor of
#' sectoral balance sheets. What is kept mirrors the composition of the 18
#' series of Koop and Korobilis -- asset prices, spreads, credit growth, surveys
#' and one uncertainty measure.
#'
#' Four regions depart from the common 21 columns. Ireland has no split of
#' household or corporate loans by maturity, so it lacks \code{hh_loans_lt},
#' \code{hh_loans_st} and \code{nfc_loans_lt}; Greece lacks
#' \code{bank_deposits}. Austria and Germany carry four columns the others do
#' not -- \code{house_prices}, \code{spread_mortgage}, \code{credit_gdp_hh} and
#' \code{credit_gdp_nfc} -- because AustriaMacroData publishes a national panel
#' for those two and neither source has the equivalents elsewhere. They are the
#' closest counterparts in this data to \code{koop}'s \code{loanhpi} and
#' \code{spread_30_mort}, and they are appended after the common columns, so
#' subsetting a block to the first 22 names recovers the strictly comparable
#' set.
#'
#' The money market is a single euro-area market, so \code{irate_short} and
#' \code{spread_money} are the same series in all eleven blocks, as is
#' \code{usd}; only \code{irate_long}, and through it \code{spread_term}, is
#' country-specific.
#'
#' Two uncertainty measures are common to every block. \code{vix} is
#' \code{koop}'s own series, the CBOE volatility index, taken from FRED
#' (\code{VIXCLS}) and averaged from daily closes over complete quarters;
#' \code{gpr} is the Caldara and Iacoviello (2022) geopolitical risk index, the
#' authors' \emph{global} one rather than the country-specific version where
#' they publish it. Both are global rather than domestic measures here, so
#' either can move any of these indices without anything having happened in that
#' economy.
#'
#' Neither is a good anchor for \code{\link{fci}}'s \code{sign_on}, which is
#' worth knowing before reaching for \code{sign_on = "vix"} by analogy with the
#' paper. The VIX loads 0.06 in absolute value in Italy, 0.21 in France and 0.38
#' in Germany, and it loads with the \emph{wrong} sign in Greece, Italy and
#' Portugal, whose factor is dominated by a sovereign crisis that US equity
#' volatility did not track; signing on it there turns those three indices
#' upside down. \code{gpr} is weak too, at 0.34 in Spain. The survey balances
#' are what load strongly everywhere -- \code{esi} at no less than 0.82 -- but
#' they rise with \emph{loosening}, so \code{-fci(model, sign_on = "esi")} is
#' the robust route to a stress-oriented index, and is what
#' \code{vignette("euro-area-and-member-countries")} uses.
#'
#' National financial accounts are occasionally reclassified, and in a growth
#' rate that shows up as one quarter in which a whole sector account jumps at
#' once. Two such quarters are blanked -- Austria in 2005 Q4 and the
#' Netherlands in 2009 Q4 -- by a uniform rule: a balance-sheet observation
#' beyond five standard deviations is set to \code{NA} when at least three
#' balance-sheet series in that region exceed it in the same quarter. No part of
#' 2008--09 is caught, since the crash reaches at most three standard deviations
#' in these columns. The Dutch break matters: left in, that single quarter
#' dominates the Dutch factor, whose correlation with the euro area index is
#' then 0.10 rather than 0.85.
#'
#' @source EA-MD-QD, release 2026-04, of Barigozzi and Lissona; the Austrian and
#' German panels of \url{https://github.com/franzmohr/AustriaMacroData}; and
#' series \code{VIXCLS} from FRED,
#' \url{https://fred.stlouisfed.org/series/VIXCLS}. See \code{data-raw/eamd.R}
#' for how each column is built.
#'
#' @references
#'
#' Barigozzi, M., & Lissona, C. (2024). EA-MD-QD: Large Euro Area and Euro
#' Member Countries Datasets for Macroeconomic Research. Zenodo.
#' \doi{10.5281/zenodo.10514667}
#'
#' Caldara, D., & Iacoviello, M. (2022). Measuring geopolitical risk.
#' \emph{American Economic Review, 112}(4), 1194--1225.
#' \doi{10.1257/aer.20191823}
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
"eamd"


#' Austria over a long sample
#'
#' Macroeconomic and financial variables for Austria from 1993 Q1, built to the
#' design of \code{\link{koop}} and \code{\link{eamd}} but trading the
#' cross-country comparability of the latter for thirty more quarters of one
#' country, quarterly.
#'
#' @usage data("atlong")
#'
#' @format A list of two time-series objects:
#' \describe{
#'   \item{financial}{27 financial variables, 1993 Q1 to 2026 Q3. The block is
#'   ragged by construction: eleven series are present in the first quarter and
#'   the rest begin later, the sectoral balance sheets last, in 2000. Missing
#'   values are carried through to the filter by \code{\link{create_fcimodel}}
#'   rather than imputed, which is the same shape \code{koop} has -- its VIX
#'   starts in 1990 against a sample from 1970.}
#'   \item{macro}{Output growth (\code{y}), inflation (\code{p}) and the
#'   unemployment rate (\code{u}), 1993 Q1 to 2025 Q4, in the units
#'   \code{koop$macro} uses. Complete by construction: the macro block is the
#'   state itself and \code{create_fcimodel} does not accept gaps in it.}
#' }
#'
#' @details \code{eamd} begins in 2000 Q1 because EA-MD-QD does, for all eleven
#' of its regions, and it cannot be extended for one member without giving up
#' what makes it comparable. This dataset gives that up deliberately. Austria is
#' one of the two regions where it is possible, because AustriaMacroData
#' publishes a national panel reaching back to the 1950s; the other nine regions
#' have nothing before 2000 in either source. Use \code{eamd} for anything
#' cross-country.
#'
#' What sets the start is the unemployment rate, which neither source carries
#' before 1993 Q1 in a form that also runs to the present. Real GDP reaches back
#' to 1960 and would allow considerably more.
#'
#' Two things differ from \code{eamd} beyond the sample. Inflation is the
#' consumer price index rather than the GDP deflator, because the national panel
#' publishes no deflator; its own CPI starts in 1996, so the twelve quarters
#' before that come from FRED's OECD series \code{AUTCPIALLQINMEI}. Growth rates
#' are spliced rather than levels, so no chain-linking is involved and no level
#' break can be introduced, and the two sources correlate 0.92 in quarterly
#' inflation over their 116-quarter overlap. And where a concept exists in both
#' sources the national one is used throughout rather than switching in 2000, so
#' that a column means one thing over the whole sample -- which makes the short
#' rate Austria's own money-market rate before 1999 and the euro rate after.
#'
#' The two are close where they overlap despite those differences: the index of
#' Koop and Korobilis (2014) estimated on this block correlates 0.98 with the
#' one estimated on \code{eamd$at} over the 103 quarters they share, and both
#' peak in 2009 Q1. The extra years read as the ERM crisis and the Austrian
#' recession of 1993, which is the tightest reading in the sample outside 2008
#' and 2001.
#'
#' One quarter is blanked, Austria 2005 Q4, by the same simultaneity rule
#' \code{data-raw/eamd.R} applies and for the same reason; see \code{?eamd}.
#'
#' @source AustriaMacroData, \url{https://github.com/franzmohr/AustriaMacroData};
#' EA-MD-QD, release 2026-04, of Barigozzi and Lissona, for the sectoral balance
#' sheets and the money-market spread from 2000; and series \code{VIXCLS} and
#' \code{AUTCPIALLQINMEI} from FRED. See \code{data-raw/atlong.R} for how each
#' column is built.
#'
#' @references
#'
#' Barigozzi, M., & Lissona, C. (2024). EA-MD-QD: Large Euro Area and Euro
#' Member Countries Datasets for Macroeconomic Research. Zenodo.
#' \doi{10.5281/zenodo.10514667}
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
"atlong"
