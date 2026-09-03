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
