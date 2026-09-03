#' Financial Conditions Index (Koop and Korobilis, 2014)
#'
#' @description
#' Deprecated. Superseded by the model object workflow. \code{kk_algorithm} was never
#' finished and returned the principal components of \code{x} rather than an
#' estimated index; it is kept only so that old scripts fail with a message that
#' says what to call instead.
#'
#' @param y a time-series object of macroeconomic variables.
#' @param x a time-series object of financial variables.
#' @param ... further arguments, ignored.
#'
#' @details Use the three step workflow instead:
#'
#' \preformatted{
#' model <- create_fcimodel(y = y, x = x, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#' index <- fci(model)
#' }
#'
#' @return Nothing. The function always throws.
#'
#' @seealso \code{\link{create_fcimodel}}, \code{\link{fci}}
#'
#' @keywords internal
#' @export
kk_algorithm <- function(y, x, ...) {
  stop("'kk_algorithm' was removed in fincond 0.1.0 and never produced an ",
       "index. Use create_fcimodel(y = y, x = x, p = 4, n = 1), then ",
       "add_priors(), then add_posterior_coefficients(), then fci().")
}
