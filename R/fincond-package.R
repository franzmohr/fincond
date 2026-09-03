#' fincond: Financial Conditions Indices
#'
#' Estimates a financial conditions index (FCI) from the time varying parameter
#' FAVAR of Koop and Korobilis (2014).
#'
#' @section The workflow:
#'
#' The package follows the conventions of \pkg{bvartools} and uses its generics,
#' so a model is built up in the same three steps:
#'
#' \preformatted{
#' model <- create_fcimodel(y = macro, x = financial, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#' index <- fci(model)
#' }
#'
#' @section What the estimator is, and is not:
#'
#' The model is a FAVAR whose loadings and transition coefficients both drift as
#' random walks. It is estimated by two Kalman passes -- the parameters given the
#' factors, then the factors given the parameters -- and not by posterior
#' simulation.
#'
#' That makes it \emph{approximately} Bayesian, and the approximation is worth
#' being explicit about. The Kalman recursions themselves are exact Bayesian
#' updating and the initial conditions are genuine priors. What is not Bayesian
#' is everything to do with the variances: the covariance matrices of the two
#' error terms are exponentially weighted moving averages of past residuals, and
#' the variance of the drift in the two coefficient blocks is not a parameter at
#' all but whatever the forgetting factors imply. Those are plug-in point
#' estimates substituted into the filter, and their uncertainty is discarded.
#'
#' The result is one path per quantity rather than a posterior distribution over
#' paths. This package therefore returns time series, not \code{coda} chains, and
#' nothing it returns should be read as a credible interval. Removing exactly
#' those priors is what makes the estimator fast enough to run over hundreds of
#' thousands of candidate models, which is what Koop and Korobilis (2014) do and
#' what a simulation based dynamic factor model cannot do.
#'
#' For a fully Bayesian dynamic factor model, with priors on the innovation
#' variances and posterior draws to go with them, see the \pkg{dfmtools} package.
#'
#' @author Franz X. Mohr
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
#' Raftery, A. E., Kárný, M., & Ettler, P. (2010). Online prediction under model
#' uncertainty via dynamic model averaging. \emph{Technometrics, 52}(1), 52--66.
#'
#' West, M., & Harrison, J. (1997). \emph{Bayesian forecasting and dynamic
#' models} (2nd ed.). New York: Springer.
#'
#' @useDynLib fincond, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom bvartools add_priors add_posterior_coefficients
#' @keywords internal
"_PACKAGE"
