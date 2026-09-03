#' Estimate a Financial Conditions Index Model
#'
#' Runs the two-step algorithm of Koop and Korobilis (2014) and attaches the
#' estimated paths to the model object.
#'
#' @param object an object of class \code{"fcimodel"}, usually the result of a
#' call to \code{\link{add_priors.fcimodel}}.
#' @param posterior_function a function that estimates the model. Used to
#' override the built-in estimator; it is passed \code{object} and must return
#' the same shape.
#' @param ... not used.
#'
#' @details The name follows the \pkg{bvartools} generic it implements, but the
#' estimator does not simulate and what it returns is not a posterior sample. Each
#' element is one path of length \eqn{T}, the mean of the smoothed distribution
#' at each period conditional on the discounted variances -- one number per
#' period, not a chain. They are returned as time-series objects rather than as
#' \code{coda} chains for that reason, and no interval computed from them is a
#' credible interval. See \code{\link{fincond}} for what the approximation
#' involves.
#'
#' The two steps are the parameters given the factors and then the factors given
#' the parameters, each a Kalman filter followed by a backward smoothing pass.
#' For the two coefficient blocks that backward pass simplifies exactly: where
#' the predicted state covariance is the previous one divided by a forgetting
#' factor \eqn{\kappa}, the Rauch-Tung-Striebel gain is \eqn{\kappa} times the
#' identity, so smoothing a random walk under discounting is backward exponential
#' smoothing and no covariance path has to be stored.
#'
#' @return The model object, with the additional element \code{posterior}, a list
#' containing
#' \item{factors}{A time-series object of the estimated factors, \eqn{T \times n}.
#' The financial conditions index is the first column; see \code{\link{fci}}.}
#' \item{lambda}{A time-series object of the loadings, one row per period. Within
#' a row the elements run along the rows of the \eqn{n_x \times n_z} loading
#' matrix.}
#' \item{beta}{A time-series object of the transition coefficients, one row per
#' period, vectorised by column as in \pkg{bvartools}.}
#' \item{u_sigma}{A time-series object of the diagonal of the measurement error
#' covariance.}
#' \item{v_sigma}{A time-series object of the diagonal of the transition error
#' covariance.}
#' \item{loglik}{A time-series object of the one step ahead predictive log
#' density of the observed data. This is what a dynamic model averaging layer
#' weights competing specifications by.}
#'
#' @examples
#'
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
#' @export
add_posterior_coefficients.fcimodel <- function(object, posterior_function = NULL,
                                                ...) {

  class_of_object <- class(object)

  # Copy in case the estimation fails
  model <- object
  if ("posterior" %in% names(model)) {
    model[["posterior"]] <- NULL
  }

  if (is.null(posterior_function)) {
    object <- try({
      if (is.null(object[["priors"]])) {
        stop("Element 'priors' is missing. Did you call add_priors?")
      }
      algorithm <- object[["model"]][["algorithm"]]
      if (is.null(algorithm)) {
        stop("Element 'model$algorithm' is missing. Was the object produced by ",
             "create_fcimodel?")
      }
      if (algorithm != "FavarKk2014") {
        stop("Algorithm '", algorithm, "' not supported.")
      }

      object <- .kk_favar(object)

      # Every estimated quantity is a path over the sample, so every one of them
      # is given the sample's dates.
      tsp_y <- stats::tsp(model[["data"]][["y"]])
      for (i in c("factors", "lambda", "beta", "u_sigma", "v_sigma", "loglik")) {
        object[["posterior"]][[i]] <-
          stats::ts(object[["posterior"]][[i]],
                    start = tsp_y[1], frequency = tsp_y[3])
      }
      dimnames(object[["posterior"]][["factors"]])[[2]] <-
        paste0("factor", seq_len(object[["model"]][["n"]]))
      dimnames(object[["posterior"]][["u_sigma"]])[[2]] <-
        dimnames(model[["data"]][["x"]])[[2]]

      object
    })
  } else {
    object <- try(posterior_function(object))
  }

  # Produce something if the estimation fails, as the bvartools methods do: a
  # single unusable specification should not take a whole list of them down.
  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }

  class(object) <- class_of_object

  return(object)
}
