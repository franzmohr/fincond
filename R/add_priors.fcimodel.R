#' Add Priors to a Financial Conditions Index Model
#'
#' Adds the initial conditions and forgetting factors of the model of Koop and
#' Korobilis (2014), and the principal components the algorithm starts from.
#'
#' @param object an object of class \code{"fcimodel"}, usually the result of a
#' call to \code{\link{create_fcimodel}}.
#' @param lambda a named list of the prior on the loadings of the period before
#' the sample. Element \code{mu} is the mean and \code{v} the variance of every
#' element, both scalars.
#' @param beta a named list of the prior on the transition coefficients of the
#' period before the sample. If \code{mu} and \code{v} are \code{NULL}, a
#' Minnesota prior with shrinkage \code{kappa1} is obtained from
#' \code{\link[bvartools]{minnesota_prior}} on the state implied by the principal
#' components. Supplying \code{mu} and \code{v} as scalars replaces it with a
#' flat prior of that mean and variance.
#' @param factors a named list of the prior on the state of the period before the
#' sample, with elements \code{mu} and \code{v} as above.
#' @param u a named list with element \code{v}, the initial value of the diagonal
#' of the covariance matrix of the measurement errors.
#' @param v a named list with element \code{v}, the initial value of the diagonal
#' of the covariance matrix of the transition errors.
#' @param decay a named list of the four forgetting factors, each in \eqn{(0, 1]}.
#' See 'Details'.
#' @param ... not used.
#'
#' @details Only initial conditions get priors here. Nothing in this function is
#' a prior on the variance of an error term or on the variance of the drift in a
#' coefficient, because the model has none: those four quantities are obtained by
#' discounting instead, which is the whole reason the estimator needs no
#' simulation.
#'
#' The elements of \code{decay} are
#' \describe{
#'  \item{\code{u}}{the decay factor of the measurement error covariance
#'  \eqn{V_t}, which is the exponentially weighted moving average
#'  \eqn{V_t = \kappa_u V_{t - 1} + (1 - \kappa_u) e_t e_t'} of past squared
#'  prediction errors. Koop and Korobilis (2014) use \eqn{0.96}.}
#'  \item{\code{v}}{the same for the transition error covariance \eqn{Q_t}.
#'  \eqn{0.96}.}
#'  \item{\code{lambda}}{the forgetting factor of the loadings. The predicted
#'  state covariance is \eqn{R_t = S_{t - 1} / \kappa_\lambda}, so the drift of a
#'  loading is a fixed fraction of how uncertain that loading currently is rather
#'  than a variance of its own. \eqn{0.99}.}
#'  \item{\code{beta}}{the same for the transition coefficients. \eqn{0.99}.}
#' }
#'
#' A value of \code{1} switches the corresponding form of time variation off: the
#' covariance stops being updated, or the coefficients stop drifting.
#'
#' @return The model object, with the additional elements
#' \item{priors}{A list of the prior specifications.}
#' \item{data$factors}{A time-series object of the principal components the
#' algorithm starts from.}
#' and the forgetting factors in \code{model$decay}.
#'
#' @examples
#'
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
#' model <- add_priors(model)
#'
#' # Coefficients that do not drift at all
#' model_const <- add_priors(model, decay = list(u = 0.96, v = 0.96,
#'                                               lambda = 1, beta = 1))
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
#' @export
add_priors.fcimodel <- function(object,
                                lambda = list(mu = 0, v = 1),
                                beta = list(mu = NULL, v = NULL, kappa1 = 0.1),
                                factors = list(mu = 0, v = 10),
                                u = list(v = 0.1),
                                v = list(v = 0.1),
                                decay = list(u = 0.96, v = 0.96,
                                             lambda = 0.99, beta = 0.99),
                                ...) {

  n <- object[["model"]][["n"]]
  p <- object[["model"]][["p"]]
  n_y <- object[["model"]][["n_y"]]
  n_x <- object[["model"]][["n_x"]]
  constant <- object[["model"]][["constant"]]
  if (is.null(n) || is.null(p) || is.null(n_y)) {
    stop("Element 'model' is incomplete. Was the object produced by create_fcimodel?")
  }

  n_z <- n_y + n
  n_s <- n_z * p
  n_beta <- n_z * (n_z * p + as.integer(constant))

  # Forgetting factors ----
  object[["model"]][["decay"]] <- .fci_decay(decay)

  # Starting values ----
  #
  # The algorithm is a pair of conditional steps, so it has to be told where to
  # start. Principal components are what Koop and Korobilis start from, and they
  # are also what the Minnesota prior below is built on: a prior for the
  # transition equation needs a state to be a transition of.
  x <- object[["data"]][["x"]]
  y <- object[["data"]][["y"]]
  f <- .fci_principal_components(x, n)
  f <- stats::ts(f, start = stats::start(x), frequency = stats::frequency(x))
  dimnames(f)[[2]] <- paste0("factor", seq_len(n))
  object[["data"]][["factors"]] <- f

  # Loadings ----
  object[["priors"]][["lambda"]] <-
    list(mu = matrix(.fci_scalar(lambda, "mu", "lambda"), n_z),
         v = diag(.fci_positive(lambda, "v", "lambda"), n_z))

  # Transition coefficients ----
  if (is.null(beta[["mu"]]) && is.null(beta[["v"]])) {
    kappa1 <- if (is.null(beta[["kappa1"]])) 0.1 else beta[["kappa1"]]
    if (kappa1 <= 0) {
      stop("Argument 'beta$kappa1' must be positive.")
    }
    z <- stats::ts(cbind(y, f), start = stats::start(y),
                   frequency = stats::frequency(y))
    # The transition equation is a VAR in the state, so bvartools can produce
    # its Minnesota prior directly. The coefficient ordering the two packages
    # use is the same -- lags first, deterministic terms last, vectorised by
    # column -- which is what makes this reuse and not a coincidence.
    var_z <- bvartools::create_bvarmodel(
      z, p = p, deterministic = if (constant) "const" else "none")
    mn <- bvartools::minnesota_prior(var_z, kappa1 = kappa1)
    beta_mu <- matrix(mn[["mu"]], n_beta)
    beta_v <- solve(mn[["v_inv"]])
  } else {
    if (is.null(beta[["mu"]]) || is.null(beta[["v"]])) {
      stop("Arguments 'beta$mu' and 'beta$v' must either both be supplied or ",
           "both be NULL, in which case a Minnesota prior is used.")
    }
    beta_mu <- matrix(.fci_scalar(beta, "mu", "beta"), n_beta)
    beta_v <- diag(.fci_positive(beta, "v", "beta"), n_beta)
  }
  object[["priors"]][["beta"]] <- list(mu = beta_mu, v = beta_v)

  # State ----
  object[["priors"]][["factors"]] <-
    list(mu = matrix(.fci_scalar(factors, "mu", "factors"), n_s),
         v = diag(.fci_positive(factors, "v", "factors"), n_s))

  # Initial covariance matrices ----
  #
  # Only the diagonals: both covariance matrices are diagonal throughout, since
  # an exponentially weighted moving average of outer products is only used on
  # its diagonal here.
  object[["priors"]][["u"]] <- list(v = matrix(.fci_positive(u, "v", "u"), n_x))
  object[["priors"]][["v"]] <- list(v = matrix(.fci_positive(v, "v", "v"), n_z))

  return(object)
}

# The four forgetting factors, checked and put in the order the estimator reads
# them. Each has to be in (0, 1]: at one the corresponding quantity stops moving,
# and at zero the filter would forget everything it has seen.
.fci_decay <- function(decay) {
  required <- c("u", "v", "lambda", "beta")
  missing_fields <- setdiff(required, names(decay))
  if (length(missing_fields) > 0) {
    stop("Argument 'decay' is missing the forgetting factor",
         if (length(missing_fields) > 1) "s" else "", " ",
         paste0("'", missing_fields, "'", collapse = ", "), ".")
  }
  out <- list()
  for (field in required) {
    value <- decay[[field]]
    if (length(value) != 1 || !is.finite(value) || value <= 0 || value > 1) {
      stop("Argument 'decay$", field, "' must be a single value in (0, 1].")
    }
    out[[field]] <- as.numeric(value)
  }
  out
}

# One scalar out of a prior specification, with a message that names the
# argument it came from rather than the helper.
.fci_scalar <- function(spec, field, name) {
  value <- spec[[field]]
  if (is.null(value)) {
    stop("Argument '", name, "$", field, "' is missing.")
  }
  if (length(value) != 1 || !is.finite(value)) {
    stop("Argument '", name, "$", field, "' must be a single finite value.")
  }
  as.numeric(value)
}

.fci_positive <- function(spec, field, name) {
  value <- .fci_scalar(spec, field, name)
  if (value <= 0) {
    stop("Argument '", name, "$", field, "' must be larger than 0.")
  }
  value
}

# The first `n` principal components, scaled the way Koop and Korobilis scale
# them. The scale of a factor is arbitrary; this one is the starting value the
# filter runs from, so it is worth matching theirs.
.fci_principal_components <- function(x, n) {
  x <- as.matrix(x)
  n_x <- ncol(x)
  if (n > n_x) {
    stop("Argument 'n' cannot exceed the number of columns of 'x'.")
  }
  # Series that start late are held at their own mean here, which is zero once
  # the data are normalised. That imputation touches the starting value only:
  # the filter itself drops a missing observation from the period it is missing
  # in rather than making one up.
  for (j in seq_len(n_x)) {
    missing <- is.na(x[, j])
    if (any(missing)) {
      x[missing, j] <- mean(x[, j], na.rm = TRUE)
    }
  }
  load <- sqrt(n_x) * stats::prcomp(x)[["rotation"]][, seq_len(n), drop = FALSE]
  x %*% load / n_x
}
