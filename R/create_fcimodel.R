#' Create a Financial Conditions Index Model
#'
#' Produces the input for the estimation of the time varying parameter FAVAR of
#' Koop and Korobilis (2014).
#'
#' @param y a time-series object of macroeconomic variables. These enter the
#' model as observed elements of the state vector and are measured without
#' error.
#' @param x a time-series object of financial variables, from which the index is
#' extracted.
#' @param p an integer vector of the lag order of the transition equation. See
#' 'Details'.
#' @param n an integer vector of the number of unobserved factors. See 'Details'.
#' The index of Koop and Korobilis (2014) is the single factor of \code{n = 1}.
#' @param constant logical indicating whether the transition equation contains an
#' intercept. Defaults to \code{TRUE}.
#' @param normalize_x logical indicating whether each column of \code{x} should
#' be normalised to zero mean and unity variance. Defaults to \code{TRUE}. A
#' factor extracted from series of different scales is a factor of the largest of
#' them, so this should only be switched off for data that are already
#' comparable.
#' @param iterations an integer of passes of the two-step algorithm. Defaults to
#' \code{1}, which is the estimator as published. This is \emph{not} a number of
#' MCMC draws: the model is not simulated. See 'Details'.
#'
#' @details The function produces the variable matrices of a factor augmented
#' VAR with time varying parameters. Let \eqn{y_t} be the \eqn{n_y \times 1}
#' vector of macroeconomic variables and \eqn{x_t} the \eqn{n_x \times 1} vector
#' of financial variables. The state is \eqn{z_t = (y_t', f_t')'}, where
#' \eqn{f_t} is the \eqn{n \times 1} vector of unobserved factors, and the
#' measurement equation is
#' \deqn{\left[\begin{array}{c} y_t \\ x_t \end{array}\right] =
#' \left[\begin{array}{cc} I & 0 \\ \lambda^y_t & \lambda^f_t \end{array}\right]
#' \left[\begin{array}{c} y_t \\ f_t \end{array}\right] +
#' \left[\begin{array}{c} 0 \\ e_t \end{array}\right],}
#' with \eqn{e_t \sim N(0, V_t)} and \eqn{V_t} diagonal. The macroeconomic block
#' is measured without error, which is what makes the corresponding block of the
#' state equal to the data.
#'
#' The transition equation is
#' \deqn{z_t = \sum_{i = 1}^{p} B_{i, t} z_{t - i} + u_t,}
#' with \eqn{u_t \sim N(0, Q_t)}.
#'
#' Both coefficient blocks follow random walks, \eqn{\lambda_t = \lambda_{t - 1} +
#' \epsilon_t} and \eqn{\beta_t = \beta_{t - 1} + \eta_t}. The variances of those
#' innovations are not estimated. They, and the two error covariances, are
#' obtained from four forgetting factors set in \code{\link{add_priors.fcimodel}}
#' -- see there for what each of them does. This is the substantive difference
#' from a fully Bayesian time varying dynamic factor model, and it is why the
#' estimator needs no simulation.
#'
#' Argument \code{iterations} repeats the two steps of the algorithm -- the
#' parameters given the factors, then the factors given the parameters. Koop and
#' Korobilis (2014) run each exactly once, starting from principal components,
#' which is the default here. Values above one iterate towards a fixed point of
#' the two conditional problems. They do not produce draws and do not need a
#' burn-in: nothing in this package is simulated.
#'
#' The two blocks are cut to the periods they share. Missing values are then
#' treated asymmetrically, because the model treats the two blocks
#' asymmetrically: a gap in \code{y} is refused, since that block of the state is
#' the data and is measured without error, while a gap in \code{x} is the
#' ordinary case -- the series of a financial index rarely all start on the same
#' date -- and is carried through. The filter drops a series from the periods it
#' is missing in rather than imputing anything, so a series that starts halfway
#' through contributes to the index from the date it starts and not before. This
#' differs from the code published with the paper, which replaces missing values
#' with the sample mean throughout. The one place a mean is still substituted is
#' the principal components that the algorithm starts from, which need a complete
#' matrix and only set a starting value.
#'
#' If integer vectors are provided as arguments \code{p} or \code{n}, the
#' function produces a distinct model for all combinations of those
#' specifications, in a list of class \code{"modellist"}.
#'
#' @return An object of class \code{"fcimodel"}, containing the elements
#' \item{data}{A list with the time-series objects \code{y} and \code{x}, the
#' latter normalised if \code{normalize_x} was \code{TRUE}.}
#' \item{model}{A list of model specifications.}
#'
#' @examples
#'
#' # Load the data of Koop and Korobilis (2014)
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial,
#'                          p = 4, n = 1)
#'
#' # Several lag orders at once
#' models <- create_fcimodel(y = koop$macro[[96]], x = koop$financial,
#'                           p = 1:4, n = 1)
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
#' @export
create_fcimodel <- function(y, x, p = 4, n = 1, constant = TRUE,
                            normalize_x = TRUE, iterations = 1) {

  # Input checks ----
  if (!"ts" %in% class(y)) {
    stop("Argument 'y' must be an object of class 'ts'.")
  }
  if (!"ts" %in% class(x)) {
    stop("Argument 'x' must be an object of class 'ts'.")
  }
  if (any(p < 1)) {
    stop("Argument 'p' must be at least 1.")
  }
  if (any(n < 1)) {
    stop("Argument 'n' must be at least 1.")
  }
  if (!"logical" %in% class(constant)) {
    stop("Argument 'constant' must be of class 'logical'.")
  }
  if (!"logical" %in% class(normalize_x)) {
    stop("Argument 'normalize_x' must be of class 'logical'.")
  }
  if (length(iterations) != 1 || iterations < 1) {
    stop("Argument 'iterations' must be a single integer of at least 1.")
  }

  # Data preparation ----

  y <- .as_named_mts(y, "y")
  x <- .as_named_mts(x, "x")

  # The two blocks must cover the same periods, and the model has no way of
  # saying which observation belongs to which date once they are matrices.
  # Disjoint series come back as NULL with a warning, which is not something to
  # let through to a reshape further down.
  common <- suppressWarnings(stats::ts.intersect(y, x))
  if (is.null(common) || nrow(common) == 0) {
    stop("Arguments 'y' and 'x' have no periods in common.")
  }
  n_y <- NCOL(y)
  n_x <- NCOL(x)
  tsp_common <- stats::tsp(common)
  names_y <- dimnames(y)[[2]]
  names_x <- dimnames(x)[[2]]

  y <- stats::ts(common[, seq_len(n_y), drop = FALSE],
                 start = tsp_common[1], frequency = tsp_common[3])
  x <- stats::ts(common[, n_y + seq_len(n_x), drop = FALSE],
                 start = tsp_common[1], frequency = tsp_common[3])
  dimnames(y)[[2]] <- names_y
  dimnames(x)[[2]] <- names_x

  # A missing macro observation is fatal: that block of the state is the data,
  # measured without error, so there is nothing for the filter to fall back on.
  # A missing financial observation is not -- it is the ordinary case, since the
  # series of an index rarely all start on the same date -- and the filter skips
  # it rather than imputing anything. See 'Details'.
  if (anyNA(y)) {
    stop("Argument 'y' must not contain missing values after alignment. ",
         "Use stats::na.omit or impute before calling.")
  }
  n_missing <- sum(is.na(x))
  if (all(is.na(x))) {
    stop("Argument 'x' contains no observations.")
  }
  short <- colSums(!is.na(x)) < 2
  if (any(short)) {
    stop("Every column of 'x' needs at least two observations over the periods ",
         "it shares with 'y'. Too short: ",
         paste0("'", dimnames(x)[[2]][short], "'", collapse = ", "),
         ". Drop them, or widen the sample 'y' covers.")
  }

  if (normalize_x) {
    x <- normalise(x)
  }

  tt <- nrow(y)
  if (tt <= max(p) + 1) {
    stop("Too few periods for the requested lag order.")
  }

  model <- NULL
  model$type <- "FCI"
  # The estimator add_posterior_coefficients dispatches on. Named rather than
  # derived at the point of use, so that the object says what produced it.
  model$algorithm <- "FavarKk2014"
  model$n_y <- n_y
  model$n_x <- n_x
  model$n_missing <- n_missing
  model$n <- 0
  model$p <- 0
  model$constant <- constant
  model$iterations <- as.integer(iterations)

  result <- NULL
  for (j in n) {
    for (i in p) {
      model_i <- model
      model_i$n <- as.integer(j)
      model_i$p <- as.integer(i)

      result_i <- list("data" = list("y" = y, "x" = x),
                       "model" = model_i)

      class(result_i) <- append("fcimodel", class(result_i))

      result <- c(result, list(result_i))
    }
  }

  if (length(result) == 1) {
    result <- result[[1]]
  } else {
    class(result) <- append("modellist", class(result))
  }

  return(result)
}

# A time-series object as a named matrix, so that every downstream label has
# something to come from. `name` is what a message calls the argument.
.as_named_mts <- function(x, name) {
  tsp_temp <- stats::tsp(x)
  out <- stats::ts(as.matrix(x), start = tsp_temp[1], frequency = tsp_temp[3])
  if (is.null(dimnames(out)[[2]])) {
    dimnames(out)[[2]] <- if (NCOL(out) == 1) {
      name
    } else {
      paste0(name, seq_len(NCOL(out)))
    }
  }
  out
}
