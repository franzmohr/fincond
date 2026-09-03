#' Extract a Financial Conditions Index
#'
#' Extracts the estimated index from a model object and fixes its sign.
#'
#' @param object an object containing an estimated model.
#' @param ... arguments passed forward to method.
#'
#' @export
fci <- function(object, ...) {
  UseMethod("fci")
}

#' @param factor an integer of the factor to extract. Defaults to \code{1}, which
#' is the index of Koop and Korobilis (2014) in a model with \code{n = 1}.
#' @param sign_on the name or position of a column of \code{x} whose loading on
#' the factor should be positive. Defaults to \code{NULL}, which orients the
#' index so that the average loading across all series is positive. See
#' 'Details'.
#' @param standardize logical indicating whether the index should be rescaled to
#' zero mean and unity variance. Defaults to \code{TRUE}.
#'
#' @details The sign and scale of a factor are not identified: multiplying the
#' factor by \eqn{-1} and its loadings by \eqn{-1} leaves the model unchanged.
#' Both therefore have to be fixed by convention rather than estimated, and the
#' convention has to be stated, because it decides whether a rising index means
#' tightening or loosening financial conditions.
#'
#' The default orients the index so that the average loading over the financial
#' variables is positive: the index then moves with the typical series in
#' \code{x}. That is only meaningful when those series point the same way, so
#' with mixed data -- spreads rising as prices fall -- name a reference series
#' with \code{sign_on} and read the index against it.
#'
#' @return A time-series object.
#'
#' @examples
#'
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#'
#' # Oriented so that the index rises with the VIX, i.e. so that a high value
#' # means stressed conditions
#' index <- fci(model, sign_on = "vix")
#'
#' @references
#'
#' Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
#' \emph{European Economic Review, 71}, 101--116.
#' \doi{10.1016/j.euroecorev.2014.07.002}
#'
#' @rdname fci
#' @export
fci.fcimodel <- function(object, factor = 1, sign_on = NULL,
                         standardize = TRUE, ...) {

  .fci_require_posterior(object)

  n <- object[["model"]][["n"]]
  if (length(factor) != 1 || factor < 1 || factor > n) {
    stop("Argument 'factor' must be a single integer between 1 and ", n, ".")
  }

  index <- object[["posterior"]][["factors"]][, factor]
  load <- factor_loadings(object, factor = factor)

  # The average loading over the sample, one number per financial series. Which
  # of them decides the sign is what `sign_on` chooses.
  mean_load <- colMeans(load)
  reference <- if (is.null(sign_on)) {
    mean(mean_load)
  } else {
    pos <- .fci_column(sign_on, colnames(load), "sign_on")
    mean_load[pos]
  }
  if (is.finite(reference) && reference < 0) {
    index <- -index
  }

  if (standardize) {
    scale_by <- stats::sd(index)
    if (scale_by > 0) {
      index <- (index - mean(index)) / scale_by
    }
  }

  stats::ts(index, start = stats::start(object[["posterior"]][["factors"]]),
            frequency = stats::frequency(object[["posterior"]][["factors"]]))
}

#' Extract Time Varying Factor Loadings
#'
#' Extracts the estimated loadings of the financial variables on one factor.
#'
#' Named in full rather than \code{loadings}, which \pkg{stats} already uses for
#' something else -- and uses as a plain accessor rather than as a generic, so a
#' method could not have been added to it without masking it.
#'
#' @param x an object of class \code{"fcimodel"} containing an estimated model.
#' @param factor an integer of the factor whose loadings should be extracted.
#' @param ... not used.
#'
#' @details A loading is a series' exposure to the common component, and in this
#' model it moves. Reading the path of one is how a series is seen to enter or
#' leave the index -- an exposure that fell to zero over the sample says the
#' series stopped informing the index, which a constant-loading model would have
#' had to carry as noise throughout instead.
#'
#' @return A time-series object with one column per column of \code{x}.
#'
#' @examples
#'
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#'
#' head(factor_loadings(model))
#'
#' @export
factor_loadings <- function(x, factor = 1, ...) {

  object <- x
  .fci_require_posterior(object)

  n <- object[["model"]][["n"]]
  n_y <- object[["model"]][["n_y"]]
  n_x <- object[["model"]][["n_x"]]
  n_z <- n_y + n

  if (length(factor) != 1 || factor < 1 || factor > n) {
    stop("Argument 'factor' must be a single integer between 1 and ", n, ".")
  }

  # Within a period the loadings run along the rows of the n_x x n_z matrix, so
  # the column of series i on factor `factor` is at (i - 1) * n_z + n_y + factor.
  pos <- (seq_len(n_x) - 1) * n_z + n_y + factor
  out <- object[["posterior"]][["lambda"]][, pos, drop = FALSE]
  dimnames(out)[[2]] <- dimnames(object[["data"]][["x"]])[[2]]
  out
}

#' Plot a Financial Conditions Index Model
#'
#' @param x an object of class \code{"fcimodel"} containing an estimated model.
#' @param which what to plot: \code{"fci"} for the index, \code{"loadings"} for
#' the time varying loadings of the financial variables, or \code{"sigma"} for
#' the measurement error variances.
#' @param factor an integer of the factor to plot.
#' @param ... arguments passed forward to \code{\link[stats]{plot.ts}}.
#'
#' @return \code{NULL}, invisibly. Called for its side effect.
#'
#' @examples
#'
#' data("koop")
#'
#' model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
#' model <- add_priors(model)
#' model <- add_posterior_coefficients(model)
#'
#' plot(model)
#'
#' @export
plot.fcimodel <- function(x, which = "fci", factor = 1, ...) {

  which <- match.arg(which, c("fci", "loadings", "sigma"))
  .fci_require_posterior(x)

  if (which == "fci") {
    index <- fci(x, factor = factor, ...)
    stats::plot.ts(index, ylab = "Index", xlab = "",
                   main = "Financial conditions index")
    graphics::abline(h = 0, lty = "dashed", col = "grey40")
  }
  if (which == "loadings") {
    .fci_plot_many(factor_loadings(x, factor = factor), "Factor loadings", ...)
  }
  if (which == "sigma") {
    .fci_plot_many(x[["posterior"]][["u_sigma"]],
                   "Measurement error variances", ...)
  }

  invisible(NULL)
}

#' @param x an object of class \code{"fcimodel"}.
#' @param ... not used.
#'
#' @return \code{x}, invisibly.
#'
#' @rdname create_fcimodel
#' @export
print.fcimodel <- function(x, ...) {

  model <- x[["model"]]
  cat("Financial conditions index model\n\n")
  cat("Estimator:            Koop and Korobilis (2014) TVP-FAVAR\n")
  cat("Macro variables:     ", model[["n_y"]], "\n")
  cat("Financial variables: ", model[["n_x"]], "\n")
  cat("Factors:             ", model[["n"]], "\n")
  cat("Lag order:           ", model[["p"]], "\n")
  cat("Passes:              ", model[["iterations"]], "\n")
  if (!is.null(model[["decay"]])) {
    cat("Forgetting factors:   ",
        paste0(names(model[["decay"]]), " = ", unlist(model[["decay"]]),
               collapse = ", "), "\n", sep = "")
  }
  if (isTRUE(x[["error"]])) {
    cat("\nEstimation failed.\n")
  } else if (is.null(x[["posterior"]])) {
    cat("\nNot estimated yet.\n")
  } else {
    y <- x[["data"]][["y"]]
    cat("\nEstimated over", nrow(y), "periods,",
        paste0(stats::start(y)[1], " to ", stats::end(y)[1]), "\n")
    cat("Paths, not posterior draws. See ?fincond.\n")
  }

  invisible(x)
}

# One panel per series while that stays readable, and one panel with every series
# on it once it does not. plot.ts refuses more than ten panels outright, and an
# index is regularly built from more series than that.
.fci_plot_many <- function(series, main, ...) {
  if (NCOL(series) <= 10) {
    stats::plot.ts(series, main = main, xlab = "", ...)
  } else {
    colours <- grDevices::hcl.colors(NCOL(series), palette = "Dark 3")
    stats::plot.ts(series, plot.type = "single", col = colours,
                   main = main, xlab = "", ylab = "", ...)
    graphics::abline(h = 0, lty = "dashed", col = "grey40")
    graphics::legend("topleft", legend = colnames(series), col = colours,
                     lty = 1, bty = "n", cex = 0.6, ncol = 2)
  }
  invisible(NULL)
}

.fci_require_posterior <- function(object) {
  if (isTRUE(object[["error"]])) {
    stop("The estimation of this model failed.")
  }
  if (is.null(object[["posterior"]])) {
    stop("Element 'posterior' is missing. Did you call add_posterior_coefficients?")
  }
  invisible(NULL)
}

# A column position out of a name or a number, with a message that names the
# argument it came from.
.fci_column <- function(value, names_available, name) {
  if (is.character(value)) {
    pos <- match(value, names_available)
    if (is.na(pos)) {
      stop("Argument '", name, "' is not a column of 'x'.")
    }
    return(pos)
  }
  if (length(value) != 1 || value < 1 || value > length(names_available)) {
    stop("Argument '", name, "' must be a column name or a position in 'x'.")
  }
  as.integer(value)
}
