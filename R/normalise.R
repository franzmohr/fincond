#' Normalise Series
#'
#' Normalises the columns of a time-series object to zero mean and unity
#' variance.
#'
#' @param x a time-series object.
#'
#' @details Differs from \code{\link[base]{scale}} in two ways that matter for a
#' factor model: the time-series attributes and the column names survive, and
#' each column is normalised over its own non-missing observations rather than
#' the whole matrix being dropped where any series is short.
#'
#' A factor extracted from series of different scales is a factor of the largest
#' of them, which is why \code{\link{create_fcimodel}} does this by default.
#'
#' @return A time-series object of the same shape as \code{x}.
#'
#' @examples
#'
#' data("koop")
#'
#' x <- normalise(koop$financial)
#' round(colMeans(x, na.rm = TRUE), 10)
#'
#' @export
normalise <- function(x) {
  tsp_temp <- stats::tsp(x)
  names_temp <- dimnames(as.matrix(x))[[2]]
  x <- apply(as.matrix(x), 2, .norm_helper)
  x <- stats::as.ts(x)
  stats::tsp(x) <- tsp_temp
  dimnames(x)[[2]] <- names_temp
  return(x)
}

# One column, normalised over the observations it has. A constant series has no
# scale to divide by, so it is centred and left alone rather than turned into
# NaN, which would take the whole factor with it.
.norm_helper <- function(x) {
  pos <- which(!is.na(x))
  result <- x[pos]
  scale_by <- stats::sd(result)
  if (!is.finite(scale_by) || scale_by == 0) {
    x[pos] <- result - mean(result)
  } else {
    x[pos] <- (result - mean(result)) / scale_by
  }
  return(x)
}
