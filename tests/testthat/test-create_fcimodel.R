make_data <- function(tt = 60, n_y = 2, n_x = 5, seed = 1) {
  set.seed(seed)
  y <- stats::ts(matrix(stats::rnorm(tt * n_y), tt, n_y),
                 start = c(1990, 1), frequency = 4)
  x <- stats::ts(matrix(stats::rnorm(tt * n_x), tt, n_x),
                 start = c(1990, 1), frequency = 4)
  dimnames(y)[[2]] <- paste0("m", seq_len(n_y))
  dimnames(x)[[2]] <- paste0("f", seq_len(n_x))
  list(y = y, x = x)
}

test_that("a single specification gives one fcimodel", {
  d <- make_data()
  model <- create_fcimodel(y = d$y, x = d$x, p = 2, n = 1)

  expect_s3_class(model, "fcimodel")
  expect_equal(model$model$n_y, 2)
  expect_equal(model$model$n_x, 5)
  expect_equal(model$model$n, 1)
  expect_equal(model$model$p, 2)
  expect_equal(model$model$algorithm, "FavarKk2014")
})

test_that("vectors of p and n give a modellist of every combination", {
  d <- make_data()
  models <- create_fcimodel(y = d$y, x = d$x, p = 1:3, n = 1:2)

  expect_s3_class(models, "modellist")
  expect_length(models, 6)
  expect_true(all(vapply(models, inherits, logical(1), "fcimodel")))
})

test_that("x is normalised by default and left alone when asked", {
  d <- make_data()
  model <- create_fcimodel(y = d$y, x = d$x, p = 2, n = 1)
  expect_equal(unname(colMeans(model$data$x)), rep(0, 5))
  expect_equal(unname(apply(model$data$x, 2, stats::sd)), rep(1, 5))

  raw <- create_fcimodel(y = d$y, x = d$x, p = 2, n = 1, normalize_x = FALSE)
  expect_equal(unname(as.matrix(raw$data$x)), unname(as.matrix(d$x)))
})

test_that("y and x are cut to the periods they share", {
  d <- make_data(tt = 60)
  # x starts two years late and ends a year early
  x_short <- stats::window(d$x, start = c(1992, 1), end = c(2003, 4))
  model <- create_fcimodel(y = d$y, x = x_short, p = 2, n = 1)

  expect_equal(stats::start(model$data$y), c(1992, 1))
  expect_equal(stats::end(model$data$y), stats::end(x_short))
  expect_equal(nrow(model$data$y), nrow(model$data$x))
})

test_that("column names survive", {
  d <- make_data()
  model <- create_fcimodel(y = d$y, x = d$x, p = 2, n = 1)
  expect_equal(colnames(model$data$x), paste0("f", 1:5))
  expect_equal(colnames(model$data$y), paste0("m", 1:2))
})

test_that("impossible inputs are refused", {
  d <- make_data()

  expect_error(create_fcimodel(y = unclass(d$y), x = d$x), "class 'ts'")
  expect_error(create_fcimodel(y = d$y, x = unclass(d$x)), "class 'ts'")
  expect_error(create_fcimodel(y = d$y, x = d$x, p = 0), "at least 1")
  expect_error(create_fcimodel(y = d$y, x = d$x, n = 0), "at least 1")
  expect_error(create_fcimodel(y = d$y, x = d$x, iterations = 0), "at least 1")
  expect_error(create_fcimodel(y = d$y, x = d$x, p = 500), "Too few periods")

  # No overlap at all
  x_late <- stats::ts(d$x, start = c(2100, 1), frequency = 4)
  expect_error(create_fcimodel(y = d$y, x = x_late), "no periods in common")
})

test_that("a gap in the macro block is refused", {
  d <- make_data()
  d$y[5, 2] <- NA
  expect_error(create_fcimodel(y = d$y, x = d$x, p = 2, n = 1),
               "'y' must not contain missing values")
})

test_that("a gap in the financial block is carried and counted", {
  d <- make_data()
  d$x[1:10, 2] <- NA
  model <- create_fcimodel(y = d$y, x = d$x, p = 2, n = 1)

  expect_equal(model$model$n_missing, 10)
  expect_equal(sum(is.na(model$data$x)), 10)
  # Normalisation uses the observations a series has, not the ones it has not.
  expect_equal(mean(model$data$x[, 2], na.rm = TRUE), 0)
  expect_equal(stats::sd(model$data$x[, 2], na.rm = TRUE), 1)
})

test_that("a financial series with nothing in the shared sample is named", {
  d <- make_data()
  d$x[, 3] <- NA
  expect_error(create_fcimodel(y = d$y, x = d$x, p = 2, n = 1), "'f3'")
})
