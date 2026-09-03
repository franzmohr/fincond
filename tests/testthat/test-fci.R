make_estimated <- function(tt = 100, n_x = 6, seed = 3) {
  set.seed(seed)
  f <- as.numeric(stats::filter(stats::rnorm(tt), 0.8, method = "recursive"))
  x <- outer(f, seq(0.8, 1.2, length.out = n_x)) +
    matrix(stats::rnorm(tt * n_x, sd = 0.3), tt, n_x)
  dimnames(x)[[2]] <- paste0("f", seq_len(n_x))
  y <- matrix(stats::rnorm(tt * 2), tt, 2)
  model <- create_fcimodel(
    y = stats::ts(y, start = c(1990, 1), frequency = 4),
    x = stats::ts(x, start = c(1990, 1), frequency = 4),
    p = 2, n = 1)
  add_posterior_coefficients(add_priors(model))
}

test_that("the index is a standardised time series over the sample", {
  model <- make_estimated()
  index <- fci(model)

  expect_true(stats::is.ts(index))
  expect_length(index, nrow(model$data$y))
  expect_equal(mean(index), 0)
  expect_equal(stats::sd(index), 1)
  expect_equal(stats::tsp(index), stats::tsp(model$data$y))
})

test_that("standardize = FALSE returns the factor on its own scale", {
  model <- make_estimated()
  raw <- fci(model, standardize = FALSE)
  expect_equal(abs(as.numeric(raw)),
               abs(as.numeric(model$posterior$factors[, 1])))
})

test_that("the sign convention orients the index on the chosen series", {
  model <- make_estimated()

  # Every series loads positively here, so pinning the sign to any one of them
  # gives the same orientation as the default.
  default <- fci(model)
  on_f1 <- fci(model, sign_on = "f1")
  expect_equal(as.numeric(default), as.numeric(on_f1))

  # A series is named or given by position, and both must agree.
  expect_equal(as.numeric(fci(model, sign_on = 2)),
               as.numeric(fci(model, sign_on = "f2")))

  expect_error(fci(model, sign_on = "nope"), "not a column")
  expect_error(fci(model, sign_on = 99), "column name or a position")
})

test_that("flipping every loading flips the index", {
  model <- make_estimated()
  before <- fci(model, standardize = FALSE)

  flipped <- model
  flipped$posterior$lambda <- -flipped$posterior$lambda
  after <- fci(flipped, standardize = FALSE)

  expect_equal(as.numeric(after), -as.numeric(before))
})

test_that("the loadings come back one column per financial series", {
  model <- make_estimated()
  load <- factor_loadings(model)

  expect_true(stats::is.ts(load))
  expect_equal(dim(load), c(nrow(model$data$y), 6))
  expect_equal(colnames(load), paste0("f", 1:6))
  # A strong common factor with positive exposures should be found as such.
  expect_true(all(colMeans(load) > 0) || all(colMeans(load) < 0))
})

test_that("a factor that does not exist is refused", {
  model <- make_estimated()
  expect_error(fci(model, factor = 2), "between 1 and 1")
  expect_error(factor_loadings(model, factor = 0), "between 1 and 1")
})

test_that("an unestimated model is refused rather than returning nothing", {
  set.seed(1)
  y <- stats::ts(matrix(stats::rnorm(120), 60, 2), start = c(1990, 1), frequency = 4)
  x <- stats::ts(matrix(stats::rnorm(300), 60, 5), start = c(1990, 1), frequency = 4)
  model <- create_fcimodel(y = y, x = x, p = 2, n = 1)

  expect_error(fci(model), "Did you call add_posterior_coefficients")
  expect_error(factor_loadings(model), "Did you call add_posterior_coefficients")
})

test_that("print says what the model is and whether it ran", {
  model <- make_estimated()
  expect_output(print(model), "Financial conditions index model")
  expect_output(print(model), "Paths, not posterior draws")
})

test_that("the deprecated entry point points at the new one", {
  expect_error(kk_algorithm(1, 2), "create_fcimodel")
})
