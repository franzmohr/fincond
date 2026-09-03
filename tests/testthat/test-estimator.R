make_model <- function(tt = 80, n_y = 2, n_x = 6, p = 2, n = 1, seed = 1) {
  set.seed(seed)
  y <- stats::ts(matrix(stats::rnorm(tt * n_y), tt, n_y),
                 start = c(1990, 1), frequency = 4)
  x <- stats::ts(matrix(stats::rnorm(tt * n_x), tt, n_x),
                 start = c(1990, 1), frequency = 4)
  create_fcimodel(y = y, x = x, p = p, n = n)
}

# A sample with one strong common factor in x, so that there is something for
# the estimator to find.
make_factor_model <- function(tt = 120, n_x = 8, noise = 0.3, seed = 7) {
  set.seed(seed)
  f <- stats::filter(stats::rnorm(tt), 0.8, method = "recursive")
  f <- as.numeric(f)
  load <- seq(0.8, 1.2, length.out = n_x)
  x <- outer(f, load) + matrix(stats::rnorm(tt * n_x, sd = noise), tt, n_x)
  y <- matrix(stats::rnorm(tt * 2), tt, 2)
  list(
    f = stats::ts(f, start = c(1990, 1), frequency = 4),
    model = create_fcimodel(
      y = stats::ts(y, start = c(1990, 1), frequency = 4),
      x = stats::ts(x, start = c(1990, 1), frequency = 4),
      p = 2, n = 1)
  )
}

test_that("the estimated paths have the shapes the model implies", {
  model <- add_posterior_coefficients(add_priors(make_model()))

  expect_false(isTRUE(model$error))

  tt <- nrow(model$data$y)
  n_z <- 2 + 1
  n_beta <- n_z * (n_z * 2 + 1)

  expect_equal(dim(model$posterior$factors), c(tt, 1))
  expect_equal(dim(model$posterior$lambda), c(tt, 6 * n_z))
  expect_equal(dim(model$posterior$beta), c(tt, n_beta))
  expect_equal(dim(model$posterior$u_sigma), c(tt, 6))
  expect_equal(dim(model$posterior$v_sigma), c(tt, n_z))
  expect_length(model$posterior$loglik, tt)

  for (i in c("factors", "lambda", "beta", "u_sigma", "v_sigma", "loglik")) {
    expect_true(stats::is.ts(model$posterior[[i]]), info = i)
    expect_true(all(is.finite(model$posterior[[i]])), info = i)
    expect_equal(stats::tsp(model$posterior[[i]]), stats::tsp(model$data$y),
                 info = i)
  }
})

test_that("the estimator is deterministic", {
  model <- add_priors(make_model())
  first <- add_posterior_coefficients(model)
  second <- add_posterior_coefficients(model)

  expect_identical(as.numeric(first$posterior$factors),
                   as.numeric(second$posterior$factors))
  expect_identical(as.numeric(first$posterior$beta),
                   as.numeric(second$posterior$beta))
})

test_that("a loading forgetting factor of one gives a constant loading path", {
  # With no discounting the smoother gain is the identity, so every smoothed
  # value equals the last filtered one and the path does not drift. This is the
  # exact degeneracy of the model, not an approximation, so it is checked
  # exactly.
  model <- add_priors(make_model(),
                      decay = list(u = 0.96, v = 0.96, lambda = 1, beta = 0.99))
  model <- add_posterior_coefficients(model)

  path <- as.matrix(model$posterior$lambda)
  expect_equal(max(abs(sweep(path, 2, path[1, ]))), 0)
})

test_that("a transition forgetting factor of one gives constant coefficients", {
  model <- add_priors(make_model(),
                      decay = list(u = 0.96, v = 0.96, lambda = 0.99, beta = 1))
  model <- add_posterior_coefficients(model)

  path <- as.matrix(model$posterior$beta)
  expect_equal(max(abs(sweep(path, 2, path[1, ]))), 0)
})

test_that("an error decay of one leaves the covariances at their prior", {
  model <- add_priors(make_model(),
                      decay = list(u = 1, v = 1, lambda = 0.99, beta = 0.99),
                      u = list(v = 0.25), v = list(v = 0.5))
  model <- add_posterior_coefficients(model)

  expect_equal(as.numeric(model$posterior$u_sigma),
               rep(0.25, nrow(model$data$y) * 6))
  expect_equal(as.numeric(model$posterior$v_sigma),
               rep(0.5, nrow(model$data$y) * 3))
})

test_that("the variances move when they are allowed to", {
  model <- add_posterior_coefficients(add_priors(make_model()))

  expect_gt(max(apply(model$posterior$u_sigma, 2, stats::sd)), 0)
  expect_gt(max(apply(model$posterior$v_sigma, 2, stats::sd)), 0)
})

test_that("a common factor in x is recovered", {
  sim <- make_factor_model()
  model <- add_posterior_coefficients(add_priors(sim$model))

  index <- fci(model, standardize = TRUE)
  # The sign of a factor is a convention, so the size of the correlation is what
  # is being checked here, not its direction.
  expect_gt(abs(stats::cor(as.numeric(index), as.numeric(sim$f))), 0.9)
})

test_that("more passes change the estimate but keep it sane", {
  sim <- make_factor_model()
  one <- add_posterior_coefficients(add_priors(sim$model))

  model_five <- sim$model
  model_five$model$iterations <- 5L
  five <- add_posterior_coefficients(add_priors(model_five))

  expect_true(all(is.finite(five$posterior$factors)))
  expect_gt(abs(stats::cor(as.numeric(fci(five)), as.numeric(sim$f))), 0.9)
  expect_false(identical(as.numeric(one$posterior$factors),
                         as.numeric(five$posterior$factors)))
})

test_that("estimation without priors fails softly", {
  model <- suppressWarnings(add_posterior_coefficients(make_model()))
  expect_true(isTRUE(model$error))
  expect_s3_class(model, "fcimodel")
})

test_that("a series that starts late enters the index when it starts", {
  sim <- make_factor_model()

  # The same sample with one series unobserved for its first half.
  gapped <- sim$model
  tt <- nrow(gapped$data$x)
  gapped$data$x[seq_len(tt %/% 2), 1] <- NA

  full <- add_posterior_coefficients(add_priors(sim$model))
  part <- add_posterior_coefficients(add_priors(gapped))

  expect_false(isTRUE(part$error))
  expect_true(all(is.finite(part$posterior$factors)))
  expect_true(all(is.finite(part$posterior$lambda)))

  # Over the gap the filter has nothing to update the loading with, so every
  # value there is the same prediction. The smoothed path is therefore the
  # backward blend of that one value towards the first estimate the series does
  # inform, which is monotone -- no observation before the series starts can
  # push it either way.
  load <- as.numeric(factor_loadings(part)[, 1])
  gap <- seq_len(tt %/% 2)
  seen <- (tt %/% 2 + 1):tt
  expect_true(all(diff(load[gap]) >= 0) || all(diff(load[gap]) <= 0))

  # Over the gap it travels from the prior towards the level the data imply,
  # and once the series is observed it settles at roughly what the model with
  # no gap finds for it.
  expect_gt(abs(load[max(gap)] - load[1]), abs(diff(range(load[seen]))))
  expect_equal(mean(load[seen]),
               mean(as.numeric(factor_loadings(full)[seen, 1])),
               tolerance = 0.15)

  # The index is still the index: dropping half of one of eight series should
  # not change what it is measuring.
  expect_gt(abs(stats::cor(as.numeric(fci(part)), as.numeric(fci(full)))), 0.95)
})

test_that("the predictive likelihood ignores the periods a series is missing", {
  sim <- make_factor_model()
  gapped <- sim$model
  gapped$data$x[1:20, 1] <- NA
  part <- add_posterior_coefficients(add_priors(gapped))

  expect_true(all(is.finite(part$posterior$loglik)))
})
