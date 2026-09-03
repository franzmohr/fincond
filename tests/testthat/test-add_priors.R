make_model <- function(tt = 60, n_y = 2, n_x = 5, p = 2, n = 1, seed = 1) {
  set.seed(seed)
  y <- stats::ts(matrix(stats::rnorm(tt * n_y), tt, n_y),
                 start = c(1990, 1), frequency = 4)
  x <- stats::ts(matrix(stats::rnorm(tt * n_x), tt, n_x),
                 start = c(1990, 1), frequency = 4)
  create_fcimodel(y = y, x = x, p = p, n = n)
}

test_that("the priors have the widths the model implies", {
  model <- add_priors(make_model(n_y = 2, n_x = 5, p = 2, n = 1))

  n_z <- 2 + 1
  n_beta <- n_z * (n_z * 2 + 1)

  expect_equal(dim(model$priors$lambda$mu), c(n_z, 1))
  expect_equal(dim(model$priors$lambda$v), c(n_z, n_z))
  expect_equal(dim(model$priors$beta$mu), c(n_beta, 1))
  expect_equal(dim(model$priors$beta$v), c(n_beta, n_beta))
  expect_equal(dim(model$priors$factors$mu), c(n_z * 2, 1))
  expect_equal(dim(model$priors$factors$v), c(n_z * 2, n_z * 2))
  expect_equal(nrow(model$priors$u$v), 5)
  expect_equal(nrow(model$priors$v$v), n_z)
})

test_that("the starting values are the principal components", {
  model <- add_priors(make_model())

  expect_true(!is.null(model$data$factors))
  expect_equal(nrow(model$data$factors), nrow(model$data$x))
  expect_equal(ncol(model$data$factors), 1)
  expect_true(stats::is.ts(model$data$factors))

  # Scaled the way Koop and Korobilis scale them: sqrt(n_x) times the
  # eigenvector, divided by n_x.
  x <- as.matrix(model$data$x)
  n_x <- ncol(x)
  expected <- x %*% (sqrt(n_x) * stats::prcomp(x)$rotation[, 1, drop = FALSE]) / n_x
  expect_equal(as.numeric(model$data$factors), as.numeric(expected))
})

test_that("the transition prior is the bvartools Minnesota prior by default", {
  model <- make_model(p = 2)
  model <- add_priors(model)

  # The same prior, built the long way round through bvartools on the state.
  f <- model$data$factors
  z <- stats::ts(cbind(model$data$y, f), start = stats::start(model$data$y),
                 frequency = stats::frequency(model$data$y))
  var_z <- bvartools::create_bvarmodel(z, p = 2, deterministic = "const")
  mn <- bvartools::minnesota_prior(var_z, kappa1 = 0.1)

  expect_equal(unname(model$priors$beta$mu), unname(matrix(mn$mu)))
  expect_equal(unname(model$priors$beta$v), unname(solve(mn$v_inv)))
})

test_that("an explicit transition prior replaces the Minnesota one", {
  model <- add_priors(make_model(), beta = list(mu = 0, v = 9))
  expect_true(all(model$priors$beta$mu == 0))
  expect_equal(unique(diag(model$priors$beta$v)), 9)
})

test_that("the forgetting factors default to the paper's values", {
  model <- add_priors(make_model())
  expect_equal(model$model$decay,
               list(u = 0.96, v = 0.96, lambda = 0.99, beta = 0.99))
})

test_that("impossible forgetting factors are refused", {
  model <- make_model()

  for (bad in list(0, -0.5, 1.1, NA_real_, c(0.9, 0.9))) {
    expect_error(
      add_priors(model, decay = list(u = bad, v = 0.96, lambda = 0.99, beta = 0.99)),
      "decay\\$u"
    )
  }

  expect_error(add_priors(model, decay = list(u = 0.96, v = 0.96)),
               "missing the forgetting factors")
})

test_that("impossible variances are refused", {
  model <- make_model()

  expect_error(add_priors(model, lambda = list(mu = 0, v = 0)), "larger than 0")
  expect_error(add_priors(model, u = list(v = -1)), "larger than 0")
  expect_error(add_priors(model, factors = list(mu = 0)), "factors\\$v' is missing")
  expect_error(add_priors(model, beta = list(mu = 0, v = NULL)),
               "both be supplied")
})
