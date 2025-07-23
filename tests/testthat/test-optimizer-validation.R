# Test optimizer validation and setup

test_that("Optimizer validation works correctly", {
  setup <- setup_gaussian_mixture(d = 3)

  # Test unknown string optimizer
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      optimizer = "unknown_optimizer",
      gamma = setup$gamma
    ),
    "Unknown optimizer"
  )

  # Test invalid custom optimizer (missing required arguments)
  bad_optimizer <- function(x) { x^2 }
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      optimizer = bad_optimizer,
      gamma = setup$gamma
    ),
    "Custom optimizer must have signature"
  )

  # Test valid custom optimizer
  good_optimizer <- function(fn, par, lower, upper, ...) {
    stats::optim(par = par, fn = fn, method = "Nelder-Mead", ...)
  }

  expect_silent(
    result <- computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      optimizer = good_optimizer,
      profile_options = list(grid_points = 5),  # Small for speed
      verbose = FALSE,
      gamma = setup$gamma
    )
  )
  expect_true(is.list(result))
  expect_true("profiles" %in% names(result))
})

test_that("Bounds validation in validateAndSetupOptimizer works", {
  # Test invalid bounds structure
  expect_error(
    validateAndSetupOptimizer("optim", list(), bounds = list(lower = 1)),
    "bounds must be a list with 'lower' and 'upper' elements"
  )

  # Test non-numeric bounds
  expect_error(
    validateAndSetupOptimizer("optim", list(), bounds = list(lower = "a", upper = "b")),
    "bounds\\$lower and bounds\\$upper must be numeric vectors"
  )

  # Test mismatched lengths
  expect_error(
    validateAndSetupOptimizer("optim", list(), bounds = list(lower = c(1, 2), upper = c(1))),
    "bounds\\$lower and bounds\\$upper must have the same length"
  )

  # Test lower >= upper
  expect_error(
    validateAndSetupOptimizer("optim", list(), bounds = list(lower = c(2), upper = c(1))),
    "bounds\\$lower must be less than bounds\\$upper for all elements"
  )

  # Test automatic method selection with warning
  expect_warning(
    result <- validateAndSetupOptimizer("optim", list(method = "Nelder-Mead"),
                                       bounds = list(lower = c(0), upper = c(1))),
    "Method 'Nelder-Mead' specified for optim with bounds. Changing to 'L-BFGS-B'"
  )
  expect_equal(result$extra_args$method, "L-BFGS-B")

  # Test no warning when L-BFGS-B already specified
  expect_silent(
    result <- validateAndSetupOptimizer("optim", list(method = "L-BFGS-B"),
                                       bounds = list(lower = c(0), upper = c(1)))
  )
  expect_equal(result$extra_args$method, "L-BFGS-B")
})

# Test with DEoptim (if available)
test_that("DEoptim optimizer works if available", {
  if (requireNamespace("DEoptim", quietly = TRUE)) {
    setup <- setup_gaussian_mixture(d = 2)

    result <- computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      optimizer = "deoptim",
      optim_options = list(itermax = 20, NP = 10, trace = FALSE),  # Small for speed, disable trace
      profile_options = list(grid_points = 4),
      verbose = FALSE,
      gamma = setup$gamma
    )

    expect_true(is.list(result))
    expect_equal(result$optimizer, "deoptim")
    expect_equal(length(result$profiles), 2)
  } else {
    # Test that error is thrown when DEoptim is not available
    expect_error(
      validateAndSetupOptimizer("deoptim"),
      "DEoptim package required but not available"
    )
  }
})
