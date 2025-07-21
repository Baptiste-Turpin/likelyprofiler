# Test Suite for likelyprofiler Package
# This file contains comprehensive tests for the profile likelihood computation package

library(testthat)
library(likelyprofiler)

# Test helper functions
test_that("Utility functions work correctly", {
  # Test null coalescing operator
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(5 %||% 10, 5)
})

# Test data setup
setup_gaussian_mixture <- function(d = 5, gamma = 8) {
  # Define the Gaussian mixture negative log-likelihood
  mlogf <- function(x, gamma) {
    m1 <- m2 <- rep(0, length(x))
    m1[1] <- 1
    m2[1] <- 0
    -log(exp(-gamma * sum((x - m1)^2)) + exp(-gamma * sum((x - m2)^2)))
  }
  
  # Problem parameters
  true_params <- c(1, rep(0, d - 1))
  set.seed(23)
  optimized_params <- true_params + rnorm(length(true_params), sd = 0.05)
  names(true_params) <- names(optimized_params) <- paste0("param_", 1:d)
  bounds <- list(lower = rep(-2, d), upper = rep(2, d))
  
  list(
    mlogf = mlogf,
    gamma = gamma,
    true_params = true_params,
    optimized_params = optimized_params,
    bounds = bounds,
    d = d
  )
}

setup_simple_quadratic <- function() {
  # Simple quadratic function for testing
  quadratic_likelihood <- function(params) {
    x <- params[1]
    y <- params[2]
    (x - 2)^2 + 2 * (y - 1)^2
  }
  
  true_params <- c(2.0, 1.0)
  names(true_params) <- c("param_x", "param_y")
  bounds <- list(lower = c(0, 0), upper = c(4, 3))
  
  list(
    likelihood = quadratic_likelihood,
    true_params = true_params,
    bounds = bounds
  )
}

# Test basic input validation
test_that("Input validation works correctly", {
  setup <- setup_gaussian_mixture(d = 2)
  
  # Test invalid params_current
  expect_error(
    computeLikelihoodProfiles(
      params_current = "invalid",
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      gamma = setup$gamma
    ),
    "params_current must be a numeric vector"
  )
  
  # Test invalid negLogLikelihood
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = "not_a_function",
      bounds = setup$bounds,
      gamma = setup$gamma
    ),
    "negLogLikelihood must be a function"
  )
  
  # Test invalid bounds
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = list(lower = c(1, 2)),  # missing upper
      gamma = setup$gamma
    ),
    "bounds must be a list with 'lower' and 'upper' elements"
  )
  
  # Test mismatched bounds dimensions
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = list(lower = c(1), upper = c(1, 2)),  # different lengths
      gamma = setup$gamma
    ),
    "bounds dimensions must match parameter vector length"
  )
})

# Test optimizer validation
test_that("Optimizer validation works correctly", {
  setup <- setup_gaussian_mixture(d = 2)
  
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

# Test bounds validation in validateAndSetupOptimizer
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

# Test multidimensional Gaussian mixture (as in vignette)
test_that("Multidimensional Gaussian mixture profiling works", {
  setup <- setup_gaussian_mixture(d = 5)
  
  result <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 10,  # Small for speed
      max_grid_range_multiplier = 1.0
    ),
    optimizer = "optim",
    optim_options = list(
      method = "L-BFGS-B",
      control = list(maxit = 50)
    ),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Check result structure
  expect_true(is.list(result))
  expect_true(all(c("profiles", "confidence_intervals", "summary", "options", "optimizer") %in% names(result)))
  expect_equal(length(result$profiles), setup$d)
  expect_equal(nrow(result$confidence_intervals), setup$d)
  expect_equal(nrow(result$summary), setup$d)
  expect_equal(result$optimizer, "optim")
  
  # Check confidence intervals structure
  expect_equal(ncol(result$confidence_intervals), 2)
  expect_equal(colnames(result$confidence_intervals), c("2.5%", "97.5%"))
  expect_equal(rownames(result$confidence_intervals), names(setup$optimized_params))
  
  # Check summary structure
  expect_true(all(c("parameter", "current_value", "lower_ci", "upper_ci", "ci_width", "status") %in% names(result$summary)))
  
  # Check that at least some profiles succeeded
  expect_true(any(result$summary$status == "success"))
})

# Test one-dimensional case
test_that("One-dimensional Gaussian mixture profiling works", {
  setup <- setup_gaussian_mixture(d = 1)
  
  result <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "uniform",
      grid_points = 8,  # Small for speed
      max_grid_range_multiplier = 1.5
    ),
    optimizer = "optim",
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Check result structure for 1D case
  expect_equal(length(result$profiles), 1)
  expect_equal(nrow(result$confidence_intervals), 1)
  expect_equal(nrow(result$summary), 1)
  
  # Check that single parameter profile worked
  profile <- result$profiles[[1]]
  expect_false(profile$failed)
  expect_true(length(profile$grid_values) >= 3)
  expect_true(length(profile$profile_costs) >= 3)
})

# Test simple quadratic function
test_that("Simple quadratic function profiling works", {
  setup <- setup_simple_quadratic()
  
  result <- computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "uniform",
      grid_points = 8,
      max_grid_range_multiplier = 0.8
    ),
    optimizer = "optim",
    verbose = FALSE
  )
  
  # Check that profiles completed successfully
  expect_equal(length(result$profiles), 2)
  expect_true(all(sapply(result$profiles, function(p) !p$failed)))
  expect_true(all(result$summary$status == "success"))
  
  # For quadratic function, CI should be symmetric around true values
  # (allowing for some numerical tolerance)
  for (i in 1:2) {
    ci <- result$confidence_intervals[i, ]
    if (!any(is.na(ci))) {
      center <- (ci[1] + ci[2]) / 2
      true_val <- setup$true_params[i]
      expect_true(abs(center - true_val) < 0.1, info = paste("Parameter", i, "CI not centered"))
    }
  }
})

# Test different profile options
test_that("Different profile options work", {
  setup <- setup_gaussian_mixture(d = 3)
  
  # Test uniform grid
  result_uniform <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "uniform",
      grid_points = 6,
      ll_ratio_threshold = 3.84
    ),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Test adaptive grid
  result_adaptive <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 6,
      ll_ratio_threshold = 3.84
    ),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Both should succeed
  expect_true(is.list(result_uniform))
  expect_true(is.list(result_adaptive))
  expect_equal(length(result_uniform$profiles), 3)
  expect_equal(length(result_adaptive$profiles), 3)
  
  # Test invalid grid method
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_method = "invalid"),
      verbose = FALSE,
      gamma = setup$gamma
    ),
    "Unknown grid method"
  )
})

# Test plotting functions
test_that("Plotting functions work", {
  setup <- setup_gaussian_mixture(d = 2)
  
  result <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 5),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Test basic plotting (should not error)
  expect_silent(
    plotLikelihoodProfiles(result, 
                          options = list(max_cols = 2, show_ci = TRUE),
                          true_values = setup$true_params)
  )
  
  # Test plotting with different options
  expect_silent(
    plotLikelihoodProfiles(result,
                          options = list(log_scale = TRUE, show_threshold = FALSE),
                          plot_options = list(cex = 0.8))
  )
  
  # Test invalid input
  expect_error(
    plotLikelihoodProfiles("invalid"),
    "profile_result must be a list from computeLikelihoodProfiles"
  )
})

# Test edge cases and error handling
test_that("Edge cases are handled correctly", {
  setup <- setup_gaussian_mixture(d = 2)
  
  # Test with very small grid
  result_small <- computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 3),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  expect_true(is.list(result_small))
  
  # Test with parameters at boundary
  boundary_params <- c(-1.9, 1.9)
  names(boundary_params) <- c("param_1", "param_2")
  
  result_boundary <- computeLikelihoodProfiles(
    params_current = boundary_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 5),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  expect_true(is.list(result_boundary))
  
  # Check that warnings are issued for boundary issues
  expect_true(any(sapply(result_boundary$profiles, function(p) {
    !is.null(p$exit_flags) && isTRUE(p$exit_flags$optimal_value_out_of_bounds)
  })))
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
      optim_options = list(itermax = 20, NP = 10),  # Small for speed
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

cat("All tests completed successfully!\n")
