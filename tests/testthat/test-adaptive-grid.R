# Test adaptive grid functionality

test_that("Adaptive grid basic functionality works", {
  setup = setup_simple_quadratic()
  
  result = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 20,  # Small for speed
      max_grid_range_multiplier = 1.0,
      adaptive_max_iterations = 2,
      adaptive_initial_fraction = 0.4,
      adaptive_refinement_fraction = 0.6
    ),
    optimizer = "optim",
    optim_options = list(
      method = "L-BFGS-B",
      control = list(maxit = 50)
    ),
    verbose = FALSE
  )
  
  # Check basic result structure
  expect_true(is.list(result))
  expect_true(all(c("profiles", "confidence_intervals", "summary", "options", "optimizer") %in% names(result)))
  expect_equal(length(result$profiles), 2)  # Two parameters
  expect_equal(nrow(result$confidence_intervals), 2)
  expect_equal(result$options$grid_method, "adaptive")
  
  # Check that profiles were computed successfully
  expect_true(all(!is.na(result$confidence_intervals[, 1])))
  expect_true(all(!is.na(result$confidence_intervals[, 2])))
  
  # Check that adaptive-specific options were used
  expect_equal(result$options$adaptive_max_iterations, 2)
  expect_equal(result$options$adaptive_initial_fraction, 0.4)
  expect_equal(result$options$adaptive_refinement_fraction, 0.6)
})

test_that("Adaptive grid vs quadratic grid comparison", {
  setup = setup_simple_quadratic()
  
  # Run with quadratic grid
  result_quad = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "quadratic",
      grid_points = 25,
      max_grid_range_multiplier = 1.0
    ),
    optimizer = "optim",
    optim_options = list(method = "L-BFGS-B", control = list(maxit = 50)),
    verbose = FALSE
  )
  
  # Run with adaptive grid
  result_adaptive = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 25,
      max_grid_range_multiplier = 1.0,
      adaptive_max_iterations = 3
    ),
    optimizer = "optim",
    optim_options = list(method = "L-BFGS-B", control = list(maxit = 50)),
    verbose = FALSE
  )
  
  # Both should succeed
  expect_true(all(!is.na(result_quad$confidence_intervals)))
  expect_true(all(!is.na(result_adaptive$confidence_intervals)))
  
  # Check that confidence intervals are reasonably close (should be similar for simple quadratic)
  ci_diff = abs(result_quad$confidence_intervals - result_adaptive$confidence_intervals)
  expect_true(all(ci_diff < 1.0))  # Should be within 1 unit for this simple problem
})

test_that("Adaptive grid parameters validation", {
  setup = setup_simple_quadratic()
  
  # Test invalid adaptive_max_iterations
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$true_params,
      negLogLikelihood = setup$likelihood,
      bounds = setup$bounds,
      profile_options = list(
        grid_method = "adaptive",
        grid_points = 20,
        adaptive_max_iterations = 0  # Invalid
      ),
      verbose = FALSE
    ),
    "adaptive_max_iterations"
  )
  
  # Test invalid adaptive_initial_fraction
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$true_params,
      negLogLikelihood = setup$likelihood,
      bounds = setup$bounds,
      profile_options = list(
        grid_method = "adaptive",
        grid_points = 20,
        adaptive_initial_fraction = 1.5  # Invalid (> 1)
      ),
      verbose = FALSE
    ),
    "adaptive_initial_fraction"
  )
  
  # Test invalid adaptive_refinement_fraction
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$true_params,
      negLogLikelihood = setup$likelihood,
      bounds = setup$bounds,
      profile_options = list(
        grid_method = "adaptive",
        grid_points = 20,
        adaptive_refinement_fraction = 0  # Invalid
      ),
      verbose = FALSE
    ),
    "adaptive_refinement_fraction"
  )
})

test_that("Adaptive grid works with different grid sizes", {
  setup = setup_simple_quadratic()
  
  # Test small grid
  result_small = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 10,
      adaptive_max_iterations = 2
    ),
    optimizer = "optim",
    optim_options = list(method = "L-BFGS-B", control = list(maxit = 50)),
    verbose = FALSE
  )
  
  # Test larger grid
  result_large = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 50,
      adaptive_max_iterations = 3
    ),
    optimizer = "optim",
    optim_options = list(method = "L-BFGS-B", control = list(maxit = 50)),
    verbose = FALSE
  )
  
  # Both should succeed
  expect_true(all(!is.na(result_small$confidence_intervals)))
  expect_true(all(!is.na(result_large$confidence_intervals)))
  
  # Larger grid should generally provide more precise estimates (smaller CI width)
  ci_width_small = result_small$summary$ci_width
  ci_width_large = result_large$summary$ci_width
  
  # At least one parameter should have smaller CI with larger grid
  expect_true(any(ci_width_large <= ci_width_small * 1.1))  # Allow some tolerance
})

test_that("Adaptive grid helper functions work correctly", {
  # Test createAdaptiveInitialGrid
  grid = likelyprofiler:::createAdaptiveInitialGrid(0, 10, 5, 15)
  expect_true(length(grid) <= 15)
  expect_true(all(grid >= 0 & grid <= 10))
  expect_true(5 %in% grid)  # Center value should be included
  
  # Test that grid is sorted
  expect_equal(grid, sort(grid))
})

test_that("Adaptive grid works with single iteration", {
  setup = setup_simple_quadratic()
  
  result = computeLikelihoodProfiles(
    params_current = setup$true_params,
    negLogLikelihood = setup$likelihood,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 20,
      adaptive_max_iterations = 1,  # Single iteration only
      adaptive_initial_fraction = 1.0  # Use all points in first iteration
    ),
    optimizer = "optim",
    optim_options = list(method = "L-BFGS-B", control = list(maxit = 50)),
    verbose = FALSE
  )
  
  # Should still work with single iteration
  expect_true(all(!is.na(result$confidence_intervals)))
  expect_equal(result$options$adaptive_max_iterations, 1)
})

test_that("Adaptive grid works with Gaussian mixture", {
  setup = setup_gaussian_mixture(d = 3)
  
  result = computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "adaptive",
      grid_points = 15,  # Small for speed
      max_grid_range_multiplier = 1.0,
      adaptive_max_iterations = 2
    ),
    optimizer = "optim",
    optim_options = list(
      method = "L-BFGS-B",
      control = list(maxit = 50)
    ),
    verbose = FALSE,
    gamma = setup$gamma
  )
  
  # Check that it works with more complex likelihood
  expect_true(is.list(result))
  expect_equal(length(result$profiles), 3)
  expect_equal(nrow(result$confidence_intervals), 3)
  
  # At least some profiles should succeed
  expect_true(any(result$summary$status == "success"))
})
