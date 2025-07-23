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

test_that("Profile options validation works correctly", {
  setup <- setup_gaussian_mixture(d = 2)

  # Test invalid grid_method
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_method = "invalid_method"),
      gamma = setup$gamma
    ),
    "profile_options\\$grid_method must be one of: uniform, adaptive"
  )

  # Test invalid grid_points - non-numeric
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = "not_numeric"),
      gamma = setup$gamma
    ),
    "profile_options\\$grid_points must be a single numeric value"
  )

  # Test invalid grid_points - too small
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 2),
      gamma = setup$gamma
    ),
    "profile_options\\$grid_points must be at least 3"
  )

  # Test invalid grid_points - not integer
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 5.5),
      gamma = setup$gamma
    ),
    "profile_options\\$grid_points must be an integer"
  )

  # Test invalid max_grid_range_multiplier - non-numeric
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(max_grid_range_multiplier = "not_numeric"),
      gamma = setup$gamma
    ),
    "profile_options\\$max_grid_range_multiplier must be a single numeric value"
  )

  # Test invalid max_grid_range_multiplier - negative
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(max_grid_range_multiplier = -1),
      gamma = setup$gamma
    ),
    "profile_options\\$max_grid_range_multiplier must be positive"
  )

  # Test invalid ll_ratio_threshold - non-numeric
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(ll_ratio_threshold = "not_numeric"),
      gamma = setup$gamma
    ),
    "profile_options\\$ll_ratio_threshold must be a single numeric value"
  )

  # Test invalid ll_ratio_threshold - negative
  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(ll_ratio_threshold = -1),
      gamma = setup$gamma
    ),
    "profile_options\\$ll_ratio_threshold must be positive"
  )

  # Test warning for large grid_points
  expect_warning(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 1500),
      verbose = FALSE,
      gamma = setup$gamma
    ),
    "profile_options\\$grid_points is very large"
  )

  # Test warning for large max_grid_range_multiplier
  expect_warning(
    computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(max_grid_range_multiplier = 15),
      verbose = FALSE,
      gamma = setup$gamma
    ),
    "It is not useful to set profile_options\\$max_grid_range_multiplier greater than 2.0, as this already covers the whole parameter range."
  )
})
