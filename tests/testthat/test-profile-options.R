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
    "profile_options\\$grid_method must be one of: uniform, adaptive"
  )
})
