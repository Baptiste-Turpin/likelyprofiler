# Test different profile options

test_that("Different profile options work", {
  setup = setup_gaussian_mixture(d = 3)

  # Test linear grid
  expect_warning({
    result_linear = computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(
        grid_method = "linear",
        grid_points = 6,
        ll_ratio_threshold = 3.84
      ),
      verbose = FALSE,
      gamma = setup$gamma
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  # Test quadratic grid
  result_quadratic = computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(
      grid_method = "quadratic",
      grid_points = 6,
      ll_ratio_threshold = 3.84
    ),
    verbose = FALSE,
    gamma = setup$gamma
  )

  # Both should succeed
  expect_true(is.list(result_linear))
  expect_true(is.list(result_quadratic))
  expect_equal(length(result_linear$profiles), 3)
  expect_equal(length(result_quadratic$profiles), 3)

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
    "Error in validateProfileOptions: profile_options\\$grid_method must be one of: linear, quadratic"
  )
})
