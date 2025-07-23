# Test edge cases and error handling

test_that("Edge cases are handled correctly", {
  setup = setup_gaussian_mixture(d = 2)

  # Test with very small grid
  result_small = computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 3),
    verbose = FALSE,
    gamma = setup$gamma
  )

  expect_true(is.list(result_small))

  # Test with parameters beyond boundary
  boundary_params = c(-2.1, 1.9)
  expect_error(
    computeLikelihoodProfiles(
      params_current = boundary_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 5),
      verbose = FALSE,
      gamma = setup$gamma
    ),
    "Error in computeLikelihoodProfiles: params_current must be within the specified bounds"
  )

  # Test with parameters exactly on the boundary (should be allowed)
  boundary_exact_params = c(-2.0, 2.0)  # Exactly on bounds

  result_boundary_exact = computeLikelihoodProfiles(
    params_current = boundary_exact_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 5, ll_ratio_threshold = 50),
    verbose = FALSE,
    gamma = setup$gamma
  )

  expect_true(is.list(result_boundary_exact))
  expect_equal(length(result_boundary_exact$profiles), 2)
})
