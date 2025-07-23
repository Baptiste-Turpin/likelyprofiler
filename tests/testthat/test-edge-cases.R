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
  boundary_params <- c(-2.1, 1.9)
  names(boundary_params) <- c("param_1", "param_2")

  # Setting high ll_ratio threshold as params_current are at the boundary and
  # therefore not optimized (will trigger warning)
  result_boundary <- computeLikelihoodProfiles(
    params_current = boundary_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 5, ll_ratio_threshold = 80),
    verbose = FALSE,
    gamma = setup$gamma
  )

  expect_true(is.list(result_boundary))

  # Check that warnings are issued for boundary issues for first parameter
  vec_logical_outofbounds = sapply(result_boundary$profiles, function(p) {
    !is.null(p$exit_flags) && isTRUE(p$exit_flags$optimal_value_out_of_bounds)
  })
  expect_true(vec_logical_outofbounds[1])  # First profile should have out of bounds warning
  expect_false(any(vec_logical_outofbounds[-1])) # But not the others
})
