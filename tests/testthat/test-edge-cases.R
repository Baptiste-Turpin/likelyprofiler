# Test edge cases and error handling

test_that("Edge cases are handled correctly", {
  setup = setup_gaussian_mixture(d = 2)

  expect_warning({
    # Test with very small grid
    result_small = computeLikelihoodProfiles(
      params_current = setup$optimized_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 3),
      verbose = FALSE,
      gamma = setup$gamma
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
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

  expect_warning({
    result_boundary_exact = computeLikelihoodProfiles(
      params_current = boundary_exact_params,
      negLogLikelihood = setup$mlogf,
      bounds = setup$bounds,
      profile_options = list(grid_points = 5, ll_ratio_threshold = 50),
      verbose = FALSE,
      gamma = setup$gamma
    )},
    "Grid too narrow to capture confidence interval bounds for parameter\\(s\\): .* Consider increasing max_grid_range_multiplier",
    perl = TRUE
  )

  expect_true(is.list(result_boundary_exact))
  expect_equal(length(result_boundary_exact$profiles), 2)
})

test_that("Grid too narrow warning is properly issued", {
  setup = setup_gaussian_mixture(d = 2, gamma = 0.1)  # Very flat likelihood

  # Capture warnings during profiling with small grid and flat likelihood
  expect_warning(
    {
      result_narrow = computeLikelihoodProfiles(
        params_current = setup$optimized_params,
        negLogLikelihood = setup$mlogf,
        bounds = setup$bounds,
        profile_options = list(
          grid_points = 5,  # Small grid
          max_grid_range_multiplier = 0.1,  # Very narrow grid range
          ll_ratio_threshold = 3.84
        ),
        verbose = FALSE,
        gamma = setup$gamma
      )},
    "Grid too narrow to capture confidence interval bounds for parameter\\(s\\): .* Consider increasing max_grid_range_multiplier",
    perl = TRUE
  )

  # Verify result structure is still valid
  expect_true(is.list(result_narrow))
  expect_true(all(c("profiles", "confidence_intervals", "summary") %in% names(result_narrow)))

  # Check that exit flags are properly set
  expect_true(all(sapply(result_narrow$profiles, function(p) isTRUE(p$exit_flags$grid_too_narrow))))

  # Confidence intervals should still be returned (as grid bounds)
  expect_true(all(!is.na(result_narrow$confidence_intervals)))
})
