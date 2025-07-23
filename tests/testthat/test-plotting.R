# Test plotting functions

test_that("Plotting functions work", {
  setup = setup_gaussian_mixture(d = 2)

  result = computeLikelihoodProfiles(
    params_current = setup$optimized_params,
    negLogLikelihood = setup$mlogf,
    bounds = setup$bounds,
    profile_options = list(grid_points = 25),
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
    "Error in plotLikelihoodProfiles: profile_result must be a list from computeLikelihoodProfiles"
  )
})
