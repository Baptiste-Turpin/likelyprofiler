# Test core profiling functionality

test_that("Multidimensional Gaussian mixture profiling works", {
  setup = setup_gaussian_mixture(d = 5)

  result = computeLikelihoodProfiles(
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

test_that("One-dimensional Gaussian mixture profiling works", {
  setup = setup_gaussian_mixture(d = 1)

  result = computeLikelihoodProfiles(
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
  profile = result$profiles[[1]]
  expect_false(profile$failed)
  expect_true(length(profile$grid_values) >= 3)
  expect_true(length(profile$profile_costs) >= 3)
})

test_that("Simple quadratic function profiling works", {
  setup = setup_simple_quadratic()

  result = computeLikelihoodProfiles(
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
    ci = result$confidence_intervals[i, ]
    if (!any(is.na(ci))) {
      center = (ci[1] + ci[2]) / 2
      true_val = setup$true_params[i]
      expect_true(abs(center - true_val) < 0.1, info = paste("Parameter", i, "CI not centered"))
    }
  }
})

# test_that("Parallel computation works with cluster", {
#   # Skip if parallel package is not available
#   skip_if_not_installed("parallel")
#
#   setup = setup_gaussian_mixture(d = 3)
#
#   # Create a simple cluster with 2 cores
#   cluster = parallel::makeCluster(2)
#   on.exit(parallel::stopCluster(cluster))
#
#   # Test parallel execution
#   result_parallel = computeLikelihoodProfiles(
#     params_current = setup$optimized_params,
#     negLogLikelihood = setup$mlogf,
#     bounds = setup$bounds,
#     profile_options = list(grid_points = 5),  # Small for speed
#     cluster = cluster,
#     verbose = FALSE,
#     gamma = setup$gamma
#   )
#
#   # Check that results are valid
#   expect_true(is.list(result_parallel))
#   expect_equal(length(result_parallel$profiles), 3)
#   expect_true(all(c("profiles", "confidence_intervals", "summary") %in% names(result_parallel)))
#
#   # Compare with sequential computation
#   result_sequential = computeLikelihoodProfiles(
#     params_current = setup$optimized_params,
#     negLogLikelihood = setup$mlogf,
#     bounds = setup$bounds,
#     profile_options = list(grid_points = 5),  # Same parameters
#     cluster = NULL,
#     verbose = FALSE,
#     gamma = setup$gamma
#   )
#
#   # Results should have same structure (though values may differ slightly due to numerical precision)
#   expect_equal(length(result_parallel$profiles), length(result_sequential$profiles))
#   expect_equal(dim(result_parallel$confidence_intervals), dim(result_sequential$confidence_intervals))
#   expect_equal(nrow(result_parallel$summary), nrow(result_sequential$summary))
# })
