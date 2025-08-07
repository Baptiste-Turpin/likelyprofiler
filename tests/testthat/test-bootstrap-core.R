# Test core bootstrap functionality

test_that("getBootstrapLRTThreshold returns correct structure", {
  # Setup simple test case
  setup = setup_normal_bootstrap()

  # Test non-bootstrap fallback
  result_fixed = getBootstrapLRTThreshold(
    optimal_params = setup$ml_params,
    param_index = 1,
    is_success = TRUE,
    fixed_value = setup$ml_params[1],
    profile_options = list(threshold_method = "fixed", ll_ratio_threshold = 3.84),
    negLogLikelihood = setup$negLogLik,
    bounds = setup$bounds,
    optimizer_info = list(name = "optim", type = "builtin", extra_args = list()),
    likelihood_args = list(dataset = setup$data),
    verbose = FALSE
  )

  expect_true(is.list(result_fixed))
  expect_true(all(c("threshold", "lrt_stats") %in% names(result_fixed)))
  expect_equal(result_fixed$threshold, 3.84)
  expect_null(result_fixed$lrt_stats)

  # Test bootstrap method with very small sample for speed
  result_bootstrap = getBootstrapLRTThreshold(
    optimal_params = setup$ml_params,
    param_index = 1,
    is_success = TRUE,
    fixed_value = setup$ml_params[1],
    profile_options = list(
      threshold_method = "bootstrap",
      n_bootstrap = 10,  # Very small for speed
      bootstrap_conf_level = 0.95
    ),
    negLogLikelihood = setup$negLogLik,
    bounds = setup$bounds,
    optimizer_info = list(name = "optim", type = "builtin", extra_args = list(method = "L-BFGS-B")),
    likelihood_args = list(dataset = setup$data),
    generateData = setup$generateData,
    verbose = FALSE
  )

  expect_true(is.list(result_bootstrap))
  expect_true(all(c("threshold", "lrt_stats") %in% names(result_bootstrap)))
  expect_true(is.numeric(result_bootstrap$threshold))
  expect_true(is.numeric(result_bootstrap$lrt_stats))
  expect_equal(length(result_bootstrap$lrt_stats), 10)
})

test_that("Bootstrap integration with computeLikelihoodProfiles works", {
  setup = setup_normal_bootstrap()

  # Test with minimal settings for speed
  expect_warning({
    result_bootstrap = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 5,  # Very small for speed
        threshold_method = "bootstrap",
        bootstrap_conf_level = 0.95,
        n_bootstrap = 10  # Very small for speed
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  # Test structure
  expect_true(is.list(result_bootstrap))
  expect_true("bootstrap_stats" %in% names(result_bootstrap))
  expect_true(all(c("profiles", "confidence_intervals", "summary", "options", "optimizer") %in% names(result_bootstrap)))

  # Test bootstrap_stats array dimensions [n_bootstrap, n_params, n_grid]
  expect_equal(dim(result_bootstrap$bootstrap_stats), c(10, 2, 5))

  # Test that some bootstrap stats are computed (not all NA)
  expect_true(sum(!is.na(result_bootstrap$bootstrap_stats)) > 0)

  # Test that confidence intervals are computed
  expect_false(any(is.na(result_bootstrap$confidence_intervals)))
  expect_equal(nrow(result_bootstrap$confidence_intervals), 2)
})

test_that("Bootstrap statistics collection works", {
  setup = setup_normal_bootstrap()

  expect_warning({
    result = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 3,  # Minimal for speed
        threshold_method = "bootstrap",
        n_bootstrap = 5   # Minimal for speed
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  # Test array structure
  expect_true(is.array(result$bootstrap_stats))
  expect_equal(length(dim(result$bootstrap_stats)), 3)
  expect_equal(dim(result$bootstrap_stats)[1], 5)  # n_bootstrap
  expect_equal(dim(result$bootstrap_stats)[2], 2)  # n_params
  expect_equal(dim(result$bootstrap_stats)[3], 3)  # n_grid

  # Test that bootstrap stats are numeric
  expect_true(is.numeric(result$bootstrap_stats))
})

test_that("generateData function integration works", {
  setup = setup_normal_bootstrap()

  # Test that generateData is called and works
  test_data = setup$generateData(setup$ml_params)
  expect_true(is.numeric(test_data))
  expect_equal(length(test_data), setup$n)

  # Test that generated data can be used with negLogLik
  test_cost = setup$negLogLik(setup$ml_params, dataset = test_data)
  expect_true(is.numeric(test_cost))
  expect_true(is.finite(test_cost))
})

test_that("Bootstrap works with fixed threshold method for comparison", {
  setup = setup_normal_bootstrap()

  expect_warning({
    # Fixed method
    result_fixed = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 5,
        threshold_method = "fixed",
        ll_ratio_threshold = qchisq(0.95, 1)
      ),
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  expect_warning({
    # Bootstrap method
    result_bootstrap = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 5,
        threshold_method = "bootstrap",
        bootstrap_conf_level = 0.95,
        n_bootstrap = 10
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  # Both should have similar structure
  expect_equal(names(result_fixed), names(result_bootstrap))
  expect_equal(dim(result_fixed$confidence_intervals), dim(result_bootstrap$confidence_intervals))

  # Only bootstrap should have bootstrap_stats
  expect_null(result_fixed$bootstrap_stats)
  expect_true(!is.null(result_bootstrap$bootstrap_stats))
})
