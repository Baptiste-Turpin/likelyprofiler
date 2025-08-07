# Test bootstrap validation and edge cases

test_that("Bootstrap validation catches missing generateData", {
  setup = setup_normal_bootstrap()

  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        threshold_method = "bootstrap",
        n_bootstrap = 10
      ),
      dataset = setup$data,
      verbose = FALSE,
      # Missing generateData
    ),
    "generateData.*required.*bootstrap"
  )
})

test_that("Bootstrap validation catches invalid n_bootstrap", {
  setup = setup_normal_bootstrap()

  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        threshold_method = "bootstrap",
        n_bootstrap = 0  # Invalid
      ),
      generateData = setup$generateData,
      dataset = setup$data
    ),
    "n_bootstrap.*positive"
  )

  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        threshold_method = "bootstrap",
        n_bootstrap = -5  # Invalid
      ),
      generateData = setup$generateData,
      dataset = setup$data
    ),
    "n_bootstrap.*positive"
  )
})

test_that("Bootstrap validation catches invalid confidence level", {
  setup = setup_normal_bootstrap()

  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        threshold_method = "bootstrap",
        n_bootstrap = 10,
        bootstrap_conf_level = 0  # Invalid
      ),
      generateData = setup$generateData,
      dataset = setup$data
    ),
    "bootstrap_conf_level.*between 0 and 1"
  )

  expect_error(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        threshold_method = "bootstrap",
        n_bootstrap = 10,
        bootstrap_conf_level = 1.5  # Invalid
      ),
      generateData = setup$generateData,
      dataset = setup$data
    ),
    "bootstrap_conf_level.*between 0 and 1"
  )
})

test_that("generateData function validation works", {
  setup = setup_normal_bootstrap()

  # Test that generateData is called with correct parameters
  mock_generateData = function(params) {
    expect_equal(length(params), 2)
    expect_true(all(is.finite(params)))
    expect_true(params[2] > 0)  # sigma > 0
    return(rnorm(setup$n, params[1], params[2]))
  }

  expect_silent(
    computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 10,
        threshold_method = "bootstrap",
        n_bootstrap = 5
      ),
      generateData = mock_generateData,
      dataset = setup$data,
      verbose = FALSE
    )
  )
})

test_that("Edge case: optimization failure handling", {
  setup = setup_normal_bootstrap()

  # Test with bad generateData that could cause optimization issues
  bad_generateData = function(params) {
    # Generate extreme data that might cause optimization issues
    rep(1e10, setup$n)
  }

  # Should still run without crashing (though may have poor results)
  expect_silent({
    result = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 10,
        threshold_method = "bootstrap",
        n_bootstrap = 5
      ),
      generateData = bad_generateData,
      dataset = setup$data,
      verbose = FALSE
    )
  })

  # Should still return proper structure even with bad bootstrap data
  expect_true(is.list(result))
  expect_true("bootstrap_stats" %in% names(result))
})

test_that("Bootstrap with different confidence levels", {
  setup = setup_normal_bootstrap()

  expect_warning({
    # Test 90% confidence level
    result_90 = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 3,
        threshold_method = "bootstrap",
        bootstrap_conf_level = 0.90,
        n_bootstrap = 10
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  expect_warning({
    # Test 99% confidence level
    result_99 = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 3,
        threshold_method = "bootstrap",
        bootstrap_conf_level = 0.99,
        n_bootstrap = 10
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  # Both should work and have proper structure
  expect_true(is.list(result_90))
  expect_true(is.list(result_99))
  expect_equal(dim(result_90$bootstrap_stats), dim(result_99$bootstrap_stats))

  # 99% intervals should generally be wider than 90% intervals
  width_90 = result_90$confidence_intervals[,2] - result_90$confidence_intervals[,1]
  width_99 = result_99$confidence_intervals[,2] - result_99$confidence_intervals[,1]
  expect_true(all(width_99 >= width_90))
})

test_that("Edge case: very small bootstrap sample", {
  setup = setup_normal_bootstrap()

  expect_warning({
    # Test with minimal bootstrap sample (should still work)
    result = computeLikelihoodProfiles(
      params_current = setup$ml_params,
      negLogLikelihood = setup$negLogLik,
      bounds = setup$bounds,
      profile_options = list(
        grid_points = 3,
        threshold_method = "bootstrap",
        n_bootstrap = 3  # Very small
      ),
      generateData = setup$generateData,
      dataset = setup$data,
      verbose = FALSE
    )},
    "Warning in computeLikelihoodProfiles: Few grid points \\(<3\\) within threshold for parameter\\(s\\): .* Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
    perl = TRUE
  )

  expect_true(is.list(result))
  expect_true("bootstrap_stats" %in% names(result))
  expect_equal(dim(result$bootstrap_stats)[1], 3)
  expect_false(any(is.na(result$confidence_intervals)))
})
