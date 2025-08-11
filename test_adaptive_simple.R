# Simple standalone test for adaptive grid functionality
# This script can be run independently to test basic adaptive grid features


library(likelyprofiler)
# Simple quadratic function
quadratic_likelihood <- function(params) {
  x <- params[1]
  y <- params[2]
  (x - 2)^2 + 2 * (y - 1)^2
}

# Setup parameters
true_params <- c(2.0, 1.0)
names(true_params) <- c("param_x", "param_y")
bounds <- list(lower = c(0, 0), upper = c(4, 3))

result_compare = computeLikelihoodProfiles(
  params_current = true_params,
  negLogLikelihood = quadratic_likelihood,
  bounds = bounds,
  profile_options = list(
    grid_method = "quadratic",
    grid_points = 30,
    max_grid_range_multiplier = 1.0,
    adaptive_max_iterations = 2,
    adaptive_initial_fraction = 0.4,
    adaptive_refinement_fraction = 0.6
  ),
  optimizer = "optim",
  optim_options = list(
    method = "L-BFGS-B",
    control = list(maxit = 500)
  ),
  verbose = TRUE
)


# Test adaptive grid
result = computeLikelihoodProfiles(
  params_current = true_params,
  negLogLikelihood = quadratic_likelihood,
  bounds = bounds,
  profile_options = list(
    grid_method = "adaptive",
    grid_points = 30,
    max_grid_range_multiplier = 1.0,
    adaptive_max_iterations = 2,
    adaptive_initial_fraction = 0.4,
    adaptive_refinement_fraction = 0.6
  ),
  optimizer = "optim",
  optim_options = list(
    method = "L-BFGS-B",
    control = list(maxit = 500)
  ),
  verbose = TRUE
)

plotLikelihoodProfiles(result)
plotLikelihoodProfiles(result_compare)


cat("Confidence intervals:\n")
print(result$confidence_intervals)
cat("Summary:\n")
print(result$summary)
