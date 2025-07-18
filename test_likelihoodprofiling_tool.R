

mlogf = function(x, gamma){
  m1 = m2 = rep(0, length(x))
  m1[1] = 1
  m2[1] = -1
  - log(exp(- gamma * (sum((x-m1)^2))) + exp(- gamma * (sum((x-m2)^2))))
}
#Problem parameters
gamma = 8   #power law for the Gaussian mixture
d = 5

# Set up the problem
true_params = c(1, rep(0, d-1))
names(true_params) = paste0("param_", 1:d)
bounds = list(lower = rep(-2, d), upper = rep(2, d))

# Compute profile likelihood with small grid for speed
result = computeLikelihoodProfiles(
  params_current = true_params,
  negLogLikelihood = mlogf,
  bounds = bounds,
  profile_options = list(
    grid_method = "adaptive",
    grid_points = 50,
    max_grid_range_multiplier = 1.0
  ),
  optimizer = "optim",
  verbose = FALSE,  # Suppress output for clean example
  method = "L-BFGS-B",
  control = list(maxit = 100),
  gamma = gamma
)

# View results
print(result$summary)
print(result$confidence_intervals)

# Plot the profile likelihood
plotLikelihoodProfiles(result,
                       options = list(log_scale = FALSE, max_cols = 2),
                       true_values = true_params
)
