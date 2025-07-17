

m_scale = 1
irregular_likelihood <- function(params, m_scale) {
  x <- params[1]
  y <- params[2]

  thr = 0.5
  multiplicator = ((m_scale-m_scale*y)^3 +  m_scale*y^3) * (x<thr) + 1 * (x>=thr)
  res = abs(sin(x*2*pi)) * multiplicator
  - res
}

# Set up the problem
true_params <- c(0.25, 1.0)
names(true_params) <- c("param_x", "param_y")
bounds <- list(lower = c(0, 0), upper = c(1, 1))

# Compute profile likelihood with small grid for speed
result <- compute_all_parameter_profiles(
  params_current = true_params,
  cost_to_optimize = irregular_likelihood,
  bounds = bounds,
  profile_options = list(
    grid_method = "adaptive",
    grid_points = 15,  # Small grid for quick execution
    max_grid_range_multiplier = 1.0
  ),
  optimizer = "optim",
  verbose = FALSE,  # Suppress output for clean example
  method = "L-BFGS-B",
  control = list(maxit = 100),
  m_scale = m_scale
)

# View results
print(result$summary)
print(result$confidence_intervals)

# Plot the profile likelihood
plot(result$profiles[[1]]$grid_values, result$profiles[[1]]$profile_costs, type = "l")
plot(result$profiles[[2]]$grid_values, result$profiles[[2]]$profile_costs, type = "l")
