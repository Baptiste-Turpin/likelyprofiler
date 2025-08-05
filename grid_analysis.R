# Test script to analyze grid methods
# Testing createUniformGrid
cat("=== UNIFORM GRID ANALYSIS ===\n")

# Test case 1: Symmetric range
grid_lower <- 0
grid_upper <- 10
current_value <- 5
n_points <- 11

# Simulate the uniform grid function
createUniformGrid_test <- function(grid_lower, grid_upper, current_value, n_points) {
  res_grid <- c(seq(grid_lower, current_value, length.out = floor(n_points/2) + 1),
               seq(current_value, grid_upper, length.out = floor(n_points/2) + 1))
  unique(res_grid)
}

uniform_grid <- createUniformGrid_test(grid_lower, grid_upper, current_value, n_points)
cat("Symmetric case (0 to 10, center=5, n=11):\n")
cat("Grid points:", paste(round(uniform_grid, 3), collapse=", "), "\n")
cat("Number of points:", length(uniform_grid), "\n")
cat("Spacing pattern:", paste(round(diff(uniform_grid), 3), collapse=", "), "\n\n")

# Test case 2: Asymmetric case
current_value <- 2
uniform_grid2 <- createUniformGrid_test(grid_lower, grid_upper, current_value, n_points)
cat("Asymmetric case (0 to 10, center=2, n=11):\n")
cat("Grid points:", paste(round(uniform_grid2, 3), collapse=", "), "\n")
cat("Number of points:", length(uniform_grid2), "\n")
cat("Spacing pattern:", paste(round(diff(uniform_grid2), 3), collapse=", "), "\n\n")

cat("=== ADAPTIVE GRID ANALYSIS ===\n")

# Simulate the adaptive grid function
createAdaptiveGrid_test <- function(grid_lower, grid_upper, center_value, n_points) {
  n_half <- floor(n_points/2) + 1
  t_onesided <- seq(0, 1, length.out = n_half)^2
  left_side <- center_value - t_onesided*(center_value - grid_lower)
  right_side <- center_value + t_onesided*(grid_upper - center_value)
  
  res_grid <- c(left_side, right_side)
  unique(sort(res_grid))
}

# Test case 1: Symmetric range
current_value <- 5
adaptive_grid <- createAdaptiveGrid_test(grid_lower, grid_upper, current_value, n_points)
cat("Symmetric case (0 to 10, center=5, n=11):\n")
cat("Grid points:", paste(round(adaptive_grid, 3), collapse=", "), "\n")
cat("Number of points:", length(adaptive_grid), "\n")
cat("Spacing pattern:", paste(round(diff(adaptive_grid), 3), collapse=", "), "\n\n")

# Test case 2: Asymmetric case
current_value <- 2
adaptive_grid2 <- createAdaptiveGrid_test(grid_lower, grid_upper, current_value, n_points)
cat("Asymmetric case (0 to 10, center=2, n=11):\n")
cat("Grid points:", paste(round(adaptive_grid2, 3), collapse=", "), "\n")
cat("Number of points:", length(adaptive_grid2), "\n")
cat("Spacing pattern:", paste(round(diff(adaptive_grid2), 3), collapse=", "), "\n\n")

cat("=== ANALYSIS SUMMARY ===\n")
cat("Uniform Grid Behavior:\n")
cat("- Creates two separate linear sequences\n")
cat("- Left: grid_lower to current_value\n") 
cat("- Right: current_value to grid_upper\n")
cat("- May have uneven spacing when current_value is not centered\n\n")

cat("Adaptive Grid Behavior:\n")
cat("- Uses quadratic transformation t^2 for spacing\n")
cat("- Creates denser points near center (current_value)\n")
cat("- Spacing grows quadratically as distance from center increases\n")
cat("- More sophisticated density control\n")
