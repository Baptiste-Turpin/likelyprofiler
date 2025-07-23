# Test helper functions and setup utilities

# Test data setup for Gaussian mixture
setup_gaussian_mixture <- function(d = 5, gamma = 8) {
  # Define the Gaussian mixture negative log-likelihood
  mlogf <- function(x, gamma) {
    m1 <- m2 <- rep(0, length(x))
    m1[1] <- 1
    m2[1] <- 0
    -log(exp(-gamma * sum((x - m1)^2)) + exp(-gamma * sum((x - m2)^2)))
  }

  # Problem parameters
  true_params <- c(1, rep(0, d - 1))
  set.seed(23)
  optimized_params <- true_params + rnorm(length(true_params), sd = 0.05)
  names(true_params) <- names(optimized_params) <- paste0("param_", 1:d)
  bounds <- list(lower = rep(-2, d), upper = rep(2, d))

  list(
    mlogf = mlogf,
    gamma = gamma,
    true_params = true_params,
    optimized_params = optimized_params,
    bounds = bounds,
    d = d
  )
}

# Simple quadratic function setup for testing
setup_simple_quadratic <- function() {
  # Simple quadratic function for testing
  quadratic_likelihood <- function(params) {
    x <- params[1]
    y <- params[2]
    (x - 2)^2 + 2 * (y - 1)^2
  }

  true_params <- c(2.0, 1.0)
  names(true_params) <- c("param_x", "param_y")
  bounds <- list(lower = c(0, 0), upper = c(4, 3))

  list(
    likelihood = quadratic_likelihood,
    true_params = true_params,
    bounds = bounds
  )
}
