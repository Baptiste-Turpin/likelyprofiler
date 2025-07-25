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

# Normal distribution setup for bootstrap testing
setup_normal_bootstrap <- function() {
  # Generate synthetic normal data
  n <- 50  # Small for fast tests
  true_mu <- 2.5
  true_sigma <- 1.2
  set.seed(42)  # Fixed seed for reproducible tests
  data_normal <- rnorm(n, true_mu, true_sigma)
  
  # Set up bounds
  bounds <- list(lower = c(-1, 0.1), upper = c(6, 4))
  
  # Define negative log-likelihood for normal distribution
  negLogLikNormal <- function(params, dataset = data_normal) {
    mu <- params[1]
    sigma <- params[2]
    if (sigma <= 0) return(1e10)  # Constraint: sigma > 0
    -sum(dnorm(dataset, mu, sigma, log = TRUE))
  }
  
  # Find ML estimates
  ml_result <- optim(c(mean(data_normal), sd(data_normal)), negLogLikNormal, 
                     dataset = data_normal, method = "L-BFGS-B", 
                     lower = bounds$lower, upper = bounds$upper)
  ml_params <- ml_result$par
  names(ml_params) <- c("mu", "sigma")
  
  # Create bootstrap data generation function
  generateData <- function(params) {
    rnorm(n, params[1], params[2])
  }
  
  list(
    data = data_normal,
    bounds = bounds,
    negLogLik = negLogLikNormal,
    ml_params = ml_params,
    generateData = generateData,
    n = n,
    true_mu = true_mu,
    true_sigma = true_sigma
  )
}
