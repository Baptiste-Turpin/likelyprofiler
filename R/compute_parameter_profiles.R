#' Compute Profile Likelihood for All Parameters
#'
#' Compute profile likelihood confidence intervals for all parameters in parallel
#' using a grid-based approach with customizable optimization methods.
#'
#' @param params_current Numeric vector of population parameters at optimum
#' @param cost_to_optimize Function that takes parameter vector and returns scalar cost
#' @param bounds List with elements 'lower' and 'upper' containing parameter bounds
#' @param profile_options List of profiling options (see Details)
#' @param optimizer Character string or function specifying optimizer (default: "optim")
#' @param verbose Logical indicating whether to print progress messages
#' @param n_cores Integer number of cores for parallel computation
#' @param ... Additional arguments passed to the optimizer
#'
#' @details
#' The \code{profile_options} list can contain:
#' \itemize{
#'   \item \code{grid_method}: "uniform" or "adaptive" (default: "uniform")
#'   \item \code{grid_points}: Number of grid points per parameter (default: 25)
#'   \item \code{max_grid_range_multiplier}: Grid range multiplier (default: 2.0)
#'   \item \code{ll_ratio_threshold}: Likelihood ratio threshold for CI (default: 3.84)
#' }
#'
#' The \code{optimizer} can be:
#' \itemize{
#'   \item \code{"optim"}: Use base R optim() with method and control via ...
#'   \item \code{"deoptim"}: Use DEoptim package (must be installed)
#'   \item A function with signature: function(fn, par, lower, upper, ...)
#' }
#'
#' @return List containing:
#' \itemize{
#'   \item \code{profiles}: List of profile data for each parameter
#'   \item \code{confidence_intervals}: Matrix of [lower, upper] bounds for each parameter
#'   \item \code{summary}: Data frame summarizing results
#' }
#'
#' @examples
#' # Minimal working example with irregular likelihood surface
#'
#' # Define a multimodal likelihood function with local minima
#' irregular_likelihood <- function(params) {
#'   x <- params[1]
#'   y <- params[2]
#'
#'   # Main quadratic bowl with global minimum at (2, 1)
#'   main_cost <- (x - 2)^2 + 2 * (y - 1)^2
#'
#'   # Add local minima and ridges to create irregularities
#'   ridge1 <- 0.5 * exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.1)
#'   ridge2 <- 0.3 * exp(-((x - 3.5)^2 + (y - 1.8)^2) / 0.15)
#'
#'   # Add some noise to make optimization challenging
#'   noise <- 0.1 * sin(10 * x) * cos(8 * y)
#'
#'   return(main_cost - ridge1 - ridge2 + noise)
#' }
#'
#' # Set up the problem
#' true_params <- c(2.0, 1.0)
#' names(true_params) <- c("param_x", "param_y")
#' bounds <- list(lower = c(0, 0), upper = c(4, 3))
#'
#' # Compute profile likelihood with small grid for speed
#' result <- compute_all_parameter_profiles(
#'   params_current = true_params,
#'   cost_to_optimize = irregular_likelihood,
#'   bounds = bounds,
#'   profile_options = list(
#'     grid_method = "adaptive",
#'     grid_points = 15,  # Small grid for quick execution
#'     max_grid_range_multiplier = 1.0
#'   ),
#'   optimizer = "optim",
#'   verbose = FALSE,  # Suppress output for clean example
#'   method = "L-BFGS-B",
#'   control = list(maxit = 100)
#' )
#'
#' # View results
#' print(result$summary)
#' print(result$confidence_intervals)
#'
#' \dontrun{
#' # More comprehensive example with different optimizers
#' result_deoptim <- compute_all_parameter_profiles(
#'   params_current = true_params,
#'   cost_to_optimize = irregular_likelihood,
#'   bounds = bounds,
#'   optimizer = "deoptim",
#'   itermax = 100,
#'   NP = 20
#' )
#'
#' # Custom optimizer example
#' my_optimizer <- function(fn, par, lower, upper, ...) {
#'   optim(par = par, fn = fn, method = "Nelder-Mead", ...)
#' }
#'
#' result_custom <- compute_all_parameter_profiles(
#'   params_current = true_params,
#'   cost_to_optimize = irregular_likelihood,
#'   bounds = bounds,
#'   optimizer = my_optimizer
#' )
#' }
#'
#' @export
compute_all_parameter_profiles <- function(params_current,
                                           cost_to_optimize,
                                           bounds,
                                           profile_options = list(),
                                           optimizer = "optim",
                                           verbose = TRUE,
                                           n_cores = 1,
                                           ...) {

  # Validate inputs
  if (!is.numeric(params_current)) {
    stop("params_current must be a numeric vector")
  }

  if (!is.function(cost_to_optimize)) {
    stop("cost_to_optimize must be a function")
  }

  if (!is.list(bounds) || !all(c("lower", "upper") %in% names(bounds))) {
    stop("bounds must be a list with 'lower' and 'upper' elements")
  }

  n_params <- length(params_current)

  if (length(bounds$lower) != n_params || length(bounds$upper) != n_params) {
    stop("bounds dimensions must match parameter vector length")
  }

  # Set default profile options
  default_options <- list(
    grid_method = "uniform",
    grid_points = 25,
    max_grid_range_multiplier = 2.0,
    ll_ratio_threshold = 3.84  # 95% CI for chi-square with 1 df
  )

  profile_options <- modifyList(default_options, profile_options)

  # Validate optimizer
  optimizer_info <- validate_and_setup_optimizer(optimizer, ...)

  if (verbose) {
    cat(sprintf("Computing profile likelihood for %d parameters using %s grid...\n",
                n_params, profile_options$grid_method))
    cat(sprintf("Grid points per parameter: %d\n", profile_options$grid_points))
    cat(sprintf("Using optimizer: %s\n", optimizer_info$name))
    if (n_params > 1 && n_cores > 1) {
      cat(sprintf("Using parallel computing with %d cores...\n", n_cores))
    }
  }

  # Initialize output arrays
  profiles <- vector("list", n_params)
  confidence_intervals <- matrix(NA, nrow = n_params, ncol = 2)
  rownames(confidence_intervals) <- names(params_current) %||% paste0("param_", 1:n_params)
  colnames(confidence_intervals) <- c("2.5%", "97.5%")

  # Record computation start time
  start_time <- Sys.time()

  # Compute profiles for all parameters
  if (n_cores > 1 && n_params > 1) {
    # Parallel computation
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("compute_parameter_profile", "create_parameter_grid",
                                  "create_uniform_grid", "create_adaptive_grid",
                                  "optimize_conditional", "evaluate_conditional_cost",
                                  "extract_confidence_interval", "find_confidence_bounds",
                                  "cost_to_optimize", "bounds", "profile_options",
                                  "optimizer_info", "params_current"),
                            envir = environment())

    # Load required packages on cluster
    parallel::clusterEvalQ(cl, {
      if (exists("optimizer_info") && optimizer_info$name == "deoptim") {
        library(DEoptim)
      }
    })

    # Run profiling in parallel
    results <- parallel::parLapply(cl, 1:n_params, function(i) {
      compute_parameter_profile(
        param_index = i,
        params_current = params_current,
        cost_to_optimize = cost_to_optimize,
        bounds = bounds,
        profile_options = profile_options,
        optimizer_info = optimizer_info,
        verbose = FALSE  # Silent in parallel
      )
    })

    # Extract results
    for (i in 1:n_params) {
      profiles[[i]] <- results[[i]]$profile
      confidence_intervals[i, ] <- results[[i]]$confidence_interval
    }

  } else {
    # Sequential computation
    for (i in 1:n_params) {
      if (verbose) {
        cat(sprintf("Computing profile for parameter %d/%d...\n", i, n_params))
      }

      result <- compute_parameter_profile(
        param_index = i,
        params_current = params_current,
        cost_to_optimize = cost_to_optimize,
        bounds = bounds,
        profile_options = profile_options,
        optimizer_info = optimizer_info,
        verbose = verbose
      )

      profiles[[i]] <- result$profile
      confidence_intervals[i, ] <- result$confidence_interval
    }
  }

  if (verbose) {
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("Profile likelihood computation completed in %.2f seconds.\n", elapsed_time))

    # Report success rate
    successful_profiles <- sum(sapply(profiles, function(p) !isTRUE(p$failed)))
    cat(sprintf("Successfully computed %d/%d parameter profiles.\n",
                successful_profiles, n_params))

    # Issue consolidated warnings
    issue_consolidated_warnings(profiles, n_params)
  }

  # Validate results
  validate_profile_results(profiles, confidence_intervals, n_params, verbose)

  # Create summary
  summary_df <- create_profile_summary(profiles, confidence_intervals, params_current)

  if (verbose) {
    cat("Profile likelihood uncertainty quantification completed successfully.\n")
  }

  return(list(
    profiles = profiles,
    confidence_intervals = confidence_intervals,
    summary = summary_df,
    options = profile_options,
    optimizer = optimizer_info$name
  ))
}

#' Validate and Setup Optimizer
#'
#' @param optimizer Character string or function specifying optimizer
#' @param ... Additional arguments for optimizer
#' @return List with optimizer information
#' @keywords internal
validate_and_setup_optimizer <- function(optimizer, ...) {

  if (is.character(optimizer)) {
    if (optimizer == "optim") {
      return(list(
        name = "optim",
        type = "builtin",
        extra_args = list(...)
      ))
    } else if (optimizer == "deoptim") {
      # Check if DEoptim is available
      if (!requireNamespace("DEoptim", quietly = TRUE)) {
        stop("DEoptim package required but not available. Install with: install.packages('DEoptim')")
      }
      return(list(
        name = "deoptim",
        type = "builtin",
        extra_args = list(...)
      ))
    } else {
      stop(sprintf("Unknown optimizer: '%s'. Use 'optim', 'deoptim', or provide a custom function.", optimizer))
    }
  } else if (is.function(optimizer)) {
    # Validate custom optimizer signature
    arg_names <- names(formals(optimizer))
    required_args <- c("fn", "par", "lower", "upper")

    if (!all(required_args %in% arg_names)) {
      stop(sprintf("Custom optimizer must have signature: function(fn, par, lower, upper, ...) \nMissing arguments: %s",
                   paste(setdiff(required_args, arg_names), collapse = ", ")))
    }

    return(list(
      name = "custom",
      type = "function",
      func = optimizer,
      extra_args = list(...)
    ))
  } else {
    stop("optimizer must be a character string ('optim', 'deoptim') or a function")
  }
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Compute Profile for Single Parameter
#'
#' @param param_index Index of parameter to profile
#' @param params_current Current parameter values
#' @param cost_to_optimize Cost function to optimize
#' @param bounds Parameter bounds
#' @param profile_options Profiling options
#' @param optimizer_info Optimizer configuration
#' @param verbose Logical for verbose output
#' @return List with profile data and confidence interval
#' @keywords internal
compute_parameter_profile <- function(param_index, params_current, cost_to_optimize,
                                      bounds, profile_options, optimizer_info, verbose = FALSE) {

  tryCatch({
    # Create grid for this parameter
    grid_values <- create_parameter_grid(
      param_index = param_index,
      params_current = params_current,
      bounds = bounds,
      profile_options = profile_options
    )

    n_grid <- length(grid_values)
    profile_costs <- numeric(n_grid)
    optimizer_exit_flags <- character(n_grid)

    if (verbose) {
      cat(sprintf("  Evaluating %d grid points...\n", n_grid))
    }

    # Evaluate cost at each grid point
    for (j in 1:n_grid) {
      result <- optimize_conditional(
        param_index = param_index,
        fixed_value = grid_values[j],
        params_current = params_current,
        cost_to_optimize = cost_to_optimize,
        bounds = bounds,
        optimizer_info = optimizer_info
      )

      profile_costs[j] <- result$cost
      optimizer_exit_flags[j] <- result$exit_flag
    }

    # Create profile data structure
    profile_data <- list(
      param_index = param_index,
      grid_values = grid_values,
      profile_costs = profile_costs,
      optimizer_exit_flags = optimizer_exit_flags,
      optimal_cost = min(profile_costs, na.rm = TRUE),
      failed = FALSE
    )

    # Extract confidence interval
    confidence_interval <- extract_confidence_interval(
      profile_data = profile_data,
      profile_options = profile_options
    )

    return(list(
      profile = profile_data,
      confidence_interval = confidence_interval
    ))

  }, error = function(e) {
    warning(sprintf("Failed to compute profile for parameter %d: %s", param_index, e$message))

    profile_data <- list(
      param_index = param_index,
      grid_values = numeric(0),
      profile_costs = numeric(0),
      optimizer_exit_flags = character(0),
      optimal_cost = NA,
      failed = TRUE,
      error_message = e$message
    )

    return(list(
      profile = profile_data,
      confidence_interval = c(NA, NA)
    ))
  })
}

#' Create Parameter Grid
#'
#' @param param_index Index of parameter to profile
#' @param params_current Current parameter values
#' @param bounds Parameter bounds
#' @param profile_options Profiling options
#' @return Numeric vector of grid values
#' @keywords internal
create_parameter_grid <- function(param_index, params_current, bounds, profile_options) {

  current_value <- params_current[param_index]
  lower_bound <- bounds$lower[param_index]
  upper_bound <- bounds$upper[param_index]

  # Determine grid range
  param_range <- upper_bound - lower_bound
  max_range <- profile_options$max_grid_range_multiplier * param_range

  # Center grid around current value, but respect bounds
  grid_lower <- max(lower_bound, current_value - max_range/2)
  grid_upper <- min(upper_bound, current_value + max_range/2)

  # Adjust if one side hits a bound
  if (grid_lower == lower_bound && grid_upper < upper_bound) {
    grid_upper <- min(upper_bound, grid_lower + max_range)
  } else if (grid_upper == upper_bound && grid_lower > lower_bound) {
    grid_lower <- max(lower_bound, grid_upper - max_range)
  }

  if (profile_options$grid_method == "uniform") {
    return(create_uniform_grid(grid_lower, grid_upper, profile_options$grid_points))
  } else if (profile_options$grid_method == "adaptive") {
    return(create_adaptive_grid(grid_lower, grid_upper, current_value, profile_options$grid_points))
  } else {
    stop(sprintf("Unknown grid method: %s", profile_options$grid_method))
  }
}

#' Create Uniform Grid
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param n_points Number of grid points
#' @return Numeric vector of uniformly spaced grid values
#' @keywords internal
create_uniform_grid <- function(grid_lower, grid_upper, n_points) {
  return(seq(grid_lower, grid_upper, length.out = n_points))
}

#' Create Adaptive Grid
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param center_value Central value for denser sampling
#' @param n_points Number of grid points
#' @return Numeric vector of adaptively spaced grid values
#' @keywords internal
create_adaptive_grid <- function(grid_lower, grid_upper, center_value, n_points) {

  # Create quadratic spacing for denser sampling near center
  t <- seq(-1, 1, length.out = n_points)

  # Quadratic transformation: denser near center (t=0)
  spacing_factor <- 0.7  # Controls density near center
  adjusted_t <- sign(t) * abs(t)^spacing_factor

  # Map to parameter range
  grid_values <- grid_lower + (adjusted_t + 1) / 2 * (grid_upper - grid_lower)

  # Ensure center value is included
  center_idx <- which.min(abs(grid_values - center_value))
  grid_values[center_idx] <- center_value

  return(sort(grid_values))
}

#' Optimize with Fixed Parameter
#'
#' @param param_index Index of parameter to fix
#' @param fixed_value Value to fix parameter at
#' @param params_current Current parameter values
#' @param cost_to_optimize Cost function to optimize
#' @param bounds Parameter bounds
#' @param optimizer_info Optimizer configuration
#' @return List with optimized cost and exit flag
#' @keywords internal
optimize_conditional <- function(param_index, fixed_value, params_current,
                                 cost_to_optimize, bounds, optimizer_info) {

  # Create conditional cost function
  conditional_cost <- function(free_params, ...) {
    evaluate_conditional_cost(free_params, param_index, fixed_value,
                              params_current, cost_to_optimize, ...)
  }

  # Set up free parameters and bounds
  free_indices <- setdiff(1:length(params_current), param_index)
  initial_free_params <- params_current[free_indices]
  lower_free <- bounds$lower[free_indices]
  upper_free <- bounds$upper[free_indices]

  # Handle single parameter case
  if (length(free_indices) == 0) {
    cost <- do.call(cost_to_optimize, c(list(params_current), optimizer_info$extra_args))
    return(list(cost = cost, exit_flag = "direct_evaluation"))
  }

  # Optimize based on optimizer type
  if (optimizer_info$type == "builtin") {
    if (optimizer_info$name == "optim") {

      extra_args <- optimizer_info$extra_args

      result <- do.call(optim, c(
        list(
          par = initial_free_params,
          fn = conditional_cost,
          lower = lower_free,
          upper = upper_free
        ),
        extra_args
      ))

      cost <- result$value
      exit_flag <- if (result$convergence == 0) "success" else paste("convergence_", result$convergence, sep="")

    } else if (optimizer_info$name == "deoptim") {

      result <- do.call(DEoptim::DEoptim, c(
        list(
          fn = conditional_cost,
          lower = lower_free,
          upper = upper_free
        ),
        optimizer_info$extra_args
      ))

      cost <- result$optim$bestval
      exit_flag <- "success"  # DEoptim doesn't provide convergence codes
    }

  } else if (optimizer_info$type == "function") {

    result <- do.call(optimizer_info$func, c(
      list(
        fn = conditional_cost,
        par = initial_free_params,
        lower = lower_free,
        upper = upper_free
      ),
      optimizer_info$extra_args
    ))

    # Extract cost and exit flag based on expected structure
    if (is.list(result)) {
      cost <- result$value %||% result$minimum %||% result$cost %||% result$objective
      exit_flag <- if (!is.null(result$convergence)) {
        if (result$convergence == 0) "success" else paste("convergence_", result$convergence, sep="")
      } else {
        "success"
      }
    } else {
      cost <- result
      exit_flag <- "success"
    }
  }

  return(list(cost = cost, exit_flag = exit_flag))
}

#' Evaluate Conditional Cost Function
#'
#' @param free_params Free parameter values
#' @param param_index Index of fixed parameter
#' @param fixed_value Value of fixed parameter
#' @param params_current Template parameter vector
#' @param cost_to_optimize Original cost function
#' @param extra_args Additional arguments to pass to cost function
#' @return Scalar cost value
#' @keywords internal
evaluate_conditional_cost <- function(free_params, param_index, fixed_value,
                                      params_current, cost_to_optimize, ...) {

  # Reconstruct full parameter vector
  full_params <- params_current
  full_params[param_index] <- fixed_value

  free_indices <- setdiff(1:length(params_current), param_index)
  full_params[free_indices] <- free_params

  # Call cost function with additional arguments
  cost_to_optimize(full_params, ...)
}

#' Extract Confidence Interval from Profile
#'
#' @param profile_data Profile data structure
#' @param profile_options Profiling options
#' @return Numeric vector [lower, upper] confidence bounds
#' @keywords internal
extract_confidence_interval <- function(profile_data, profile_options) {

  if (profile_data$failed || length(profile_data$profile_costs) == 0) {
    return(c(NA, NA))
  }

  # Calculate likelihood ratio threshold
  threshold_cost <- profile_data$optimal_cost + profile_options$ll_ratio_threshold

  # Find confidence bounds
  confidence_bounds <- find_confidence_bounds(
    grid_values = profile_data$grid_values,
    profile_costs = profile_data$profile_costs,
    threshold_cost = threshold_cost
  )

  return(confidence_bounds)
}

#' Find Confidence Bounds by Interpolation
#'
#' @param grid_values Grid parameter values
#' @param profile_costs Profile cost values
#' @param threshold_cost Threshold cost for confidence level
#' @return Numeric vector [lower, upper] confidence bounds
#' @keywords internal
find_confidence_bounds <- function(grid_values, profile_costs, threshold_cost) {

  if (length(grid_values) < 2) {
    return(c(NA, NA))
  }

  # Sort by grid values
  sorted_indices <- order(grid_values)
  sorted_grid <- grid_values[sorted_indices]
  sorted_costs <- profile_costs[sorted_indices]

  # Find crossings of threshold
  below_threshold <- sorted_costs <= threshold_cost

  # Find lower bound (first crossing from left)
  lower_bound <- NA
  for (i in 1:(length(sorted_costs) - 1)) {
    if (below_threshold[i] != below_threshold[i + 1]) {
      # Linear interpolation
      if (below_threshold[i] && !below_threshold[i + 1]) {
        # Crossing from below to above threshold (lower bound)
        t <- (threshold_cost - sorted_costs[i]) / (sorted_costs[i + 1] - sorted_costs[i])
        lower_bound <- sorted_grid[i] + t * (sorted_grid[i + 1] - sorted_grid[i])
        break
      }
    }
  }

  # Find upper bound (last crossing from right)
  upper_bound <- NA
  for (i in length(sorted_costs):2) {
    if (below_threshold[i] != below_threshold[i - 1]) {
      # Linear interpolation
      if (below_threshold[i] && !below_threshold[i - 1]) {
        # Crossing from above to below threshold (upper bound)
        t <- (threshold_cost - sorted_costs[i - 1]) / (sorted_costs[i] - sorted_costs[i - 1])
        upper_bound <- sorted_grid[i - 1] + t * (sorted_grid[i] - sorted_grid[i - 1])
        break
      }
    }
  }

  # If no crossings found, use grid extremes where cost is below threshold
  if (is.na(lower_bound)) {
    valid_indices <- which(below_threshold)
    if (length(valid_indices) > 0) {
      lower_bound <- min(sorted_grid[valid_indices])
    }
  }

  if (is.na(upper_bound)) {
    valid_indices <- which(below_threshold)
    if (length(valid_indices) > 0) {
      upper_bound <- max(sorted_grid[valid_indices])
    }
  }

  return(c(lower_bound, upper_bound))
}

#' Issue Consolidated Warnings
#'
#' @param profiles List of profile results
#' @param n_params Number of parameters
#' @keywords internal
issue_consolidated_warnings <- function(profiles, n_params) {

  # Count failures and convergence issues
  failed_params <- which(sapply(profiles, function(p) isTRUE(p$failed)))

  convergence_issues <- list()
  for (i in 1:n_params) {
    if (!profiles[[i]]$failed) {
      flags <- profiles[[i]]$optimizer_exit_flags
      non_success <- flags[flags != "success" & flags != "direct_evaluation"]
      if (length(non_success) > 0) {
        convergence_issues[[length(convergence_issues) + 1]] <- list(
          param = i,
          issues = unique(non_success)
        )
      }
    }
  }

  # Issue warnings
  if (length(failed_params) > 0) {
    warning(sprintf("Profile computation failed for parameter(s): %s",
                    paste(failed_params, collapse = ", ")), call. = FALSE)
  }

  if (length(convergence_issues) > 0) {
    for (issue in convergence_issues) {
      warning(sprintf("Convergence issues for parameter %d: %s",
                      issue$param, paste(issue$issues, collapse = ", ")), call. = FALSE)
    }
  }
}

#' Validate Profile Results
#'
#' @param profiles List of profile results
#' @param confidence_intervals Matrix of confidence intervals
#' @param n_params Number of parameters
#' @param verbose Logical for verbose output
#' @keywords internal
validate_profile_results <- function(profiles, confidence_intervals, n_params, verbose) {

  # Check for any valid results
  valid_profiles <- sum(!sapply(profiles, function(p) isTRUE(p$failed)))

  if (valid_profiles == 0) {
    stop("All parameter profiles failed to compute. Check cost function and bounds.")
  }

  if (verbose && valid_profiles < n_params) {
    cat(sprintf("Warning: Only %d/%d profiles computed successfully.\n",
                valid_profiles, n_params))
  }

  # Check confidence interval validity
  invalid_ci <- is.na(confidence_intervals[, 1]) | is.na(confidence_intervals[, 2])

  if (verbose && any(invalid_ci)) {
    invalid_params <- which(invalid_ci)
    cat(sprintf("Warning: Confidence intervals could not be determined for parameter(s): %s\n",
                paste(invalid_params, collapse = ", ")))
  }
}

#' Create Profile Summary
#'
#' @param profiles List of profile results
#' @param confidence_intervals Matrix of confidence intervals
#' @param params_current Current parameter values
#' @return Data frame with summary statistics
#' @keywords internal
create_profile_summary <- function(profiles, confidence_intervals, params_current) {

  n_params <- length(profiles)
  param_names <- names(params_current) %||% paste0("param_", 1:n_params)

  summary_df <- data.frame(
    parameter = param_names,
    current_value = params_current,
    lower_ci = confidence_intervals[, 1],
    upper_ci = confidence_intervals[, 2],
    ci_width = confidence_intervals[, 2] - confidence_intervals[, 1],
    stringsAsFactors = FALSE
  )

  # Add status information
  summary_df$status <- sapply(1:n_params, function(i) {
    if (profiles[[i]]$failed) {
      "failed"
    } else if (is.na(confidence_intervals[i, 1]) || is.na(confidence_intervals[i, 2])) {
      "no_ci"
    } else {
      "success"
    }
  })

  return(summary_df)
}
