#' Compute Profile Likelihood for All Parameters
#'
#' Compute profile likelihood confidence intervals for all parameters in parallel
#' using a grid-based approach with customizable optimization methods.
#'
#' @param params_current Numeric vector of parameters at optimum
#' @param negLogLikelihood Function that takes parameter vector and returns scalar cost
#' @param bounds List with elements 'lower' and 'upper' containing parameter bounds
#' @param profile_options List of profiling options (see Details)
#' @param optimizer Character string or function specifying optimizer (default: "optim")
#' @param verbose Logical indicating whether to print progress messages
#' @param cluster Optional cluster object for parallel computation (default: NULL)
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
#' # Basic example with simple quadratic function
#'
#' # Define a simple quadratic likelihood function
#' quadratic_likelihood = function(params) {
#'   x = params[1]
#'   y = params[2]
#'
#'   # Simple quadratic bowl with minimum at (2, 1)
#'   cost = (x - 2)^2 + 2 * (y - 1)^2
#'
#'   return(cost)
#' }
#'
#' # Set up the problem
#' true_params = c(2.0, 1.0)
#' names(true_params) = c("param_x", "param_y")
#' bounds = list(lower = c(0, 0), upper = c(4, 3))
#'
#' # Compute profile likelihood with small grid for speed
#' result = computeLikelihoodProfiles(
#'   params_current = true_params,
#'   negLogLikelihood = quadratic_likelihood,
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
#' result_deoptim = computeLikelihoodProfiles(
#'   params_current = true_params,
#'   negLogLikelihood = quadratic_likelihood,
#'   bounds = bounds,
#'   optimizer = "deoptim",
#'   itermax = 100,
#'   NP = 20
#' )
#'
#' # Custom optimizer example
#' my_optimizer = function(fn, par, lower, upper, ...) {
#'   optim(par = par, fn = fn, method = "Nelder-Mead", ...)
#' }
#'
#' result_custom = computeLikelihoodProfiles(
#'   params_current = true_params,
#'   negLogLikelihood = quadratic_likelihood,
#'   bounds = bounds,
#'   optimizer = my_optimizer
#' )
#' }
#'
#' @export
computeLikelihoodProfiles = function(params_current,
                                     negLogLikelihood,
                                     bounds,
                                     profile_options = list(),
                                     optimizer = "optim",
                                     verbose = TRUE,
                                     cluster = NULL,
                                     ...) {

  # Validate inputs
  if (!is.numeric(params_current)) {
    stop("params_current must be a numeric vector")
  }

  if (!is.function(negLogLikelihood)) {
    stop("negLogLikelihood must be a function")
  }

  if (!is.list(bounds) || !all(c("lower", "upper") %in% names(bounds))) {
    stop("bounds must be a list with 'lower' and 'upper' elements")
  }

  n_params = length(params_current)

  if (length(bounds$lower) != n_params || length(bounds$upper) != n_params) {
    stop("bounds dimensions must match parameter vector length")
  }

  # Set default profile options
  default_options = list(
    grid_method = "uniform",
    grid_points = 25,
    max_grid_range_multiplier = 2.0,
    ll_ratio_threshold = 3.84  # 95% CI for chi-square with 1 df
  )

  profile_options = modifyList(default_options, profile_options)

  # Validate optimizer
  optimizer_info = validateAndSetupOptimizer(optimizer, ...)

  if (verbose) {
    cat(sprintf("Computing profile likelihood for %d parameters using %s grid...\n",
                n_params, profile_options$grid_method))
    cat(sprintf("Grid points per parameter: %d\n", profile_options$grid_points))
    cat(sprintf("Using optimizer: %s\n", optimizer_info$name))
    if (!is.null(cluster)) {
      cat("Using parallel computation with provided cluster...\n")
    }
  }

  # Initialize output arrays
  profiles = vector("list", n_params)
  confidence_intervals = matrix(NA, nrow = n_params, ncol = 2)
  rownames(confidence_intervals) = names(params_current) %||% paste0("param_", 1:n_params)
  colnames(confidence_intervals) = c("2.5%", "97.5%")

  # Record computation start time
  start_time = Sys.time()

  # Compute profiles for all parameters
  if (!is.null(cluster)) {
    # Load required packages on cluster
    parallel::clusterEvalQ(cluster, {
      if (exists("optimizer_info") && optimizer_info$name == "deoptim") {
        library(DEoptim)
      }
    })

    # Run profiling in parallel
    results = parallel::parLapply(cluster, 1:n_params, function(i) {
      computeParameterProfile(
        param_index = i,
        params_current = params_current,
        negLogLikelihood = negLogLikelihood,
        bounds = bounds,
        profile_options = profile_options,
        optimizer_info = optimizer_info,
        verbose = FALSE  # Silent in parallel
      )
    })

    # Extract results
    for (i in 1:n_params) {
      profiles[[i]] = results[[i]]$profile
      confidence_intervals[i, ] = results[[i]]$confidence_interval
    }

  } else {
    # Sequential computation
    for (i in 1:n_params) {
      if (verbose) {
        cat(sprintf("Computing profile for parameter %d/%d...\n", i, n_params))
      }

      result = computeParameterProfile(
        param_index = i,
        params_current = params_current,
        negLogLikelihood = negLogLikelihood,
        bounds = bounds,
        profile_options = profile_options,
        optimizer_info = optimizer_info,
        verbose = verbose
      )

      profiles[[i]] = result$profile
      confidence_intervals[i, ] = result$confidence_interval
    }
  }

  if (verbose) {
    elapsed_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("Profile likelihood computation completed in %.2f seconds.\n", elapsed_time))

    # Report success rate
    successful_profiles = sum(sapply(profiles, function(p) !isTRUE(p$failed)))
    cat(sprintf("Successfully computed %d/%d parameter profiles.\n",
                successful_profiles, n_params))

    # Issue consolidated warnings
    issueConsolidatedWarnings(profiles, n_params)
  }

  # Validate results
  validateProfileResults(profiles, confidence_intervals, n_params, verbose)

  # Create summary
  summary_df = createProfileSummary(profiles, confidence_intervals, params_current)

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
validateAndSetupOptimizer = function(optimizer, ...) {

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
    arg_names = names(formals(optimizer))
    required_args = c("fn", "par", "lower", "upper")

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
`%||%` = function(x, y) if (is.null(x)) y else x

#' Compute Profile for Single Parameter
#'
#' @param param_index Index of parameter to profile
#' @param params_current Current parameter values
#' @param negLogLikelihood Cost function to optimize
#' @param bounds Parameter bounds
#' @param profile_options Profiling options
#' @param optimizer_info Optimizer configuration
#' @param verbose Logical for verbose output
#' @return List with profile data and confidence interval
#' @keywords internal
computeParameterProfile = function(param_index, params_current, negLogLikelihood,
                                   bounds, profile_options, optimizer_info, verbose = FALSE) {

  tryCatch({
    # Create grid for this parameter
    grid_values = createParameterGrid(
      param_index = param_index,
      params_current = params_current,
      bounds = bounds,
      profile_options = profile_options
    )

    n_grid = length(grid_values)
    profile_costs = numeric(n_grid)
    optimizer_exit_flags = character(n_grid)
    optimal_params_matrix = matrix(NA, nrow = n_grid, ncol = length(params_current))

    if (verbose) {
      cat(sprintf("  Evaluating %d grid points...\n", n_grid))
    }

    # Find the index of the optimal parameter value in the grid
    optimal_index = which.min(abs(grid_values - params_current[param_index]))

    # Initialize warm start parameters
    warm_start_params = params_current

    # Evaluate cost at each grid point - Right side (from optimal outward)
    for (j in optimal_index:n_grid) {
      tryCatch({
        result = optimizeConditional(
          param_index = param_index,
          fixed_value = grid_values[j],
          warm_start_params = warm_start_params,
          negLogLikelihood = negLogLikelihood,
          bounds = bounds,
          optimizer_info = optimizer_info
        )

        profile_costs[j] = result$cost
        optimizer_exit_flags[j] = result$exit_flag
        optimal_params_matrix[j, ] = result$optimal_params

        # Update warm start for next iteration (warm start strategy)
        if (result$exit_flag == "success" || result$exit_flag == "direct_evaluation") {
          warm_start_params = result$optimal_params
        }
        # On failure, keep using previous warm start parameters

      }, error = function(e) {
        warning(sprintf("Optimization failed for parameter %d at grid point %d: %s",
                       param_index, j, e$message))
        profile_costs[j] = NA
        optimizer_exit_flags[j] = "error"
        optimal_params_matrix[j, ] = warm_start_params  # Use previous warm start
        # Don't update warm start on error
      })
    }

    # Reset warm start for left side optimization
    warm_start_params = params_current

    # Evaluate cost at each grid point - Left side (from optimal outward)
    for (j in (optimal_index-1):1) {
      if (j < 1) break  # Safety check

      tryCatch({
        result = optimizeConditional(
          param_index = param_index,
          fixed_value = grid_values[j],
          warm_start_params = warm_start_params,
          negLogLikelihood = negLogLikelihood,
          bounds = bounds,
          optimizer_info = optimizer_info
        )

        profile_costs[j] = result$cost
        optimizer_exit_flags[j] = result$exit_flag
        optimal_params_matrix[j, ] = result$optimal_params

        # Update warm start for next iteration (warm start strategy)
        if (result$exit_flag == "success" || result$exit_flag == "direct_evaluation") {
          warm_start_params = result$optimal_params
        }
        # On failure, keep using previous warm start parameters

      }, error = function(e) {
        warning(sprintf("Optimization failed for parameter %d at grid point %d: %s",
                       param_index, j, e$message))
        profile_costs[j] = NA
        optimizer_exit_flags[j] = "error"
        optimal_params_matrix[j, ] = warm_start_params  # Use previous warm start
        # Don't update warm start on error
      })
    }

    # Create profile data structure
    profile_data = list(
      param_index = param_index,
      optimal_param_value = params_current[param_index],
      grid_values = grid_values,
      profile_costs = profile_costs,
      optimizer_exit_flags = optimizer_exit_flags,
      optimal_params_matrix = optimal_params_matrix,
      optimal_cost = min(profile_costs, na.rm = TRUE),
      failed = FALSE
    )

    # Extract confidence interval
    confidence_interval = extractConfidenceInterval(
      profile_data = profile_data,
      profile_options = profile_options
    )

    return(list(
      profile = profile_data,
      confidence_interval = confidence_interval
    ))

  }, error = function(e) {
    warning(sprintf("Failed to compute profile for parameter %d: %s", param_index, e$message))

    profile_data = list(
      param_index = param_index,
      optimal_param_value = params_current[param_index],
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
createParameterGrid = function(param_index, params_current, bounds, profile_options) {

  current_value = params_current[param_index]
  lower_bound = bounds$lower[param_index]
  upper_bound = bounds$upper[param_index]

  # Determine grid range
  param_range = upper_bound - lower_bound
  max_range = profile_options$max_grid_range_multiplier * param_range

  # Center grid around current value, but respect bounds
  grid_lower = max(lower_bound, current_value - max_range/2)
  grid_upper = min(upper_bound, current_value + max_range/2)

  # Adjust if one side hits a bound
  if (grid_lower == lower_bound && grid_upper < upper_bound) {
    grid_upper = min(upper_bound, grid_lower + max_range)
  } else if (grid_upper == upper_bound && grid_lower > lower_bound) {
    grid_lower = max(lower_bound, grid_upper - max_range)
  }

  if (profile_options$grid_method == "uniform") {
    return(createUniformGrid(grid_lower, grid_upper, current_value, profile_options$grid_points))
  } else if (profile_options$grid_method == "adaptive") {
    return(createAdaptiveGrid(grid_lower, grid_upper, current_value, profile_options$grid_points))
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
createUniformGrid = function(grid_lower, grid_upper, current_value, n_points) {
  res_grid = c(seq(grid_lower, current_value, length.out = floor(n_points/2)) + 1,
               seq(current_value, grid_upper, length.out = floor(n_points/2) + 1))

  unique(res_grid)
}

#' Create Adaptive Grid
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param center_value Central value for denser sampling
#' @param n_points Number of grid points
#' @return Numeric vector of adaptively spaced grid values
#' @keywords internal
createAdaptiveGrid = function(grid_lower, grid_upper, center_value, n_points) {

  # Create quadratic spacing for denser sampling near center
  n_half = floor(n_points/2) + 1
  t_onesided = seq(0, 1, length.out = n_half)^2
  left_side = center_value - t_onesided*(center_value - grid_lower)
  right_side = center_value + t_onesided*(grid_upper - center_value)

  res_grid = c(left_side, right_side)
  unique(sort(res_grid))
}

#' Optimize with Fixed Parameter
#'
#' @param param_index Index of parameter to fix
#' @param fixed_value Value to fix parameter at
#' @param warm_start_params Warm start parameter values for optimization
#' @param negLogLikelihood Cost function to optimize
#' @param bounds Parameter bounds
#' @param optimizer_info Optimizer configuration
#' @return List with optimized cost, exit flag, and optimal parameters
#' @keywords internal
optimizeConditional = function(param_index, fixed_value, warm_start_params,
                               negLogLikelihood, bounds, optimizer_info) {

  # Create conditional cost function
  conditional_cost = function(free_params, ...) {
    evaluateConditionalCost(free_params, param_index, fixed_value,
                            warm_start_params, negLogLikelihood, ...)
  }

  # Set up free parameters and bounds
  free_indices = setdiff(1:length(warm_start_params), param_index)
  initial_free_params = warm_start_params[free_indices]
  lower_free = bounds$lower[free_indices]
  upper_free = bounds$upper[free_indices]

  # Handle single parameter case
  if (length(free_indices) == 0) {
    browser()
    cost = do.call(negLogLikelihood, c(list(warm_start_params), optimizer_info$extra_args))
    return(list(cost = cost, exit_flag = "direct_evaluation", optimal_params = warm_start_params))
  }

  # Optimize based on optimizer type
  if (optimizer_info$type == "builtin") {
    if (optimizer_info$name == "optim") {

      extra_args = optimizer_info$extra_args

      result = do.call(optim, c(
        list(
          par = initial_free_params,
          fn = conditional_cost,
          lower = lower_free,
          upper = upper_free
        ),
        extra_args
      ))

      cost = result$value
      exit_flag = if (result$convergence == 0) "success" else paste("convergence_", result$convergence, sep="")

    } else if (optimizer_info$name == "deoptim") {

      result = do.call(DEoptim::DEoptim, c(
        list(
          fn = conditional_cost,
          lower = lower_free,
          upper = upper_free
        ),
        optimizer_info$extra_args
      ))

      cost = result$optim$bestval
      exit_flag = "success"  # DEoptim doesn't provide convergence codes
    }

  } else if (optimizer_info$type == "function") {

    result = do.call(optimizer_info$func, c(
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
      cost = result$value %||% result$minimum %||% result$cost %||% result$objective
      exit_flag = if (!is.null(result$convergence)) {
        if (result$convergence == 0) "success" else paste("convergence_", result$convergence, sep="")
      } else {
        "success"
      }
    } else {
      cost = result
      exit_flag = "success"
    }
  }

  # Reconstruct full parameter vector with optimized free parameters
  optimal_params = warm_start_params
  if (optimizer_info$name == "optim") {
    optimal_params[free_indices] = result$par
  } else if (optimizer_info$name == "deoptim") {
    optimal_params[free_indices] = result$optim$bestmem
  } else if (optimizer_info$type == "function") {
    if (is.list(result)) {
      optimal_params[free_indices] = result$par %||% result$minimum %||% result$solution %||% result$x
    }
  }
  optimal_params[param_index] = fixed_value  # Ensure fixed parameter stays fixed

  return(list(cost = cost, exit_flag = exit_flag, optimal_params = optimal_params))
}

#' Evaluate Conditional Cost Function
#'
#' @param free_params Free parameter values
#' @param param_index Index of fixed parameter
#' @param fixed_value Value of fixed parameter
#' @param warm_start_params Template parameter vector
#' @param negLogLikelihood Original cost function
#' @param extra_args Additional arguments to pass to cost function
#' @return Scalar cost value
#' @keywords internal
evaluateConditionalCost = function(free_params, param_index, fixed_value,
                                   warm_start_params, negLogLikelihood, ...) {

  # Reconstruct full parameter vector
  full_params = warm_start_params
  full_params[param_index] = fixed_value

  free_indices = setdiff(1:length(warm_start_params), param_index)
  full_params[free_indices] = free_params

  # Call cost function with additional arguments
  negLogLikelihood(full_params, ...)
}

#' Extract Confidence Interval from Profile
#'
#' @param profile_data Profile data structure
#' @param profile_options Profiling options
#' @param optimal_param_value Optimal parameter value from params_current
#' @return Numeric vector [lower, upper] confidence bounds
#' @keywords internal
extractConfidenceInterval = function(profile_data, profile_options) {

  if (profile_data$failed || length(profile_data$profile_costs) == 0) {
    return(c(NA, NA))
  }

  # Calculate likelihood ratio threshold
  threshold_cost = profile_data$optimal_cost + profile_options$ll_ratio_threshold

  # Find confidence bounds
  confidence_bounds = findConfidenceBounds(
    grid_values = profile_data$grid_values,
    profile_costs = profile_data$profile_costs,
    threshold_cost = threshold_cost,
    optimal_param_value = profile_data$optimal_param_value
  )

  return(confidence_bounds)
}

#' Find Confidence Bounds by Interpolation
#'
#' @param grid_values Grid parameter values
#' @param profile_costs Profile cost values
#' @param threshold_cost Threshold cost for confidence level
#' @param optimal_param_value Optimal parameter value from params_current
#' @return Numeric vector [lower, upper] confidence bounds
#' @keywords internal
findConfidenceBounds = function(grid_values, profile_costs, threshold_cost, optimal_param_value) {

  if (length(grid_values) < 2) {
    return(c(NA, NA))
  }

  # browser()

  # Sort by grid values
  sorted_indices = order(grid_values)
  sorted_grid = grid_values[sorted_indices]
  sorted_costs = profile_costs[sorted_indices]

  n = length(sorted_grid)

  # Convert costs to likelihood ratios (following MATLAB approach)
  ll_ratios = 2 * (sorted_costs - min(sorted_costs, na.rm = TRUE))

  # Find crossings where likelihood ratio crosses threshold
  crossings = c()

  for (i in 1:(n-1)) {
    # Check if threshold is crossed between points i and i+1
    if ((ll_ratios[i] <= threshold_cost && ll_ratios[i+1] > threshold_cost) ||
        (ll_ratios[i] > threshold_cost && ll_ratios[i+1] <= threshold_cost)) {

      # Linear interpolation to find exact crossing point
      if (abs(ll_ratios[i+1] - ll_ratios[i]) > 1e-10) {
        t = (threshold_cost - ll_ratios[i]) / (ll_ratios[i+1] - ll_ratios[i])
        crossing_point = sorted_grid[i] + t * (sorted_grid[i+1] - sorted_grid[i])
        crossings = c(crossings, crossing_point)
      }
    }
  }

  # Handle edge cases
  if (length(crossings) == 0) {
    # No crossings found - check if all points are within threshold
    within_threshold = ll_ratios <= threshold_cost
    if (sum(within_threshold) == 0) {
      return(c(NA, NA))
    } else if (all(within_threshold)) {
      # All points within threshold - use grid extremes
      return(c(min(sorted_grid), max(sorted_grid)))
    } else {
      # Some points within threshold - use their extremes
      return(c(min(sorted_grid[within_threshold]), max(sorted_grid[within_threshold])))
    }
  }

  # Find regions within threshold
  within_threshold = ll_ratios <= threshold_cost

  if (sum(within_threshold) == 0) {
    return(c(NA, NA))
  }

  # Determine bounds based on crossings and threshold regions
  sorted_crossings = sort(crossings)

  # Strategy: find the largest continuous region within threshold
  # that contains the minimum likelihood point (optimal parameter)
  lower_bound = min(sorted_grid[within_threshold])
  upper_bound = max(sorted_grid[within_threshold])

  # Refine bounds using crossings if they exist
  if (length(sorted_crossings) > 0) {
    # Find crossings that bracket the optimal parameter value
    lower_crossings = sorted_crossings[sorted_crossings <= optimal_param_value]
    upper_crossings = sorted_crossings[sorted_crossings >= optimal_param_value]

    if (length(lower_crossings) > 0) {
      lower_bound = max(lower_crossings)
    }

    if (length(upper_crossings) > 0) {
      upper_bound = min(upper_crossings)
    }
  }

  # Final validation
  if (lower_bound > upper_bound) {
    temp = lower_bound
    lower_bound = upper_bound
    upper_bound = temp
  }

  return(c(lower_bound, upper_bound))
}

#' Issue Consolidated Warnings
#'
#' @param profiles List of profile results
#' @param n_params Number of parameters
#' @keywords internal
issueConsolidatedWarnings = function(profiles, n_params) {

  # Count failures and convergence issues
  failed_params = which(sapply(profiles, function(p) isTRUE(p$failed)))

  convergence_issues = list()
  for (i in 1:n_params) {
    if (!profiles[[i]]$failed) {
      flags = profiles[[i]]$optimizer_exit_flags
      non_success = flags[flags != "success" & flags != "direct_evaluation"]
      if (length(non_success) > 0) {
        convergence_issues[[length(convergence_issues) + 1]] = list(
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
validateProfileResults = function(profiles, confidence_intervals, n_params, verbose) {

  # Check for any valid results
  valid_profiles = sum(!sapply(profiles, function(p) isTRUE(p$failed)))

  if (valid_profiles == 0) {
    stop("All parameter profiles failed to compute. Check cost function and bounds.")
  }

  if (verbose && valid_profiles < n_params) {
    cat(sprintf("Warning: Only %d/%d profiles computed successfully.\n",
                valid_profiles, n_params))
  }

  # Check confidence interval validity
  invalid_ci = is.na(confidence_intervals[, 1]) | is.na(confidence_intervals[, 2])

  if (verbose && any(invalid_ci)) {
    invalid_params = which(invalid_ci)
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
createProfileSummary = function(profiles, confidence_intervals, params_current) {

  n_params = length(profiles)
  param_names = names(params_current) %||% paste0("param_", 1:n_params)

  summary_df = data.frame(
    parameter = param_names,
    current_value = params_current,
    lower_ci = confidence_intervals[, 1],
    upper_ci = confidence_intervals[, 2],
    ci_width = confidence_intervals[, 2] - confidence_intervals[, 1],
    stringsAsFactors = FALSE
  )

  # Add status information
  summary_df$status = sapply(1:n_params, function(i) {
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
