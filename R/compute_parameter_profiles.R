#' Compute Profile Likelihood for All Parameters
#'
#' Compute profile likelihood confidence intervals for all parameters in
#' parallel using a grid-based approach with customizable optimization methods.
#'
#' @param params_current Numeric vector of parameters at optimum
#' @param negLogLikelihood Function that takes parameter vector and returns
#'   scalar cost
#' @param bounds List with elements 'lower' and 'upper' containing parameter
#'   bounds
#' @param profile_options List of profiling options (see Details)
#' @param optimizer Character string or function specifying optimizer (default:
#'   "optim")
#' @param optim_options List of optimizer-specific arguments (default: list())
#' @param verbose Logical indicating whether to print progress messages
#' @param cluster Optional cluster object for parallel computation (default:
#'   NULL)
#' @param ... Additional arguments passed to the `negLogLikelihood` function
#'
#' @details The \code{profile_options} list can contain:
#' \itemize{
#'   \item \code{grid_method}: "uniform" or "adaptive" (default: "adaptive").
#'   "uniform" uses evenly spaced grid points, "adaptive" adjusts the density of grid points
#'   to be more dense near the optimum using a quadratic spacing.
#'   \item \code{grid_points}: Number of grid points per parameter (default: 100)
#'   \item \code{max_grid_range_multiplier}: Controls the width of the profiling grid as a
#'   fraction of the total parameter range defined by bounds. The grid is centered around
#'   the current parameter value and extends max_range/2 in each direction, where
#'   max_range = max_grid_range_multiplier * (upper_bound - lower_bound). A value of 1.0
#'   covers the full parameter range for a parameter in the center of the range,
#'   while 2.0 (default) attempts to cover twice the
#'   range but is constrained by the bounds, this ensures that the whole range is covered even if the
#'   parameter is near the boundary.
#'   \item \code{ll_ratio_threshold}: Likelihood ratio threshold for CI (default: 3.84)
#' }
#'
#'   The \code{optimizer} can be:
#' \itemize{
#'   \item \code{"optim"}: Use base R optim() with options via optim_options
#'   \item \code{"deoptim"}: Use DEoptim package (must be installed)
#'   \item A function with signature: function(fn, par, lower, upper, ...)
#' }
#'
#'   The \code{optim_options} list can contain optimizer-specific arguments:
#' \itemize{
#'   \item For optim: \code{method}, \code{control}, etc.
#'   \item For deoptim: \code{itermax}, \code{NP}, \code{trace}, etc. Default values are: `itermax = 100, trace = FALSE`.
#'     Other values are the defaults set by \code{\link[DEoptim]{DEoptim.control}}.
#'   \item For custom optimizers: any arguments the optimizer function accepts
#' }
#'
#' @return List containing:
#' \itemize{
#'   \item \code{profiles}: List of profile data for each parameter
#'   \item \code{confidence_intervals}: Matrix of \eqn{[lower, upper]} bounds for each parameter
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
#'   optim_options = list(
#'     method = "L-BFGS-B",
#'     control = list(maxit = 100)
#'   ),
#'   verbose = FALSE  # Suppress output for clean example
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
#'   optim_options = list(
#'     itermax = 100,
#'     NP = 20
#'   )
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
                                     optim_options = list(),
                                     verbose = TRUE,
                                     cluster = NULL,
                                     ...) {

  # Validate inputs
  if (!is.numeric(params_current)) {
    stop("Error in computeLikelihoodProfiles: params_current must be a numeric vector")
  }

  if (!is.function(negLogLikelihood)) {
    stop("Error in computeLikelihoodProfiles: negLogLikelihood must be a function")
  }

  if (!is.list(bounds) || !all(c("lower", "upper") %in% names(bounds))) {
    stop("Error in computeLikelihoodProfiles: bounds must be a list with 'lower' and 'upper' elements")
  }

  n_params = length(params_current)

  if (length(bounds$lower) != n_params || length(bounds$upper) != n_params) {
    stop("Error in computeLikelihoodProfiles: bounds dimensions must match parameter vector length")
  }

  if (!is.numeric(bounds$lower) || !is.numeric(bounds$upper)) {
    stop("Error in computeLikelihoodProfiles: bounds$lower and bounds$upper must be numeric vectors")
  }
  if (any(bounds$lower >= bounds$upper)) {
    stop("Error in computeLikelihoodProfiles: bounds$lower must be less than bounds$upper for all elements")
  }
  if (any(params_current < bounds$lower) || any(params_current > bounds$upper)) {
    stop("Error in computeLikelihoodProfiles: params_current must be within the specified bounds")
  }

  if (is.null(names(params_current))) {
    names(params_current) = paste0("param_", 1:n_params)
  }

  # Set default profile options
  default_options = list(
    grid_method = "adaptive",
    grid_points = 100,
    max_grid_range_multiplier = 2.0,
    ll_ratio_threshold = 3.84  # 95% CI for chi-square with 1 df
  )

  profile_options = utils::modifyList(default_options, profile_options)

  # Validate profile options
  validateProfileOptions(profile_options)

  # Validate optimizer
  optimizer_info = validateAndSetupOptimizer(optimizer, optim_options, bounds)

  # Capture likelihood function arguments
  likelihood_args = list(...)

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
  rownames(confidence_intervals) = names(params_current)
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
        likelihood_args = likelihood_args,
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
        likelihood_args = likelihood_args,
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
#' @param optim_options List of optimizer-specific arguments
#' @param bounds Optional list with elements 'lower' and 'upper' containing parameter bounds
#' @return List with optimizer information
#' @keywords internal
validateAndSetupOptimizer = function(optimizer, optim_options = list(), bounds = NULL) {

  # Validate bounds if provided
  if (!is.null(bounds)) {
    if (!is.list(bounds) || !all(c("lower", "upper") %in% names(bounds))) {
      stop("Error in validateAndSetupOptimizer: bounds must be a list with 'lower' and 'upper' elements")
    }

    if (!is.numeric(bounds$lower) || !is.numeric(bounds$upper)) {
      stop("Error in validateAndSetupOptimizer: bounds$lower and bounds$upper must be numeric vectors")
    }

    if (length(bounds$lower) != length(bounds$upper)) {
      stop("Error in validateAndSetupOptimizer: bounds$lower and bounds$upper must have the same length")
    }

    if (any(bounds$lower >= bounds$upper)) {
      stop("Error in validateAndSetupOptimizer: bounds$lower must be less than bounds$upper for all elements")
    }
  }

  if (is.character(optimizer)) {
    if (optimizer == "optim") {
      # Handle bounds for optim
      if (!is.null(bounds)) {
        # Check if method is already specified
        current_method = optim_options$method
        if (!is.null(current_method) && current_method != "L-BFGS-B") {
          warning(sprintf("Warning in validateAndSetupOptimizer: Method '%s' specified for optim with bounds. Changing to 'L-BFGS-B' to handle bounds properly.", current_method), call. = FALSE)
        }
        # Set method to L-BFGS-B for bound-constrained optimization
        optim_options$method = "L-BFGS-B"
      }

      return(list(
        name = "optim",
        type = "builtin",
        extra_args = optim_options
      ))
    } else if (optimizer == "deoptim") {
      # Check if DEoptim is available
      if (!requireNamespace("DEoptim", quietly = TRUE)) {
        stop("Error in validateAndSetupOptimizer: DEoptim package required but not available. Install with: install.packages('DEoptim')")
      }
      default_optim_options = DEoptim::DEoptim.control(itermax = 100, trace = FALSE)
      optim_options = utils::modifyList(default_optim_options, optim_options)
      return(list(
        name = "deoptim",
        type = "builtin",
        extra_args = optim_options
      ))
    } else {
      stop(sprintf("Error in validateAndSetupOptimizer: Unknown optimizer: '%s'. Use 'optim', 'deoptim', or provide a custom function.", optimizer))
    }
  } else if (is.function(optimizer)) {
    # Validate custom optimizer signature
    arg_names = names(formals(optimizer))
    required_args = c("fn", "par", "lower", "upper")

    if (!all(required_args %in% arg_names)) {
      stop(sprintf("Error in validateAndSetupOptimizer: Custom optimizer must have signature: function(fn, par, lower, upper, ...) \nMissing arguments: %s",
                   paste(setdiff(required_args, arg_names), collapse = ", ")))
    }

    return(list(
      name = "custom",
      type = "function",
      func = optimizer,
      extra_args = optim_options
    ))
  } else {
    stop("Error in validateAndSetupOptimizer: optimizer must be a character string ('optim', 'deoptim') or a function")
  }
}

#' Validate Profile Options
#'
#' @param profile_options List of profiling options to validate
#' @keywords internal
validateProfileOptions = function(profile_options) {

  # Validate grid_method
  valid_grid_methods = c("uniform", "adaptive")
  if (!profile_options$grid_method %in% valid_grid_methods) {
    stop(sprintf("Error in validateProfileOptions: profile_options$grid_method must be one of: %s. Got: '%s'",
                 paste(valid_grid_methods, collapse = ", "), profile_options$grid_method))
  }

  # Validate grid_points
  if (!is.numeric(profile_options$grid_points) || length(profile_options$grid_points) != 1) {
    stop("Error in validateProfileOptions: profile_options$grid_points must be a single numeric value")
  }
  if (profile_options$grid_points < 3) {
    stop("Error in validateProfileOptions: profile_options$grid_points must be at least 3")
  }
  if (profile_options$grid_points != round(profile_options$grid_points)) {
    stop("Error in validateProfileOptions: profile_options$grid_points must be an integer")
  }

  # Validate max_grid_range_multiplier
  if (!is.numeric(profile_options$max_grid_range_multiplier) ||
      length(profile_options$max_grid_range_multiplier) != 1) {
    stop("Error in validateProfileOptions: profile_options$max_grid_range_multiplier must be a single numeric value")
  }
  if (profile_options$max_grid_range_multiplier <= 0) {
    stop("Error in validateProfileOptions: profile_options$max_grid_range_multiplier must be positive")
  }

  # Validate ll_ratio_threshold
  if (!is.numeric(profile_options$ll_ratio_threshold) ||
      length(profile_options$ll_ratio_threshold) != 1) {
    stop("Error in validateProfileOptions: profile_options$ll_ratio_threshold must be a single numeric value")
  }
  if (profile_options$ll_ratio_threshold <= 0) {
    stop("Error in validateProfileOptions: profile_options$ll_ratio_threshold must be positive")
  }

  # Issue warnings for potentially problematic values
  if (profile_options$grid_points > 1000) {
    warning(sprintf("Warning in validateProfileOptions: profile_options$grid_points is very large (%d). This may result in long computation times.",
                    profile_options$grid_points), call. = FALSE)
  }

  if (profile_options$max_grid_range_multiplier > 2) {
    warning("Warning in validateProfileOptions: It is not useful to set profile_options$max_grid_range_multiplier greater than 2.0, as this already covers the whole parameter range.",
            call. = FALSE)
  }
}

#' Optimize Grid Direction
#'
#' Helper function to optimize along a specific direction from the optimal parameter.
#' Handles the common optimization logic for both left and right sides.
#'
#' @param grid_indices Vector of grid indices to optimize (e.g., optimal_index:n_grid for right side)
#' @param param_index Index of parameter to profile (1-based index)
#' @param grid_values Numeric vector of grid parameter values to evaluate
#' @param warm_start_params Numeric vector of initial warm start parameter values
#' @param negLogLikelihood Function that takes parameter vector and returns scalar cost
#' @param bounds List with elements 'lower' and 'upper' containing parameter bounds
#' @param optimizer_info List containing optimizer configuration from validateAndSetupOptimizer
#' @param profile_costs Numeric vector to store profile costs (modified in place)
#' @param optimizer_exit_flags Character vector to store optimizer exit flags (modified in place)
#' @param optimal_params_matrix Matrix to store optimal parameters for each grid point (modified in place)
#' @param optimization_success Logical vector to track optimization success (modified in place)
#' @param likelihood_args List of arguments to pass to negLogLikelihood function
#' @return List with updated arrays: profile_costs, optimizer_exit_flags, optimal_params_matrix, optimization_success, and success_count
#' @keywords internal
optimizeGridDirection = function(grid_indices, param_index, grid_values, warm_start_params,
                                 negLogLikelihood, bounds, optimizer_info, likelihood_args,
                                 profile_costs, optimizer_exit_flags, optimal_params_matrix,
                                 optimization_success) {

  current_warm_start = warm_start_params
  success_count = 0

  for (j in grid_indices) {
    if (j < 1 || j > length(grid_values)){
      stop("Error in optimizeGridDirection: grid index out of bounds")
    }
    browser()

    tryCatch({
      result = optimizeConditional(
        param_index = param_index,
        fixed_value = grid_values[j],
        warm_start_params = current_warm_start,
        negLogLikelihood = negLogLikelihood,
        bounds = bounds,
        optimizer_info = optimizer_info,
        likelihood_args = likelihood_args
      )

      profile_costs[j] = result$cost
      optimizer_exit_flags[j] = result$exit_flag
      optimal_params_matrix[j, ] = result$optimal_params

      # Track success for exit flags
      is_success = (result$exit_flag == "success" || result$exit_flag == "direct_evaluation")
      optimization_success[j] = is_success
      if (is_success) success_count = success_count + 1

      # Update warm start for next iteration (warm start strategy)
      if (is_success) {
        current_warm_start = result$optimal_params
      }
      # On failure, keep using previous warm start parameters

    }, error = function(e) {
      warning(sprintf("Warning in optimizeGridDirection: Optimization failed for parameter %d at grid point %d: %s",
                      param_index, j, e$message))
      profile_costs[j] = NA
      optimizer_exit_flags[j] = "error"
      optimization_success[j] = FALSE
      optimal_params_matrix[j, ] = current_warm_start  # Use previous warm start
      # Don't update warm start on error
    })
  }

  return(list(
    profile_costs = profile_costs,
    optimizer_exit_flags = optimizer_exit_flags,
    optimal_params_matrix = optimal_params_matrix,
    optimization_success = optimization_success,
    success_count = success_count
  ))
}

#' Compute Profile for Single Parameter
#'
#' @param param_index Index of parameter to profile
#' @param params_current Current parameter values
#' @param negLogLikelihood Cost function to optimize
#' @param bounds Parameter bounds
#' @param profile_options Profiling options
#' @param optimizer_info Optimizer configuration
#' @param likelihood_args List of arguments to pass to negLogLikelihood function
#' @param verbose Logical for verbose output
#' @return List with profile data and confidence interval
#' @keywords internal
computeParameterProfile = function(param_index, params_current, negLogLikelihood,
                                   bounds, profile_options, optimizer_info,
                                   likelihood_args = list(), verbose = FALSE) {

  tryCatch({
    # Create grid for this parameter
    grid_result = createParameterGrid(
      param_index = param_index,
      params_current = params_current,
      bounds = bounds,
      profile_options = profile_options
    )

    grid_values = grid_result$grid_values

    # Initialize exit flags structure
    exit_flags = list(
      low_success_rate = FALSE,
      confidence_interval_failed = FALSE,
      optimal_value_out_of_bounds = FALSE,
      grid_issue = FALSE,
      grid_too_narrow = FALSE
    )
    if (!is.null(grid_result$optimal_value_out_of_bounds)) {
      exit_flags$optimal_value_out_of_bounds = grid_result$optimal_value_out_of_bounds
    }
    if (!is.null(grid_result$grid_issue)) {
      exit_flags$grid_issue = grid_result$grid_issue
    }

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

    # Initialize success tracking
    optimization_success = logical(n_grid)
    success_count = 0

    # Optimize right side (from optimal outward)
    right_result = optimizeGridDirection(
      grid_indices = optimal_index:n_grid,
      param_index = param_index,
      grid_values = grid_values,
      warm_start_params = params_current,
      negLogLikelihood = negLogLikelihood,
      bounds = bounds,
      optimizer_info = optimizer_info,
      likelihood_args = likelihood_args,
      profile_costs = profile_costs,
      optimizer_exit_flags = optimizer_exit_flags,
      optimal_params_matrix = optimal_params_matrix,
      optimization_success = optimization_success
    )

    # Update arrays with right side results
    profile_costs = right_result$profile_costs
    optimizer_exit_flags = right_result$optimizer_exit_flags
    optimal_params_matrix = right_result$optimal_params_matrix
    optimization_success = right_result$optimization_success
    success_count = success_count + right_result$success_count

    # Optimize left side (from optimal outward) - reset warm start
    # Only happens if the optimal index is not the first grid point
    # i.e. there is a left side to explore
    if (optimal_index > 1) {
      left_result = optimizeGridDirection(
        grid_indices = seq(optimal_index-1, 1),
        param_index = param_index,
        grid_values = grid_values,
        warm_start_params = params_current,  # Reset to original optimal params
        negLogLikelihood = negLogLikelihood,
        bounds = bounds,
        optimizer_info = optimizer_info,
        likelihood_args = likelihood_args,
        profile_costs = profile_costs,
        optimizer_exit_flags = optimizer_exit_flags,
        optimal_params_matrix = optimal_params_matrix,
        optimization_success = optimization_success
      )

      # Update arrays with left side results
      profile_costs = left_result$profile_costs
      optimizer_exit_flags = left_result$optimizer_exit_flags
      optimal_params_matrix = left_result$optimal_params_matrix
      optimization_success = left_result$optimization_success
      success_count = success_count + left_result$success_count
    }

    # Calculate success rate and set exit flags
    success_rate = success_count / n_grid
    if (success_rate < 0.5) {
      exit_flags$low_success_rate = TRUE
    }
    # Check if corresponding cost is below threshold. If not set current params at this index
    # to the minimum observed cost with a warning.
    min_cost_index = which.min(profile_costs)
    min_observed_cost = profile_costs[min_cost_index]
    optimal_grid_cost = profile_costs[optimal_index]

    # If the cost at optimal parameter position is significantly higher than minimum found
    ll_ratio_optimal = 2 * (optimal_grid_cost - min_observed_cost)
    ll_threshold = profile_options$ll_ratio_threshold
    if (!is.na(optimal_grid_cost) && ll_ratio_optimal >= ll_threshold) {
      warning(sprintf("Warning in computeLikelihoodProfiles: For parameter %d, likelihood ratio test at provided optimal parameter (%.4f) exceeds threshold (%.4f). Using best parameters found during profiling.",
                      param_index, ll_ratio_optimal, ll_threshold), call. = FALSE)

      # Update the profile data to use the best parameters found
      optimal_param_value = grid_values[min_cost_index]
    } else {
      optimal_param_value = params_current[param_index]
    }

    # Create profile data structure
    profile_data = list(
      param_index = param_index,
      optimal_param_value = optimal_param_value,
      grid_values = grid_values,
      profile_costs = profile_costs,
      optimizer_exit_flags = optimizer_exit_flags,
      optimal_params_matrix = optimal_params_matrix,
      optimization_success = optimization_success,
      success_rate = success_rate,
      exit_flags = exit_flags,
      optimal_cost = min_observed_cost,
      failed = FALSE
    )

    # Extract confidence interval
    confidence_interval = extractConfidenceInterval(
      profile_data = profile_data,
      profile_options = profile_options
    )

    # Update exit flags based on CI extraction
    if (any(is.na(confidence_interval))) {
      profile_data$exit_flags$confidence_interval_failed = TRUE
    }

    return(list(
      profile = profile_data,
      confidence_interval = confidence_interval
    ))

  }, error = function(e) {
    warning(sprintf("Warning in computeParameterProfile: Failed to compute profile for parameter %d: %s", param_index, e$message))

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
#' @return List with grid_values and optional flags
#' @keywords internal
createParameterGrid = function(param_index, params_current, bounds, profile_options) {

  current_value = params_current[param_index]
  lower_bound = bounds$lower[param_index]
  upper_bound = bounds$upper[param_index]

  # Check if optimal value is within bounds
  optimal_value_out_of_bounds = current_value < lower_bound || current_value > upper_bound

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
    grid_values = createUniformGrid(grid_lower, grid_upper, current_value, profile_options$grid_points)
  } else if (profile_options$grid_method == "adaptive") {
    grid_values = createAdaptiveGrid(grid_lower, grid_upper, current_value, profile_options$grid_points)
  } else {
    stop(sprintf("Error in createParameterGrid: Unknown grid method: %s", profile_options$grid_method))
  }

  # Check for grid issues
  grid_issue = length(grid_values) < 3

  list(
    grid_values = grid_values,
    optimal_value_out_of_bounds = optimal_value_out_of_bounds,
    grid_issue = grid_issue
  )
}

#' Create Uniform Grid
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param n_points Number of grid points
#' @return Numeric vector of uniformly spaced grid values
#' @keywords internal
createUniformGrid = function(grid_lower, grid_upper, current_value, n_points) {
  res_grid = c(seq(grid_lower, current_value, length.out = floor(n_points/2) + 1),
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
#' @param likelihood_args List of arguments to pass to negLogLikelihood function
#' @return List with optimized cost, exit flag, and optimal parameters
#' @keywords internal
optimizeConditional = function(param_index, fixed_value, warm_start_params,
                               negLogLikelihood, bounds, optimizer_info, likelihood_args = list()) {

  # Create conditional cost function
  conditional_cost = function(free_params) {
    evaluateConditionalCost(free_params, param_index, fixed_value,
                            warm_start_params, negLogLikelihood, likelihood_args)
  }

  # Set up free parameters and bounds
  free_indices = setdiff(1:length(warm_start_params), param_index)
  initial_free_params = warm_start_params[free_indices]
  lower_free = bounds$lower[free_indices]
  upper_free = bounds$upper[free_indices]

  # Handle single parameter case
  if (length(free_indices) == 0) {
    # Construct full parameter vector with fixed parameter
    full_params = warm_start_params
    full_params[param_index] = fixed_value
    cost = do.call(negLogLikelihood, c(list(full_params), likelihood_args))
    return(list(cost = cost, exit_flag = "direct_evaluation", optimal_params = full_params))
  }

  # Optimize based on optimizer type
  if (optimizer_info$type == "builtin") {
    if (optimizer_info$name == "optim") {
      result = do.call(stats::optim, c(
        list(
          par = initial_free_params,
          fn = conditional_cost,
          lower = lower_free,
          upper = upper_free
        ),
        optimizer_info$extra_args
      ))

      cost = result$value
      exit_flag = if (result$convergence == 0) "success" else paste0("convergence_", result$convergence)

    } else if (optimizer_info$name == "deoptim") {
      # Check that we have valid bounds for DEoptim
      if (length(lower_free) == 0 || length(upper_free) == 0) {
        stop("Error in optimizeConditional: DEoptim requires at least one free parameter")
      }
      if (length(lower_free) != length(upper_free)) {
        stop("Error in optimizeConditional: DEoptim lower and upper bounds must have same length")
      }
      if (any(is.na(lower_free)) || any(is.na(upper_free))) {
        stop("Error in optimizeConditional: DEoptim bounds cannot contain NA values")
      }
      if (any(lower_free >= upper_free)) {
        stop("Error in optimizeConditional: DEoptim requires lower < upper for all parameters")
      }
      browser()

      result = DEoptim::DEoptim(fn = conditional_cost,
                                lower = lower_free,
                                upper = upper_free,
                                control = optimizer_info$extra_args)

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

    # Extract cost, exit flag based on expected structure
    if (is.list(result)) {
      cost = result$value
      exit_flag = if (!is.null(result$convergence)) {
        if (result$convergence == 0) "success" else paste0("convergence_", result$convergence)
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
      optimal_params[free_indices] = result$par
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
#' @param likelihood_args List of arguments to pass to negLogLikelihood function
#' @return Scalar cost value
#' @keywords internal
evaluateConditionalCost = function(free_params, param_index, fixed_value,
                                   warm_start_params, negLogLikelihood, likelihood_args = list()) {
  # Reconstruct full parameter vector
  full_params = warm_start_params
  full_params[param_index] = fixed_value

  free_indices = setdiff(1:length(warm_start_params), param_index)
  full_params[free_indices] = free_params

  # Call cost function with likelihood arguments only
  print(full_params)
  res = do.call(negLogLikelihood, c(list(full_params), likelihood_args))
  print(res)
  res
}

#' Extract Confidence Interval from Profile
#'
#' @param profile_data Profile data structure
#' @param profile_options Profiling options
#' @return Numeric vector \eqn{[lower, upper]} confidence bounds
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
#' @return Numeric vector \eqn{[lower, upper]} confidence bounds
#' @keywords internal
findConfidenceBounds = function(grid_values, profile_costs, threshold_cost, optimal_param_value) {

  if (length(grid_values) < 2) {
    return(c(NA, NA))
  }

  # Sort by grid values
  sorted_indices = order(grid_values)
  sorted_grid = grid_values[sorted_indices]
  sorted_costs = profile_costs[sorted_indices]

  n = length(sorted_grid)

  # Convert costs to likelihood ratios
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

  # Count failures and other issues
  failed_params = which(sapply(profiles, function(p) isTRUE(p$failed)))
  low_success_rate_params = c()
  ci_failed_params = c()
  bounds_issue_params = c()
  grid_issue_params = c()

  # Analyze exit flags
  for (i in 1:n_params) {
    if (!profiles[[i]]$failed && !is.null(profiles[[i]]$exit_flags)) {
      flags = profiles[[i]]$exit_flags
      if (isTRUE(flags$low_success_rate)) low_success_rate_params = c(low_success_rate_params, i)
      if (isTRUE(flags$confidence_interval_failed)) ci_failed_params = c(ci_failed_params, i)
      if (isTRUE(flags$optimal_value_out_of_bounds)) bounds_issue_params = c(bounds_issue_params, i)
      if (isTRUE(flags$grid_issue)) grid_issue_params = c(grid_issue_params, i)
    }
  }

  # Issue warnings
  if (length(failed_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Profile computation failed for parameter(s): %s",
                    paste(failed_params, collapse = ", ")), call. = FALSE)
  }

  if (length(low_success_rate_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Low optimization success rate (<50%%) for parameter(s): %s",
                    paste(low_success_rate_params, collapse = ", ")), call. = FALSE)
  }

  if (length(ci_failed_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Failed to extract confidence intervals for parameter(s): %s",
                    paste(ci_failed_params, collapse = ", ")), call. = FALSE)
  }

  if (length(bounds_issue_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Optimal values outside bounds for parameter(s): %s",
                    paste(bounds_issue_params, collapse = ", ")), call. = FALSE)
  }

  if (length(grid_issue_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Grid coverage issues for parameter(s): %s",
                    paste(grid_issue_params, collapse = ", ")), call. = FALSE)
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
    stop("Error in computeLikelihoodProfiles: All parameter profiles failed to compute. Check cost function and bounds.")
  }

  if (verbose && valid_profiles < n_params) {
    warning(sprintf("Warning in validateProfileResults: Only %d/%d profiles computed successfully.\n",
                    valid_profiles, n_params))
  }

  # Check confidence interval validity and quality
  invalid_ci = is.na(confidence_intervals[, 1]) | is.na(confidence_intervals[, 2])

  if (verbose && any(invalid_ci)) {
    invalid_params = which(invalid_ci)
    warning(sprintf("Warning in validateProfileResults: Confidence intervals could not be determined for parameter(s): %s\n",
                    paste(invalid_params, collapse = ", ")))
  }

  # Check for CI quality issues
  if (verbose) {
    valid_ci = !invalid_ci
    if (sum(valid_ci) > 0) {
      ci_widths = confidence_intervals[valid_ci, 2] - confidence_intervals[valid_ci, 1]
      extremely_wide = sum(ci_widths > 1000, na.rm = TRUE)
      extremely_narrow = sum(ci_widths < 1e-10, na.rm = TRUE)

      if (extremely_wide > 0) {
        warning(sprintf("Warning in validateProfileResults: %d parameters have very wide confidence intervals (>1000).\n", extremely_wide))
      }
      if (extremely_narrow > 0) {
        warning(sprintf("Warning in validateProfileResults: %d parameters have very narrow confidence intervals (<1e-10).\n", extremely_narrow))
      }
    }
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
