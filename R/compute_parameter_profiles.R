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
#' @param nested_cluster_call Optional expression to create a nested cluster for
#'   DEoptim parallel computation within each parameter profile (default: NULL).
#'   Only used when optimizer is `"deoptim"`. The expression should create an
#'   object named `nested_cluster` when evaluated. This enables nested
#'   parallelization: outer parallelization across parameters, inner
#'   parallelization within \code{\link[DEoptim]{DEoptim}}.
#' @param generateData Function to generate bootstrap datasets. Required when threshold_method is "bootstrap".
#'   Should take a parameter vector and return a dataset that can be passed to negLogLikelihood.
#'   The function signature should be: `function(params) { ... return(dataset) }`
#' @param ... Additional arguments passed to the `negLogLikelihood` function
#'
#' @details The \code{profile_options} list can contain:
#' \itemize{
#'   \item \code{grid_method}: "linear", "quadratic", or "adaptive" (default: "quadratic").
#'   "linear" creates linear spacing on each side of the optimum (may have different spacing
#'   on left/right sides if optimum is off-center). "quadratic" uses quadratic transformation
#'   for denser sampling near the optimum with progressively wider spacing away from it.
#'   "adaptive" uses iterative refinement to focus grid points near confidence interval boundaries.
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
#'   \item \code{threshold_method}: "fixed" or "bootstrap" (default: "fixed").
#'   When "fixed", uses \code{ll_ratio_threshold} directly. When "bootstrap",
#'   uses empirical thresholds from bootstrap datasets.
#'   \item \code{bootstrap_conf_level}: Confidence level for bootstrap thresholds (default: 0.95)
#'   \item \code{n_bootstrap}: Number of bootstrap samples for empirical threshold computation (default: 500)
#'   \item \code{adaptive_max_iterations}: Maximum refinement iterations for adaptive grids (default: 3)
#'   \item \code{adaptive_initial_fraction}: Fraction of total points for initial grid (default: 0.3)
#'   \item \code{adaptive_refinement_fraction}: Fraction of remaining points per iteration (default: 0.5)
#'
#'
#' }
#'
#'   The \code{optimizer} can be:
#' \itemize{
#'   \item \code{"optim"}: Use base R \code{\link[stats]{optim}} with options via optim_options
#'   \item \code{"deoptim"}: Use \code{\link[DEoptim]{DEoptim}} (must be installed)
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
#' @section Nested Parallel Computation: For DEoptim optimization, nested
#'   parallelization can be enabled using the \code{nested_cluster_call}
#'   parameter. This allows:
#' \itemize{
#'   \item Outer parallelization: Multiple parameter profiles computed in parallel using \code{cluster}
#'   \item Inner parallelization: DEoptim uses parallel computation within each profile using a nested cluster
#' }
#'
#'   The \code{nested_cluster_call} should be an expression (e.g., created with
#'   \code{quote()}) that, when evaluated, creates an object named
#'   \code{nested_cluster}. For example:
#' \code{nested_cluster_call = quote({nested_cluster = parallel::makeCluster(2)
#'                                    parallel::clusterEvalQ(nested_cluster, library(some_package))
#'                                   })}.
#'
#'   The object `nested_cluster` will be passed to
#'   \code{\link[DEoptim]{DEoptim.control}} as the \code{cluster} argument.
#'
#'   \strong{Important:} Ensure the nested cluster uses fewer cores than
#'   available to avoid resource conflicts with the outer cluster. The nested
#'   cluster is automatically passed to DEoptim's cluster parameter for parallel
#'   function evaluations.
#'
#' @return List containing:
#' \itemize{
#'   \item \code{profiles}: List of profile data for each parameter
#'   \item \code{confidence_intervals}: Matrix of \eqn{[lower, upper]} bounds for each parameter
#'   \item \code{summary}: Data frame summarizing results
#'   \item \code{options}: Profile options used in computation
#'   \item \code{optimizer}: Name of optimizer used
#'   \item \code{bootstrap_stats}: 3D array of bootstrap likelihood ratio statistics
#'   \eqn{[n_bootstrap, n_params, n_grid]} (only when \code{threshold_method = "bootstrap"})
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
#'     grid_method = "quadratic",
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
#' @md
computeLikelihoodProfiles = function(params_current,
                                     negLogLikelihood,
                                     bounds,
                                     profile_options = list(),
                                     optimizer = "optim",
                                     optim_options = list(),
                                     verbose = TRUE,
                                     cluster = NULL,
                                     nested_cluster_call = NULL,
                                     generateData = NULL,  # New parameter
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
    grid_method = "quadratic",
    grid_points = 100,
    max_grid_range_multiplier = 2.0,
    ll_ratio_threshold = stats::qchisq(0.95, 1),  # 95% CI for chi-square with 1 df
    threshold_method = "fixed",
    bootstrap_conf_level = 0.95,
    n_bootstrap = 500,
    adaptive_max_iterations = 3,
    adaptive_initial_fraction = 0.3,
    adaptive_refinement_fraction = 0.5
  )

  profile_options = utils::modifyList(default_options, profile_options)

  # Validate profile options
  validateProfileOptions(profile_options, generateData, negLogLikelihood)

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
  if (verbose) {
    if (!is.null(cluster)) {
      cat("Computing profiles in parallel...\n")
    } else {
      cat("Computing profiles sequentially...\n")
    }
  }

  # Use locLapply for conditional parallel/sequential execution
  results = locLapply(cluster, seq_len(n_params), function(i) {
    if (verbose && is.null(cluster)) {
      cat(sprintf("Computing profile for parameter %d/%d...\n", i, n_params))
    }

    computeParameterProfile(
      param_index = i,
      params_current = params_current,
      negLogLikelihood = negLogLikelihood,
      bounds = bounds,
      profile_options = profile_options,
      optimizer_info = optimizer_info,
      likelihood_args = likelihood_args,
      nested_cluster_call = nested_cluster_call,
      verbose = if (is.null(cluster)) verbose else FALSE,
      generateData = generateData
    )
  })

  # Extract results
  # Initialize bootstrap statistics array if using bootstrap method
  if (profile_options$threshold_method == "bootstrap") {
    n_grid_points = dim(results[[1]]$bootstrap_stats)[2] # Take the size of the first parameter's grid, assuming all grids have the same size
    bootstrap_stats = array(NA, dim = c(profile_options$n_bootstrap, n_params, n_grid_points))
  } else {
    bootstrap_stats = NULL
  }
  for (i_param in seq_len(n_params)) {
    profiles[[i_param]] = results[[i_param]]$profile
    confidence_intervals[i_param, ] = results[[i_param]]$confidence_interval

    # Collect bootstrap statistics
    if (profile_options$threshold_method == "bootstrap" && !is.null(results[[i_param]]$bootstrap_stats)) {
      bootstrap_stats[, i_param, ] = results[[i_param]]$bootstrap_stats
    }
  }

  if (verbose) {
    elapsed_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("Profile likelihood computation completed in %.2f seconds.\n", elapsed_time))

    # Report success rate
    successful_profiles = sum(sapply(profiles, function(p) !isTRUE(p$failed)))
    cat(sprintf("Successfully computed %d/%d parameter profiles.\n",
                successful_profiles, n_params))

  }
  # Issue consolidated warnings
  issueConsolidatedWarnings(profiles, n_params)

  # Validate results
  validateProfileResults(profiles, confidence_intervals, n_params, verbose)

  # Create summary
  summary_df = createProfileSummary(profiles, confidence_intervals, params_current)

  if (verbose) {
    cat("Profile likelihood uncertainty quantification completed successfully.\n")
  }

  list(
    profiles = profiles,
    confidence_intervals = confidence_intervals,
    summary = summary_df,
    options = profile_options,
    optimizer = optimizer_info$name,
    bootstrap_stats = bootstrap_stats
  )
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
        stop("Error in validateAndSetupOptimizer: DEoptim package required for 'deoptim' method but not available. Install with: install.packages('DEoptim')")
      }
      if (!requireNamespace("lhs", quietly = TRUE)) {
        stop("Error in validateAndSetupOptimizer: lhs package required for 'deoptim' method but not available. Install with: install.packages('lhs')")
      }
      default_optim_options = list(itermax = 100, trace = FALSE)
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
#' @param generateData Function to generate bootstrap datasets (required for bootstrap method)
#' @param negLogLikelihood Negative log-likelihood function (required for bootstrap validation)
#' @keywords internal
validateProfileOptions = function(profile_options, generateData = NULL, negLogLikelihood = NULL) {

  # Validate grid_method
  valid_grid_methods = c("linear", "quadratic", "adaptive")
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

  # Validate threshold_method
  valid_threshold_methods = c("fixed", "bootstrap")
  if (!profile_options$threshold_method %in% valid_threshold_methods) {
    stop(sprintf("Error in validateProfileOptions: profile_options$threshold_method must be one of: %s. Got: '%s'",
                 paste(valid_threshold_methods, collapse = ", "), profile_options$threshold_method))
  }

  # Validate bootstrap_conf_level
  if (!is.numeric(profile_options$bootstrap_conf_level) ||
      length(profile_options$bootstrap_conf_level) != 1) {
    stop("Error in validateProfileOptions: profile_options$bootstrap_conf_level must be a single numeric value")
  }
  if (profile_options$bootstrap_conf_level <= 0 || profile_options$bootstrap_conf_level >= 1) {
    stop("Error in validateProfileOptions: profile_options$bootstrap_conf_level must be between 0 and 1")
  }

  # Validate bootstrap parameters if bootstrap method is used
  if (profile_options$threshold_method == "bootstrap") {
    if (is.null(generateData)) {
      stop("Error in computeLikelihoodProfiles: generateData function is required when threshold_method is 'bootstrap'")
    }

    # Check if likelihood function accepts a dataset parameter
    likelihood_args = names(formals(negLogLikelihood))
    if (length(likelihood_args) < 2) {
      stop("Error in computeLikelihoodProfiles: When threshold_method is 'bootstrap', negLogLikelihood function must accept at least 2 arguments (parameters and dataset)")
    }
    if (!("dataset" %in% likelihood_args)) {
      stop("Error in computeLikelihoodProfiles: When threshold_method is 'bootstrap', negLogLikelihood function must accept a 'dataset' argument")
    }
    if (!is.numeric(profile_options$n_bootstrap) || length(profile_options$n_bootstrap) != 1) {
      stop("Error in validateProfileOptions: profile_options$n_bootstrap must be a single numeric value")
    }
    if (profile_options$n_bootstrap <= 0 || profile_options$n_bootstrap != round(profile_options$n_bootstrap)) {
      stop("Error in validateProfileOptions: profile_options$n_bootstrap must be a positive integer")
    }
  }

  # Validate adaptive grid parameters
  if (profile_options$grid_method == "adaptive") {
    # Validate adaptive_max_iterations
    if (!is.numeric(profile_options$adaptive_max_iterations) ||
        length(profile_options$adaptive_max_iterations) != 1) {
      stop("Error in validateProfileOptions: profile_options$adaptive_max_iterations must be a single numeric value")
    }
    if (profile_options$adaptive_max_iterations < 1 ||
        profile_options$adaptive_max_iterations != round(profile_options$adaptive_max_iterations)) {
      stop("Error in validateProfileOptions: profile_options$adaptive_max_iterations must be a positive integer")
    }

    # Validate adaptive_initial_fraction
    if (!is.numeric(profile_options$adaptive_initial_fraction) ||
        length(profile_options$adaptive_initial_fraction) != 1) {
      stop("Error in validateProfileOptions: profile_options$adaptive_initial_fraction must be a single numeric value")
    }
    if (profile_options$adaptive_initial_fraction <= 0 || profile_options$adaptive_initial_fraction > 1) {
      stop("Error in validateProfileOptions: profile_options$adaptive_initial_fraction must be between 0 and 1")
    }

    # Validate adaptive_refinement_fraction
    if (!is.numeric(profile_options$adaptive_refinement_fraction) ||
        length(profile_options$adaptive_refinement_fraction) != 1) {
      stop("Error in validateProfileOptions: profile_options$adaptive_refinement_fraction must be a single numeric value")
    }
    if (profile_options$adaptive_refinement_fraction <= 0 || profile_options$adaptive_refinement_fraction > 1) {
      stop("Error in validateProfileOptions: profile_options$adaptive_refinement_fraction must be between 0 and 1")
    }
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
#' @param profile_options Profiling options
#' @param profile_costs Numeric vector to store profile costs (modified in place)
#' @param lrt_threshold_vec Numeric vector to store likelihood ratio thresholds (modified in place)
#' @param optimizer_exit_flags Character vector to store optimizer exit flags (modified in place)
#' @param optimal_params_matrix Matrix to store optimal parameters for each grid point (modified in place)
#' @param optimization_success Logical vector to track optimization success (modified in place)
#' @param likelihood_args List of arguments to pass to negLogLikelihood function
#' @param nested_cluster Optional cluster object for nested parallel computation
#' @param generateData Function to generate bootstrap datasets
#' @param param_bootstrap_stats Array to store bootstrap statistics for this parameter
#' @param verbose Logical indicating whether to print progress messages
#' @return List with updated arrays: profile_costs, optimizer_exit_flags, optimal_params_matrix, optimization_success, and success_count
#' @keywords internal
optimizeGridDirection = function(grid_indices, param_index, grid_values, warm_start_params,
                                 negLogLikelihood, bounds, optimizer_info, likelihood_args,
                                 profile_options, profile_costs, lrt_threshold_vec,
                                 optimizer_exit_flags, optimal_params_matrix,
                                 optimization_success, nested_cluster = NULL,
                                 generateData = NULL, param_bootstrap_stats = NULL,
                                 verbose = FALSE, adaptive_config = NULL) {

  # Initialize variables for repeat loop
  iteration = 1
  current_warm_start = warm_start_params
  success_count = 0

  # For adaptive grids, track total points used
  if (!is.null(adaptive_config)) {
    points_used = 0
    current_grid = c()
    max_iterations = adaptive_config$max_iterations
    total_budget = adaptive_config$total_grid_points
  } else {
    max_iterations = 1  # Static grids execute once
    current_grid = grid_values
    total_budget = length(grid_values)
  }
  # Initialize bootstrap stats
  if (profile_options$threshold_method == "bootstrap" && is.null(param_bootstrap_stats)) {
    param_bootstrap_stats = array(NA, dim = c(profile_options$n_bootstrap, total_budget))
  }

  repeat {
    # Determine grid points for this iteration
    if (!is.null(adaptive_config)) {
      if (iteration == 1) {
        # First iteration: create initial coarse grid
        n_points_this_iter = round(adaptive_config$initial_fraction * total_budget)
        new_grid_points = createAdaptiveInitialGrid(
          adaptive_config$grid_lower,
          adaptive_config$grid_upper,
          adaptive_config$current_value,
          n_points_this_iter
        )

        if (verbose) {
          cat(sprintf("    Iteration %d: Creating initial grid with %d points\n", iteration, n_points_this_iter))
        }
      } else {
        # Subsequent iterations: refine based on current profile
        remaining_budget = total_budget - points_used
        n_points_this_iter = min(remaining_budget,
                                 round(adaptive_config$refinement_fraction * total_budget))

        if (n_points_this_iter <= 0) break

        new_grid_points = createRefinementPoints(
          current_grid,
          profile_costs,
          lrt_threshold_vec,
          adaptive_config,
          n_points_this_iter
        )

        if (verbose) {
          cat(sprintf("    Iteration %d: Adding %d refinement points\n", iteration, length(new_grid_points)))
        }
      }

      # Merge and sort grid with synchronized array reordering
      merged_grid = c(current_grid, new_grid_points)
      current_grid_size = length(merged_grid)

      # Extend arrays to accommodate new grid size before reordering
      if (current_grid_size > length(profile_costs)) {
        old_size = length(profile_costs)
        extension_size = current_grid_size - old_size

        # Extend all arrays
        profile_costs = c(profile_costs, rep(NA, extension_size))
        lrt_threshold_vec = c(lrt_threshold_vec, rep(profile_options$ll_ratio_threshold, extension_size))
        optimizer_exit_flags = c(optimizer_exit_flags, rep(NA_character_, extension_size))
        optimization_success = c(optimization_success, rep(FALSE, extension_size))

        # Extend matrix
        if (nrow(optimal_params_matrix) < current_grid_size) {
          additional_rows = current_grid_size - nrow(optimal_params_matrix)
          optimal_params_matrix = rbind(optimal_params_matrix,
                                       matrix(NA, nrow = additional_rows, ncol = ncol(optimal_params_matrix)))
        }

        # Extend bootstrap stats if needed
        if (!is.null(param_bootstrap_stats) && ncol(param_bootstrap_stats) < current_grid_size) {
          additional_cols = current_grid_size - ncol(param_bootstrap_stats)
          param_bootstrap_stats = cbind(param_bootstrap_stats,
                                       array(NA, dim = c(nrow(param_bootstrap_stats), additional_cols)))
        }
      }

      # Apply synchronized ordering to grid and ALL arrays
      grid_order = order(merged_grid)
      current_grid = merged_grid[grid_order]
      profile_costs = profile_costs[grid_order]
      lrt_threshold_vec = lrt_threshold_vec[grid_order]
      optimizer_exit_flags = optimizer_exit_flags[grid_order]
      optimization_success = optimization_success[grid_order]
      optimal_params_matrix = optimal_params_matrix[grid_order, , drop = FALSE]
      if (!is.null(param_bootstrap_stats)) {
        param_bootstrap_stats = param_bootstrap_stats[, grid_order, drop = FALSE]
      }

      # Track which grid positions need optimization (only new points)
      # Mark positions that correspond to new points as needing optimization
      old_grid_length = length(current_grid) - length(new_grid_points)

      # After sorting, identify which indices correspond to new points by checking for NA values
      # (new points will have NA in profile_costs since they haven't been optimized yet)
      grid_indices_this_iter = which(is.na(profile_costs))
      points_used = points_used + length(new_grid_points)

    } else {
      # Static grid: use provided grid_indices
      grid_indices_this_iter = grid_indices
    }

    # Execute optimization for current iteration's grid points
    for (j in grid_indices_this_iter) {
      current_param_value = current_grid[j]

      tryCatch({
        result = optimizeConditional(
          param_index = param_index,
          fixed_value = current_param_value,
          warm_start_params = current_warm_start,
          negLogLikelihood = negLogLikelihood,
          bounds = bounds,
          optimizer_info = optimizer_info,
          likelihood_args = likelihood_args,
          nested_cluster = nested_cluster
        )

        profile_costs[j] = result$cost
        optimizer_exit_flags[j] = result$exit_flag
        optimal_params_matrix[j, ] = result$optimal_params

        # Track success for exit flags
        is_success = (result$exit_flag == "success" || result$exit_flag == "direct_evaluation")
        optimization_success[j] = is_success
        if (is_success) success_count = success_count + 1

        # Update likelihood ratio threshold for this grid point
        resultBootstrapLRT = getBootstrapLRTThreshold(
          optimal_params = result$optimal_params,
          param_index = param_index,
          is_success = is_success,
          fixed_value = current_param_value,
          profile_options = profile_options,
          negLogLikelihood = negLogLikelihood,
          bounds = bounds,
          optimizer_info = optimizer_info,
          likelihood_args = likelihood_args,
          nested_cluster = nested_cluster,
          generateData = generateData,
          verbose = verbose
        )
        lrt_threshold_vec[j] = resultBootstrapLRT$threshold
        if (profile_options$threshold_method == "bootstrap" && is_success) {
          param_bootstrap_stats[, j] = resultBootstrapLRT$lrt_stats
        }

        # Update warm start for next iteration (warm start strategy)
        if (is_success) {
          current_warm_start = result$optimal_params
        }
        # On failure, keep using previous warm start parameters

      }, error = function(e) {
        stop(sprintf("Error in optimizeGridDirection: Optimization failed for parameter %d at grid point %d: %s",
                     param_index, j, e$message))
      })
    }

    # Check termination conditions
    if (is.null(adaptive_config) ||
        iteration >= max_iterations ||
        points_used >= total_budget) {
      break
    }

    iteration = iteration + 1
  }

  # Final cleanup and trimming for adaptive grids
  if (!is.null(adaptive_config)) {
    final_grid_size = length(current_grid)
    profile_costs = profile_costs[1:final_grid_size]
    lrt_threshold_vec = lrt_threshold_vec[1:final_grid_size]
    optimizer_exit_flags = optimizer_exit_flags[1:final_grid_size]
    optimal_params_matrix = optimal_params_matrix[1:final_grid_size, , drop = FALSE]
    optimization_success = optimization_success[1:final_grid_size]
    if (!is.null(param_bootstrap_stats)) {
      param_bootstrap_stats = param_bootstrap_stats[, 1:final_grid_size, drop = FALSE]
    }
  }

  # Return results with final_grid for adaptive case
  result_list = list(
    profile_costs = profile_costs,
    lrt_threshold_vec = lrt_threshold_vec,
    optimizer_exit_flags = optimizer_exit_flags,
    optimal_params_matrix = optimal_params_matrix,
    optimization_success = optimization_success,
    success_count = success_count,
    bootstrap_stats = param_bootstrap_stats
  )

  if (!is.null(adaptive_config)) {
    result_list$final_grid = current_grid
  }

  return(result_list)
}

#' Compute Bootstrap Likelihood Ratio Test Threshold
#'
#' Computes empirical likelihood ratio test thresholds using bootstrap resampling.
#' Generates bootstrap datasets from null hypothesis parameters and performs
#' H0 vs H1 optimizations to compute LRT statistics for threshold estimation.
#'
#' @param optimal_params Optimal parameter values at current grid point
#' @param param_index Index of parameter being profiled
#' @param is_success Logical indicating if optimization was successful
#' @param fixed_value Fixed value for the profiled parameter
#' @param profile_options List of profiling options including bootstrap settings
#' @param negLogLikelihood Negative log-likelihood function
#' @param bounds Parameter bounds
#' @param optimizer_info Optimizer configuration
#' @param likelihood_args Additional arguments for likelihood function
#' @param nested_cluster Optional cluster for nested parallel computation
#' @param generateData Function to generate bootstrap datasets from parameters
#' @param verbose Logical for verbose output
#' @return List with threshold value and bootstrap LRT statistics
#' @keywords internal
getBootstrapLRTThreshold = function(optimal_params, param_index, is_success, fixed_value,
                                    profile_options, negLogLikelihood, bounds,
                                    optimizer_info, likelihood_args, nested_cluster = NULL,
                                    generateData = NULL, verbose = FALSE) {

  # If not bootstrap method or optimization failed, use fixed threshold
  if (profile_options$threshold_method != "bootstrap" || !is_success) {
    return(list(threshold = profile_options$ll_ratio_threshold,
                lrt_stats = NULL))
  }

  # Generate bootstrap samples from null hypothesis (parameters at this grid point)
  n_boot = profile_options$n_bootstrap

  # Storage for LRT statistics
  lrt_stats = rep(NA, n_boot)

  # Generate bootstrap datasets and compute LRT statistics
  for (b in seq_len(n_boot)) {
    tryCatch({
      # Generate dataset from parameters at the current grid point (H0)
      # This ensures we're sampling from the null hypothesis distribution
      bootstrap_data = generateData(optimal_params)

      likelihood_args_loc = utils::modifyList(likelihood_args, list(dataset = bootstrap_data))
      # Create bootstrap-specific likelihood function
      bootstrap_likelihood = function(params) {
        # Call original likelihood with bootstrap data
        do.call(negLogLikelihood, c(list(params), likelihood_args_loc))
      }

      # 1. Optimize without constraints (H1)
      h1_result = optimizeFull(
        start_params = optimal_params,
        negLogLikelihood = bootstrap_likelihood,
        bounds = bounds,
        optimizer_info = optimizer_info,
        nested_cluster = nested_cluster
      )

      # 2. Optimize with parameter fixed (H0)
      h0_result = optimizeConditional(
        param_index = param_index,
        fixed_value = fixed_value,
        warm_start_params = optimal_params,
        negLogLikelihood = bootstrap_likelihood,
        bounds = bounds,
        optimizer_info = optimizer_info,
        nested_cluster = nested_cluster
      )

      # Calculate LRT statistic: 2 * (L_H0 - L_H1)
      lrt_stats[b] = 2 * (h0_result$cost - h1_result$cost)

    }, error = function(e) {
      stop(sprintf("Error in getBootstrapLRTThreshold: Bootstrap iteration %d failed: %s",
                   b, e$message))
    })
  }

  if (sum(!is.na(lrt_stats)) < 100 && verbose && n_boot >= 200){
    warning(paste0("Warning in getBootstrapLRTThreshold: Less than 100 valid bootstrap iterations for parameter ",
                   param_index, ". This may lead to unreliable threshold estimation."))
  }

  # Calculate empirical threshold as the quantile of bootstrap LRT statistics
  threshold = stats::quantile(lrt_stats, profile_options$bootstrap_conf_level, na.rm = TRUE)

  list(threshold = threshold,
       lrt_stats = lrt_stats)
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
#' @param nested_cluster_call Expression to create nested cluster for DEoptim
#' @param verbose Logical for verbose output
#' @param generateData Function to generate bootstrap datasets
#' @return List with profile data and confidence interval
#' @keywords internal
computeParameterProfile = function(param_index, params_current, negLogLikelihood,
                                   bounds, profile_options, optimizer_info,
                                   likelihood_args = list(),
                                   nested_cluster_call = NULL,
                                   verbose = FALSE,
                                   generateData = NULL) {

  if (!is.null(nested_cluster_call) && optimizer_info$name == "deoptim") {
    # Evaluate the expression 'nested_cluster_call'
    # It defines the object `nested_cluster`
    tryCatch({
      eval(nested_cluster_call)
      if (!exists("nested_cluster") || !inherits(nested_cluster, "cluster")) {
        stop("Error in computeParameterProfile: nested_cluster_call did not create 'nested_cluster' object, or it is not a cluster.")
      }
    }, error = function(e) {
      stop(sprintf("Error in computeParameterProfile: Failed to evaluate nested_cluster_call: %s", e$message))
    })
  } else{
    nested_cluster = NULL
  }

  # Create grid for this parameter
  grid_result = createParameterGrid(
    param_index = param_index,
    params_current = params_current,
    bounds = bounds,
    profile_options = profile_options
  )

  # Check if this is an adaptive grid configuration
  if (!is.null(grid_result$grid_method) && grid_result$grid_method == "adaptive") {
    # For adaptive grids, we'll pass the configuration to optimizeGridDirection
    current_grid = c()  # Will be built iteratively
    adaptive_config = grid_result
    n_grid = 0  # Will grow during iterations
    max_size = adaptive_config$total_grid_points
  } else {
    # Standard static grid
    current_grid = grid_result$grid_values
    adaptive_config = NULL
    n_grid = length(current_grid)
    max_size = n_grid
  }

  # Initialize exit flags structure
  exit_flags = list(
    low_success_rate = FALSE,
    confidence_interval_failed = FALSE,
    optimal_value_out_of_bounds = FALSE,
    grid_issue = FALSE,
    grid_too_narrow = FALSE,
    few_within_threshold = FALSE
  )
  if (!is.null(grid_result$optimal_value_out_of_bounds)) {
    exit_flags$optimal_value_out_of_bounds = grid_result$optimal_value_out_of_bounds
  }
  if (!is.null(grid_result$grid_issue)) {
    exit_flags$grid_issue = grid_result$grid_issue
  }

  # Initialize arrays
  profile_costs = rep(NA, max_size)
  lrt_threshold_vec = rep(profile_options$ll_ratio_threshold, max_size)
  optimizer_exit_flags = character(max_size)
  optimal_params_matrix = matrix(NA, nrow = max_size, ncol = length(params_current))
  optimization_success = logical(max_size)

  if (verbose) {
    cat(sprintf("  Evaluating %d grid points...\n", max_size))
  }

  # Find the index of the optimal parameter value in the grid
  if (!is.null(adaptive_config)) {
    optimal_index = 1  # Will be updated after first iteration
  } else {
    optimal_index = which.min(abs(current_grid - params_current[param_index]))
  }

  # Initialize warm start parameters
  warm_start_params = params_current

  # Initialize success tracking
  success_count = 0

  # For adaptive grids, use different optimization approach
  if (!is.null(adaptive_config)) {
    # Adaptive grid optimization with iterative refinement
    adaptive_result = optimizeGridDirection(
      grid_indices = NULL,  # Will be determined iteratively
      param_index = param_index,
      grid_values = current_grid,  # Use current_grid consistently
      warm_start_params = params_current,
      negLogLikelihood = negLogLikelihood,
      bounds = bounds,
      optimizer_info = optimizer_info,
      likelihood_args = likelihood_args,
      profile_options = profile_options,
      profile_costs = profile_costs,
      lrt_threshold_vec = lrt_threshold_vec,
      optimizer_exit_flags = optimizer_exit_flags,
      optimal_params_matrix = optimal_params_matrix,
      optimization_success = optimization_success,
      nested_cluster = nested_cluster,
      generateData = generateData,
      param_bootstrap_stats = NULL,
      verbose = verbose,
      adaptive_config = adaptive_config
    )

    # Extract final results from adaptive optimization
    current_grid = adaptive_result$final_grid
    profile_costs = adaptive_result$profile_costs
    lrt_threshold_vec = adaptive_result$lrt_threshold_vec
    optimizer_exit_flags = adaptive_result$optimizer_exit_flags
    optimal_params_matrix = adaptive_result$optimal_params_matrix
    optimization_success = adaptive_result$optimization_success
    success_count = adaptive_result$success_count
    param_bootstrap_stats = adaptive_result$bootstrap_stats
    n_grid = length(current_grid)

    # Update optimal index for adaptive case
    optimal_index = which.min(abs(current_grid - params_current[param_index]))

  } else {
    # Static grid optimization - helper function for shared optimization arguments
    get_optimization_args = function(grid_indices, warm_start, bootstrap_stats = NULL) {
      list(
        grid_indices = grid_indices,
        param_index = param_index,
        grid_values = current_grid,
        warm_start_params = warm_start,
        negLogLikelihood = negLogLikelihood,
        bounds = bounds,
        optimizer_info = optimizer_info,
        likelihood_args = likelihood_args,
        profile_options = profile_options,
        profile_costs = profile_costs,
        lrt_threshold_vec = lrt_threshold_vec,
        optimizer_exit_flags = optimizer_exit_flags,
        optimal_params_matrix = optimal_params_matrix,
        optimization_success = optimization_success,
        nested_cluster = nested_cluster,
        generateData = generateData,
        param_bootstrap_stats = bootstrap_stats,
        verbose = verbose
      )
    }

    # Optimize right side (from optimal outward)
    right_result = do.call(optimizeGridDirection, get_optimization_args(optimal_index:n_grid, params_current))

    # Update arrays and collect bootstrap stats
    profile_costs = right_result$profile_costs
    lrt_threshold_vec = right_result$lrt_threshold_vec
    optimizer_exit_flags = right_result$optimizer_exit_flags
    optimal_params_matrix = right_result$optimal_params_matrix
    optimization_success = right_result$optimization_success
    success_count = success_count + right_result$success_count
    param_bootstrap_stats = right_result$bootstrap_stats

    # Optimize left side if needed
    if (optimal_index > 1) {
      left_result = do.call(optimizeGridDirection, get_optimization_args(seq(optimal_index-1, 1), params_current, param_bootstrap_stats))

      # Update arrays with left side results
      profile_costs = left_result$profile_costs
      lrt_threshold_vec = left_result$lrt_threshold_vec
      optimizer_exit_flags = left_result$optimizer_exit_flags
      optimal_params_matrix = left_result$optimal_params_matrix
      optimization_success = left_result$optimization_success
      success_count = success_count + left_result$success_count
      param_bootstrap_stats = left_result$bootstrap_stats
    }
  } # End of static grid optimization

  # Calculate success rate and set exit flags
  success_rate = success_count / n_grid
  if (success_rate < 0.5) {
    exit_flags$low_success_rate = TRUE
  }

  # Check if corresponding cost is below threshold. If not set current params at this index
  # to the minimum observed cost with a warning.
  optimal_index = which.min(abs(current_grid - params_current[param_index]))  # Use current_grid consistently
  min_cost_index = which.min(profile_costs)
  min_observed_cost = profile_costs[min_cost_index]
  optimal_grid_cost = profile_costs[optimal_index]

  # If the cost at optimal parameter position is significantly higher than minimum found
  ll_ratio_optimal = 2 * (optimal_grid_cost - min_observed_cost)
  lrt_threshold_optimal = lrt_threshold_vec[optimal_index]
  if (!is.na(optimal_grid_cost) && ll_ratio_optimal >= lrt_threshold_optimal) {
    warning(sprintf("Warning in computeLikelihoodProfiles: For parameter %d, likelihood ratio test at provided optimal parameter (%.4f) exceeds threshold (%.4f). Using best parameters found during profiling.",
                    param_index, ll_ratio_optimal, lrt_threshold_optimal), call. = FALSE)

    # Update the profile data to use the best parameters found
    optimal_param_value = current_grid[min_cost_index]
  } else {
    optimal_param_value = params_current[param_index]
  }

  # Create profile data structure
  profile_data = list(
    param_index = param_index,
    optimal_param_value = optimal_param_value,
    grid_values = current_grid,  # Always use current_grid as authoritative
    profile_costs = profile_costs,
    lrt_threshold_vec = lrt_threshold_vec,
    optimizer_exit_flags = optimizer_exit_flags,
    optimal_params_matrix = optimal_params_matrix,
    optimization_success = optimization_success,
    success_rate = success_rate,
    exit_flags = exit_flags,
    optimal_cost = min_observed_cost,
    failed = FALSE
  )

  # Extract confidence interval
  ci_result = extractConfidenceInterval(
    profile_data = profile_data,
    profile_options = profile_options
  )

  confidence_interval = ci_result$confidence_bounds

  # Update exit flags based on CI extraction
  if (any(is.na(confidence_interval))) {
    profile_data$exit_flags$confidence_interval_failed = TRUE
  }

  # Update exit flags for grid too narrow condition
  if (ci_result$grid_too_narrow) {
    profile_data$exit_flags$grid_too_narrow = TRUE
  }

  if (ci_result$few_within_threshold) {
    profile_data$exit_flags$few_within_threshold = TRUE
  }

  return(list(
    profile = profile_data,
    confidence_interval = confidence_interval,
    bootstrap_stats = param_bootstrap_stats
  ))

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

  if (profile_options$grid_method == "linear") {
    grid_values = createLinearGrid(grid_lower, grid_upper, current_value, profile_options$grid_points)
  } else if (profile_options$grid_method == "quadratic") {
    grid_values = createQuadraticGrid(grid_lower, grid_upper, current_value, profile_options$grid_points)
  } else if (profile_options$grid_method == "adaptive") {
    # Return adaptive configuration instead of grid values
    return(list(
      grid_method = "adaptive",
      grid_lower = grid_lower,
      grid_upper = grid_upper,
      current_value = current_value,
      total_grid_points = profile_options$grid_points,
      max_iterations = profile_options$adaptive_max_iterations,
      initial_fraction = profile_options$adaptive_initial_fraction,
      refinement_fraction = profile_options$adaptive_refinement_fraction,
      optimal_value_out_of_bounds = optimal_value_out_of_bounds,
      grid_issue = FALSE
    ))
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

#' Create Linear Grid
#'
#' Creates a grid with linear spacing on each side of the current value (optimum).
#' Left side: linear spacing from `grid_lower` to `current_value`.
#' Right side: linear spacing from `current_value` to `grid_upper`.
#' Note: spacing may differ between left and right sides when `current_value` is off-center.
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param current_value Center value for grid construction
#' @param n_points Number of grid points
#' @return Numeric vector with linear spacing on each side of `current_value`
#' @keywords internal
createLinearGrid = function(grid_lower, grid_upper, current_value, n_points) {
  res_grid = c(seq(grid_lower, current_value, length.out = floor(n_points/2) + 1),
               seq(current_value, grid_upper, length.out = floor(n_points/2) + 1))

  unique(res_grid)
}

#' Create Quadratic Grid
#'
#' Creates a grid with quadratic spacing that is denser near the current value (optimum).
#' Uses quadratic transformation t^2 to create progressively wider spacing
#' as distance from current value increases. Provides symmetric behavior regardless
#' of current value position within bounds.
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param center_value Central value for denser sampling
#' @param n_points Number of grid points
#' @return Numeric vector with quadratic spacing (dense near `center_value`)
#' @keywords internal
createQuadraticGrid = function(grid_lower, grid_upper, center_value, n_points) {

  # Create quadratic spacing for denser sampling near center
  n_half = floor(n_points/2) + 1
  t_onesided = seq(0, 1, length.out = n_half)^2
  left_side = center_value - t_onesided*(center_value - grid_lower)
  right_side = center_value + t_onesided*(grid_upper - center_value)

  res_grid = c(left_side, right_side)
  unique(sort(res_grid))
}

#' Create Adaptive Initial Grid
#'
#' Creates the initial coarse grid for adaptive grid method using quadratic spacing.
#' This provides a good starting distribution focused around the current parameter value.
#'
#' @param grid_lower Lower bound for grid
#' @param grid_upper Upper bound for grid
#' @param center_value Central value for denser sampling
#' @param n_points Number of grid points
#' @return Numeric vector with quadratic spacing for initial adaptive grid
#' @keywords internal
createAdaptiveInitialGrid = function(grid_lower, grid_upper, center_value, n_points) {
  createQuadraticGrid(grid_lower, grid_upper, center_value, n_points)
}

#' Identify Refinement Regions for Adaptive Grid
#'
#' Analyzes current profile likelihood results to identify regions where additional
#' grid points would be most beneficial, typically near confidence interval boundaries
#' and areas with sparse coverage relative to likelihood curvature.
#'
#' @param current_grid Current grid parameter values
#' @param profile_costs Current profile likelihood costs
#' @param lrt_thresholds Likelihood ratio thresholds for confidence bounds
#' @param adaptive_config Adaptive grid configuration parameters
#' @return List of refinement regions with importance weights and suggested point counts
#' @keywords internal
identifyRefinementRegions = function(current_grid, profile_costs, lrt_thresholds, adaptive_config) {

  if (length(current_grid) < 3 || all(is.na(profile_costs))) {
    return(list(regions = list(), total_importance = 0))
  }

  # Sort by grid values for analysis
  sorted_indices = order(current_grid)
  sorted_grid = current_grid[sorted_indices]
  sorted_costs = profile_costs[sorted_indices]

  # Convert costs to likelihood ratios
  min_cost = min(sorted_costs, na.rm = TRUE)
  ll_ratios = 2 * (sorted_costs - min_cost)

  # Find approximate CI boundaries by looking for threshold crossings
  refinement_regions = list()
  region_count = 0

  # Method 1: Focus on threshold crossing regions using helper function
  threshold_vec = lrt_thresholds[sorted_indices]
  crossing_points = findCrossings(sorted_grid, ll_ratios, threshold_vec)

  # Convert crossings to refinement regions
  for (crossing in crossing_points) {
    # Find the interval that contains this crossing
    for (i in seq_len(length(sorted_grid) - 1)) {
      if (sorted_grid[i] <= crossing && crossing <= sorted_grid[i+1]) {
        region_count = region_count + 1
        refinement_regions[[region_count]] = list(
          center = crossing,  # Use the precise crossing point
          width = abs(sorted_grid[i+1] - sorted_grid[i]),
          importance = 1.0,  # High importance for CI boundaries
          reason = "threshold_crossing"
        )
        break
      }
    }
  }

  # Method 2: Focus on large gaps between points (sparse coverage)
  grid_gaps = diff(sorted_grid)
  if (length(grid_gaps) > 0) {
    large_gap_threshold = stats::quantile(grid_gaps, 0.75, na.rm = TRUE)

    for (i in seq_len(length(grid_gaps))) {
      if (grid_gaps[i] > large_gap_threshold) {
        region_count = region_count + 1
        refinement_regions[[region_count]] = list(
          center = (sorted_grid[i] + sorted_grid[i+1]) / 2,
          width = grid_gaps[i],
          importance = 0.1,  # Less importance for gap filling
          reason = "large_gap"
        )
      }
    }
  }

  list(regions = refinement_regions, total_importance = sum(sapply(refinement_regions, function(r) r$importance)))
}

#' Create Refinement Points for Adaptive Grid
#'
#' Generates additional grid points focused on the identified refinement regions.
#' Distributes points based on region importance and avoids duplicating existing points.
#'
#' @param current_grid Existing grid parameter values
#' @param profile_costs Current profile likelihood costs
#' @param lrt_thresholds Likelihood ratio thresholds
#' @param adaptive_config Adaptive grid configuration
#' @param n_new_points Number of new grid points to create
#' @return Numeric vector of new grid points for refinement
#' @keywords internal
createRefinementPoints = function(current_grid, profile_costs, lrt_thresholds,
                                  adaptive_config, n_new_points) {

  if (n_new_points <= 0) return(numeric(0))

  # Identify regions needing refinement
  refinement_analysis = identifyRefinementRegions(current_grid, profile_costs, lrt_thresholds, adaptive_config)

  if (length(refinement_analysis$regions) == 0) {
    # No specific regions identified, add points in largest gaps
    sorted_grid = sort(current_grid)
    if (length(sorted_grid) < 2) return(numeric(0))

    gaps = diff(sorted_grid)
    largest_gap_idx = which.max(gaps)
    gap_centers = (sorted_grid[largest_gap_idx] + sorted_grid[largest_gap_idx + 1]) / 2

    # Create points around the largest gap
    gap_width = gaps[largest_gap_idx]
    new_points = seq(gap_centers - gap_width/4, gap_centers + gap_width/4, length.out = n_new_points)
    return(new_points)
  }

  # Distribute new points among identified regions based on importance
  total_importance = refinement_analysis$total_importance
  new_points = numeric(0)
  points_allocated = 0

  for (region in refinement_analysis$regions) {
    if (points_allocated >= n_new_points) break

    # Allocate points based on region importance
    points_for_region = max(1, round((region$importance / total_importance) * n_new_points))
    points_for_region = min(points_for_region, n_new_points - points_allocated)

    if (points_for_region > 0) {
      # Create points in this region, avoiding existing points
      region_lower = region$center - region$width/2
      region_upper = region$center + region$width/2

      # Use quadratic spacing within the region for better distribution
      region_points = createQuadraticGrid(region_lower, region_upper, region$center, points_for_region)

      new_points = c(new_points, region_points)
      points_allocated = points_allocated + points_for_region
      # Removed min distance as it defeats the purpose as is. Might consider using alternative.
      # Filter out points too close to existing ones
      # min_distance = min(diff(sort(current_grid)), na.rm = TRUE) / 4  # Minimum separation
      # for (point in region_points) {
      #   if (all(abs(current_grid - point) > min_distance) && all(abs(new_points - point) > min_distance)) {
      #     new_points = c(new_points, point)
      #     points_allocated = points_allocated + 1
      #     if (points_allocated >= n_new_points) break
      #   }
      # }
    }
  }

  # If we still need more points, fill in remaining budget with gap-filling
  if (points_allocated < n_new_points && length(current_grid) > 1) {
    remaining = n_new_points - points_allocated
    all_points = sort(c(current_grid, new_points))
    gaps = diff(all_points)

    for (i in seq_len(min(remaining, length(gaps)))) {
      largest_gap_idx = which.max(gaps)
      if (!is.na(gaps[largest_gap_idx]) && gaps[largest_gap_idx] > 0) {
        gap_center = (all_points[largest_gap_idx] + all_points[largest_gap_idx + 1]) / 2
        new_points = c(new_points, gap_center)

        # Update gaps array to reflect the new point
        all_points = sort(c(all_points, gap_center))
        gaps = diff(all_points)
      }
    }
  }

  unique(new_points)
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
#' @param nested_cluster Optional cluster object for nested parallel computation
#' @return List with optimized cost, exit flag, and optimal parameters
#' @keywords internal
optimizeConditional = function(param_index, fixed_value, warm_start_params,
                               negLogLikelihood, bounds, optimizer_info, likelihood_args = list(),
                               nested_cluster = NULL) {

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

      optimizer_info$extra_args$initialpop = getInitPopOptim(optimizer_info$extra_args$initialpop,
                                                             n_pop = optimizer_info$extra_args$NP,
                                                             free_indices = free_indices,
                                                             lower_free = lower_free,
                                                             upper_free = upper_free,
                                                             initial_free_params = initial_free_params)
      if (!is.null(nested_cluster)){
        optimizer_info$extra_args$cluster = nested_cluster  # Pass nested cluster if available
      }
      deoptim_control = do.call(DEoptim::DEoptim.control, optimizer_info$extra_args)

      # Suppress the specific NP warning while preserving other warnings
      result = withCallingHandlers({
        DEoptim::DEoptim(fn = conditional_cost,
                         lower = lower_free,
                         upper = upper_free,
                         control = deoptim_control)
      }, warning = function(w) {
        # Suppress warnings about NP population size recommendation
        if (grepl("For many problems it is best to set 'NP'", w$message, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })

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

  return(list(cost = cost,
              exit_flag = exit_flag,
              optimal_params = optimal_params))
}

#' Evaluate Conditional Cost Function
#'
#' Evaluates the cost function with one parameter fixed at a specified value.
#' Reconstructs the full parameter vector and calls the original likelihood function.
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
  res = do.call(negLogLikelihood, c(list(full_params), likelihood_args))
  res
}

#' Find Threshold Crossings in Profile Data
#'
#' Helper function that finds where likelihood ratios cross thresholds.
#' Uses interpolation for accurate crossing point estimation.
#'
#' @param sorted_grid Sorted grid parameter values
#' @param ll_ratios Likelihood ratio values (same order as sorted_grid)
#' @param lrt_threshold_vec Threshold values for each grid point (same order as sorted_grid)
#' @return Numeric vector of crossing points
#' @keywords internal
findCrossings = function(sorted_grid, ll_ratios, lrt_threshold_vec) {

  if (length(sorted_grid) < 2) return(numeric(0))

  n_grid_values = length(sorted_grid)
  crossings = c()

  for (i in seq_len(n_grid_values-1)) {
    # Check if threshold is crossed between points i and i+1
    current_threshold = lrt_threshold_vec[i]
    next_threshold = lrt_threshold_vec[i+1]

    if (!is.na(ll_ratios[i]) && !is.na(ll_ratios[i+1]) &&
        !is.na(current_threshold) && !is.na(next_threshold)) {

      if ((ll_ratios[i] <= current_threshold && ll_ratios[i+1] > next_threshold) ||
          (ll_ratios[i] > current_threshold && ll_ratios[i+1] <= next_threshold)) {

        # For bootstrap case, use simple midpoint interpolation
        # For fixed threshold case, use linear interpolation
        if (current_threshold == next_threshold && abs(ll_ratios[i+1] - ll_ratios[i]) > 1e-10) {
          # Original linear interpolation for fixed threshold
          t = (current_threshold - ll_ratios[i]) / (ll_ratios[i+1] - ll_ratios[i])
          crossing_point = sorted_grid[i] + t * (sorted_grid[i+1] - sorted_grid[i])
        } else {
          # Simple midpoint interpolation for bootstrap case
          crossing_point = (sorted_grid[i] + sorted_grid[i+1]) / 2
        }
        crossings = c(crossings, crossing_point)
      }
    }
  }

  return(crossings)
}

#' Extract Confidence Interval from Profile Data
#'
#' Extracts confidence interval bounds from profile likelihood data by finding
#' where the likelihood ratio crosses the specified threshold.
#'
#' @param profile_data Profile data structure
#' @param profile_options Profiling options
#' @return List with confidence_bounds and flags grid_too_narrow and few_within_threshold
#' @keywords internal
extractConfidenceInterval = function(profile_data, profile_options) {

  if (profile_data$failed || length(profile_data$profile_costs) == 0) {
    return(list(confidence_bounds = c(NA, NA), grid_too_narrow = FALSE, few_within_threshold = FALSE))
  }

  # Find confidence bounds
  result = findConfidenceBounds(
    grid_values = profile_data$grid_values,
    profile_costs = profile_data$profile_costs,
    lrt_threshold_vec = profile_data$lrt_threshold_vec,
    optimal_param_value = profile_data$optimal_param_value
  )

  return(result)
}

#' Find Confidence Bounds Using Threshold Crossing Detection
#'
#' Finds confidence interval bounds by detecting where likelihood ratios cross
#' thresholds, using interpolation for accurate boundary estimation.
#'
#' @param grid_values Grid parameter values
#' @param profile_costs Profile cost values
#' @param lrt_threshold_vec Threshold likelihood ratios for confidence level
#' @param optimal_param_value Optimal parameter value from params_current
#' @return List with confidence_bounds vector and flags: grid_too_narrow, few_within_threshold
#' @keywords internal
findConfidenceBounds = function(grid_values, profile_costs, lrt_threshold_vec, optimal_param_value) {

  few_within_threshold = FALSE
  if (length(grid_values) < 2) {
    return(list(confidence_bounds = c(NA, NA), grid_too_narrow = FALSE, few_within_threshold = few_within_threshold))
  }

  # Sort by grid values
  sorted_indices = order(grid_values)
  sorted_grid = grid_values[sorted_indices]
  sorted_costs = profile_costs[sorted_indices]

  n_grid_values = length(sorted_grid)

  # Convert costs to likelihood ratios
  ll_ratios = 2 * (sorted_costs - min(sorted_costs, na.rm = TRUE))

  # Check how many points are within the threshold
  within_threshold = ll_ratios <= lrt_threshold_vec
  within_threshold_count = sum(within_threshold, na.rm = TRUE)
  if (within_threshold_count < 3) {
    few_within_threshold = TRUE
  }

  # Check if all points are within threshold (grid too narrow condition)
  if (all(within_threshold, na.rm = TRUE)) {
    # All points within threshold - grid too narrow
    return(list(
      confidence_bounds = c(min(sorted_grid), max(sorted_grid)),
      grid_too_narrow = TRUE,
      few_within_threshold = few_within_threshold
    ))
  }

  # Find crossings where likelihood ratio crosses threshold
  crossings = findCrossings(sorted_grid, ll_ratios, lrt_threshold_vec)

  # Handle edge cases
  if (length(crossings) == 0) {
    # No crossings found - check if all points are within threshold
    if (sum(within_threshold, na.rm = TRUE) == 0) {
      return(list(confidence_bounds = c(NA, NA), grid_too_narrow = FALSE, few_within_threshold = few_within_threshold))
    } else {
      # Some points within threshold - use their extremes
      within_indices = which(within_threshold & !is.na(within_threshold))
      return(list(confidence_bounds = c(min(sorted_grid[within_indices]), max(sorted_grid[within_indices])),
                  grid_too_narrow = FALSE,
                  few_within_threshold = few_within_threshold))
    }
  }

  # Find regions within threshold
  if (sum(within_threshold, na.rm = TRUE) == 0) {
    return(list(confidence_bounds = c(NA, NA), grid_too_narrow = FALSE, few_within_threshold = few_within_threshold))
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

  return(list(confidence_bounds = c(lower_bound, upper_bound), grid_too_narrow = FALSE, few_within_threshold = few_within_threshold))
}

#' Issue Consolidated Warnings for Profile Results
#'
#' Analyzes profile results and issues consolidated warnings for common issues
#' such as failed profiles, low success rates, and confidence interval problems.
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
  grid_too_narrow_params = c()
  few_within_threshold_params = c()

  # Analyze exit flags
  for (i in 1:n_params) {
    if (!profiles[[i]]$failed && !is.null(profiles[[i]]$exit_flags)) {
      flags = profiles[[i]]$exit_flags
      if (isTRUE(flags$low_success_rate)) low_success_rate_params = c(low_success_rate_params, i)
      if (isTRUE(flags$confidence_interval_failed)) ci_failed_params = c(ci_failed_params, i)
      if (isTRUE(flags$optimal_value_out_of_bounds)) bounds_issue_params = c(bounds_issue_params, i)
      if (isTRUE(flags$grid_issue)) grid_issue_params = c(grid_issue_params, i)
      if (isTRUE(flags$grid_too_narrow)) grid_too_narrow_params = c(grid_too_narrow_params, i)
      if (isTRUE(flags$few_within_threshold)) few_within_threshold_params = c(few_within_threshold_params, i)
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

  if (length(grid_too_narrow_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Grid too narrow to capture confidence interval bounds for parameter(s): %s. Consider increasing max_grid_range_multiplier.",
                    paste(grid_too_narrow_params, collapse = ", ")), call. = FALSE)
  }

  if (length(few_within_threshold_params) > 0) {
    warning(sprintf("Warning in computeLikelihoodProfiles: Few grid points (<3) within threshold for parameter(s): %s. Confidence intervals may be unreliable. Consider using a finer grid or a different grid method.",
                    paste(few_within_threshold_params, collapse = ", ")), call. = FALSE)
  }
}

#' Validate Profile Computation Results
#'
#' Validates the results of profile likelihood computation, checking for
#' successful profiles and confidence interval quality. Issues warnings
#' for problematic results.
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

#' Create Profile Likelihood Summary Table
#'
#' Creates a summary data frame with parameter names, current values,
#' confidence intervals, and computation status for all profiled parameters.
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

#' Generate Initial Population for DEoptim using Latin Hypercube Sampling
#'
#' Creates an initial population matrix for DEoptim optimization using Latin hypercube
#' sampling to ensure good coverage of the parameter space. Always includes the
#' provided initial parameters as the first member of the population.
#'
#' @param initialpop Optional pre-existing initial population matrix. If provided,
#'   this is returned unchanged after validation
#' @param n_pop Population size (number of individuals in the population)
#' @param free_indices Indices of free parameters being optimized
#' @param lower_free Lower bounds for free parameters
#' @param upper_free Upper bounds for free parameters
#' @param initial_free_params Initial values for free parameters
#' @return Matrix with n_pop rows and length(free_indices) columns containing
#'   the initial population, with initial_free_params as the first row
#' @keywords internal
getInitPopOptim = function(initialpop, n_pop, free_indices, lower_free, upper_free, initial_free_params) {

  # If initialpop is provided, validate and return it
  if (!is.null(initialpop)) {
    if (!is.matrix(initialpop)) {
      stop("Error in getInitPopOptim: initialpop must be a matrix")
    }

    if (ncol(initialpop) != length(free_indices)) {
      stop(sprintf("Error in getInitPopOptim: initialpop must have %d columns (one for each free parameter), got %d",
                   length(free_indices), ncol(initialpop)))
    }

    # Validate bounds
    if (any(t(initialpop) < lower_free) ||
        any(t(initialpop) > upper_free)) {
      stop("Error in getInitPopOptim: All values in initialpop must be within bounds [lower_free, upper_free]")
    }

    initialpop = rbind(initial_free_params, initialpop)
    return(initialpop)
  }

  n_dim = length(free_indices)
  if (is.null(n_pop)){
    n_pop = 10 * n_dim
  }
  # Validate inputs
  if (n_pop < 4) {
    stop("Error in getInitPopOptim: n_pop must be at least 4 if provided")
  }

  # Generate Latin hypercube sample for remaining population members
  n_lhs = n_pop - 1  # Reserve first spot for initial_free_params
  lhs_sample = lhs::randomLHS(n_lhs, n_dim)
  lhs_sample = t(lower_free + (upper_free - lower_free) * t(lhs_sample))

  # Combine initial parameters with LHS sample
  initialpop = rbind(initial_free_params, lhs_sample)
  return(initialpop)
}

#' Local lapply function for conditional parallel execution
#' @param cluster Cluster object or NULL for sequential execution
#' @param X Vector to apply function over
#' @param FUN Function to apply
#' @param ... Additional arguments to FUN
#' @return List of results
#' @keywords internal
locLapply = function(cluster, X, FUN, ...) {
  if (!is.null(cluster)) {
    parallel::parLapply(cluster, X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

#' Optimize without parameter constraints
#'
#' @param start_params Starting parameter values
#' @param negLogLikelihood Cost function to optimize
#' @param bounds Parameter bounds
#' @param optimizer_info Optimizer configuration
#' @param nested_cluster Optional cluster object for nested parallel computation
#' @return List with optimized cost and optimal parameters
#' @keywords internal
optimizeFull = function(start_params, negLogLikelihood, bounds, optimizer_info, nested_cluster = NULL) {

  # Optimize based on optimizer type
  if (optimizer_info$type == "builtin") {
    if (optimizer_info$name == "optim") {
      result = do.call(stats::optim, c(
        list(
          par = start_params,
          fn = negLogLikelihood,
          lower = bounds$lower,
          upper = bounds$upper
        ),
        optimizer_info$extra_args
      ))

      cost = result$value
      optimal_params = result$par

    } else if (optimizer_info$name == "deoptim") {
      # Handle DEoptim optimization (similar to optimizeConditional)
      optimizer_info$extra_args$initialpop = getInitPopOptim(optimizer_info$extra_args$initialpop,
                                                             n_pop = optimizer_info$extra_args$NP,
                                                             free_indices = seq_along(start_params),
                                                             lower_free = bounds$lower,
                                                             upper_free = bounds$upper,
                                                             initial_free_params = start_params)

      if (!is.null(nested_cluster)) {
        optimizer_info$extra_args$cluster = nested_cluster
      }

      deoptim_control = do.call(DEoptim::DEoptim.control, optimizer_info$extra_args)

      result = withCallingHandlers({
        DEoptim::DEoptim(fn = negLogLikelihood,
                         lower = bounds$lower,
                         upper = bounds$upper,
                         control = deoptim_control)
      }, warning = function(w) {
        if (grepl("For many problems it is best to set 'NP'", w$message, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })

      cost = result$optim$bestval
      optimal_params = result$optim$bestmem
    }
  } else if (optimizer_info$type == "function") {
    # Custom optimizer
    result = do.call(optimizer_info$func, c(
      list(
        fn = negLogLikelihood,
        par = start_params,
        lower = bounds$lower,
        upper = bounds$upper
      ),
      optimizer_info$extra_args
    ))

    if (is.list(result)) {
      cost = result$value
      optimal_params = result$par
    } else {
      cost = result
      optimal_params = start_params  # Default if not returned by optimizer
    }
  }

  return(list(cost = cost, optimal_params = optimal_params))
}

