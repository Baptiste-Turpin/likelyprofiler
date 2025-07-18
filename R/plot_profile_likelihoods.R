#' Plot Profile Likelihood Results
#'
#' Create plots of profile likelihood results using base R graphics.
#' Displays profile likelihood curves for all parameters in a grid layout.
#'
#' @param profile_result Result from computeAllParameterProfiles()
#' @param options List of plotting options (see Details)
#' @param true_values Optional vector of true parameter values to display as vertical lines
#'
#' @details
#' The \code{options} list can contain:
#' \itemize{
#'   \item \code{max_cols}: Maximum number of columns in grid layout (default: 3)
#'   \item \code{show_ci}: Show confidence interval bounds as vertical lines (default: TRUE)
#'   \item \code{show_optimal}: Show optimal parameter value as vertical line (default: TRUE)
#'   \item \code{show_threshold}: Show likelihood threshold as horizontal line (default: TRUE)
#'   \item \code{log_scale}: Plot in log-likelihood scale (default: FALSE)
#'   \item \code{main_title}: Overall plot title (default: auto-generated)
#'   \item \code{cex}: Character expansion factor for text (default: 1.0)
#'   \item \code{lwd}: Line width for curves (default: 1.5)
#'   \item \code{mar}: Plot margins (default: c(4, 4, 2, 1))
#' }
#'
#' @return Invisibly returns the original plotting parameters (for restoration)
#'
#' @examples
#' # Basic usage with profile result
#' result = computeAllParameterProfiles(...)
#' plotLikelihoodProfiles(result)
#'
#' # With options
#' plotLikelihoodProfiles(result,
#'   options = list(max_cols = 2, log_scale = TRUE, show_ci = FALSE))
#'
#' # With true values
#' plotLikelihoodProfiles(result, true_values = c(2.0, 1.0))
#'
#' @export
plotLikelihoodProfiles = function(profile_result, options = list(), true_values = NULL) {

  # Validate input
  if (!is.list(profile_result) || !all(c("profiles", "confidence_intervals", "options") %in% names(profile_result))) {
    stop("profile_result must be a list from computeAllParameterProfiles()")
  }

  # Set default options
  default_options = list(
    max_cols = 3,
    show_ci = TRUE,
    show_optimal = TRUE,
    show_threshold = TRUE,
    log_scale = FALSE,
    main_title = NULL,
    cex = 1.0,
    lwd = 1.5,
    mar = c(4, 4, 2, 1)
  )

  # Merge with user options
  plot_options = modifyList(default_options, options)

  # Get data
  profiles = profile_result$profiles
  confidence_intervals = profile_result$confidence_intervals
  profile_opts = profile_result$options
  n_params = length(profiles)

  # Get parameter names
  param_names = rownames(confidence_intervals)
  if (is.null(param_names)) {
    param_names = paste0("param_", 1:n_params)
  }

  # Validate true_values if provided
  if (!is.null(true_values)) {
    if (length(true_values) != n_params) {
      warning("Length of true_values does not match number of parameters. Ignoring true_values.")
      true_values = NULL
    }
  }

  # Store original par settings
  old_par = par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Setup grid layout
  n_cols = min(plot_options$max_cols, n_params)
  n_rows = ceiling(n_params / n_cols)

  # Set up plotting area
  par(mfrow = c(n_rows, n_cols), mar = plot_options$mar, cex = plot_options$cex)

  # Create main title if not provided
  if (is.null(plot_options$main_title)) {
    scale_text = if (plot_options$log_scale) "Log-Scale" else "Linear Scale"
    plot_options$main_title = sprintf("Profile Likelihood Results (%s, %s grid, %d points)",
                                      scale_text,
                                      profile_opts$grid_method,
                                      profile_opts$grid_points)
  }

  # Plot each parameter
  for (i in 1:n_params) {
    true_val = if (!is.null(true_values)) true_values[i] else NULL

    plotSingleProfile(
      profile_data = profiles[[i]],
      param_name = param_names[i],
      param_index = i,
      confidence_interval = confidence_intervals[i, ],
      profile_options = profile_opts,
      plot_options = plot_options,
      true_value = true_val,
      show_legend = (n_params == 1)
    )
  }

  # Add overall title if multiple plots
  if (n_params > 1) {
    mtext(plot_options$main_title, side = 3, line = -1, outer = TRUE, cex = plot_options$cex * 1.2)
  }

  invisible(old_par)
}

#' Plot Single Parameter Profile
#'
#' Internal function to plot a single parameter's profile likelihood.
#'
#' @param profile_data Profile data for single parameter
#' @param param_name Parameter name for plot title
#' @param param_index Parameter index
#' @param confidence_interval Confidence interval bounds [lower, upper]
#' @param profile_options Profile computation options
#' @param plot_options Plotting options
#' @param true_value Optional true parameter value
#' @param show_legend Whether to show legend on this plot
#'
#' @keywords internal
plotSingleProfile = function(profile_data, param_name, param_index, confidence_interval,
                             profile_options, plot_options, true_value = NULL, show_legend = FALSE) {

  # Check if profile computation failed
  if (profile_data$failed) {
    # Create empty plot with error message
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Parameter Value", ylab = "Likelihood",
         main = param_name, cex.main = plot_options$cex)
    text(0.5, 0.5, "Profile computation\nfailed",
         cex = plot_options$cex, col = "red", adj = 0.5)
    return(invisible())
  }

  # Extract data
  grid_values = profile_data$grid_values
  profile_costs = profile_data$profile_costs
  optimal_cost = profile_data$optimal_cost

  # Check for valid data
  if (length(grid_values) < 2) {
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Parameter Value", ylab = "Likelihood",
         main = param_name, cex.main = plot_options$cex)
    text(0.5, 0.5, "Insufficient\nvalid data",
         cex = plot_options$cex, col = "red", adj = 0.5)
    return(invisible())
  }

  # Filter out invalid values
  valid_indices = is.finite(profile_costs)
  if (sum(valid_indices) < 2) {
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Parameter Value", ylab = "Likelihood",
         main = param_name, cex.main = plot_options$cex)
    text(0.5, 0.5, "No valid\ndata points",
         cex = plot_options$cex, col = "red", adj = 0.5)
    return(invisible())
  }

  valid_grid = grid_values[valid_indices]
  valid_costs = profile_costs[valid_indices]

  # Sort for smooth plotting
  sorted_indices = order(valid_grid)
  sorted_grid = valid_grid[sorted_indices]
  sorted_costs = valid_costs[sorted_indices]

  # Transform to likelihood values
  cost_diff = sorted_costs - optimal_cost

  if (plot_options$log_scale) {
    # Log-likelihood difference
    y_values = -cost_diff
    y_label = "Log-Likelihood Difference"
    threshold_value = -profile_options$ll_ratio_threshold / 2
  } else {
    # Relative likelihood
    y_values = exp(-cost_diff)
    y_label = "Likelihood Ratio"
    threshold_value = exp(-profile_options$ll_ratio_threshold / 2)
  }

  # Set up plot limits
  x_range = range(sorted_grid)
  x_margin = diff(x_range) * 0.05
  xlim = c(x_range[1] - x_margin, x_range[2] + x_margin)

  y_range = range(y_values, na.rm = TRUE)
  if (plot_options$show_threshold) {
    y_range = range(c(y_range, threshold_value), na.rm = TRUE)
  }

  if (plot_options$log_scale) {
    y_margin = diff(y_range) * 0.1
    ylim = c(y_range[1] - y_margin, max(y_range[2] + y_margin, 0.1))
  } else {
    ylim = c(0, max(y_range[2] * 1.1, 1.0))
  }

  # Create base plot
  plot(sorted_grid, y_values, type = "l", lwd = plot_options$lwd, col = "blue",
       xlim = xlim, ylim = ylim, xlab = "Parameter Value", ylab = y_label,
       main = param_name, cex.main = plot_options$cex)

  # Add grid
  grid(col = "lightgray", lty = 2)

  # Re-plot the main curve on top of grid
  lines(sorted_grid, y_values, lwd = plot_options$lwd, col = "blue")

  # Add optimal parameter value
  if (plot_options$show_optimal) {
    # Find optimal parameter value (where cost is minimal)
    optimal_idx = which.min(sorted_costs)
    optimal_param = sorted_grid[optimal_idx]
    abline(v = optimal_param, col = "green", lwd = 2, lty = 2)
  }

  # Add true value if provided
  if (!is.null(true_value)) {
    abline(v = true_value, col = "black", lwd = 2, lty = 3)
  }

  # Add likelihood threshold
  if (plot_options$show_threshold) {
    abline(h = threshold_value, col = "red", lwd = 1, lty = 2)
  }

  # Add confidence interval bounds
  if (plot_options$show_ci && !any(is.na(confidence_interval))) {
    abline(v = confidence_interval[1], col = "magenta", lwd = 1, lty = 2)
    abline(v = confidence_interval[2], col = "magenta", lwd = 1, lty = 2)
  }

  # Add legend for single plot only
  if (show_legend) {
    legend_items = c("Profile")
    legend_colors = c("blue")
    legend_lty = c(1)
    legend_lwd = c(plot_options$lwd)

    if (plot_options$show_optimal) {
      legend_items = c(legend_items, "Optimal")
      legend_colors = c(legend_colors, "green")
      legend_lty = c(legend_lty, 2)
      legend_lwd = c(legend_lwd, 2)
    }

    if (!is.null(true_value)) {
      legend_items = c(legend_items, "True Value")
      legend_colors = c(legend_colors, "black")
      legend_lty = c(legend_lty, 3)
      legend_lwd = c(legend_lwd, 2)
    }

    if (plot_options$show_threshold) {
      legend_items = c(legend_items, "Threshold")
      legend_colors = c(legend_colors, "red")
      legend_lty = c(legend_lty, 2)
      legend_lwd = c(legend_lwd, 1)
    }

    if (plot_options$show_ci && !any(is.na(confidence_interval))) {
      legend_items = c(legend_items, "CI Bounds")
      legend_colors = c(legend_colors, "magenta")
      legend_lty = c(legend_lty, 2)
      legend_lwd = c(legend_lwd, 1)
    }

    legend("topright", legend = legend_items, col = legend_colors,
           lty = legend_lty, lwd = legend_lwd, cex = plot_options$cex * 0.8)
  }
}
