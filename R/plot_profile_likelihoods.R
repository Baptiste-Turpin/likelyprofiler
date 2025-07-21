#' Plot Profile Likelihood Results
#'
#' Create plots of profile likelihood results using base R graphics.
#' Displays profile likelihood curves for all parameters in a grid layout.
#'
#' @param profile_result Result from computeLikelihoodProfiles()
#' @param options List of general plotting options (see Details)
#' @param plot_options List of parameters passed to \code{\link[graphics]{par}} (see Details)
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
#'   \item \code{lwd}: Line width for curves (default: 1.5)
#'   \item \code{line_x}: Parameter \code{line} in \code{\link[graphics]{mtext}} for x-axis label (default: 0.5)
#'   \item \code{line_y}: Parameter \code{line} in \code{\link[graphics]{mtext}} for y-axis label (default: 1.5)
#'   \item \code{cex.lab}: Relative size of axis labels (default: NULL)
#' }
#'
#' The \code{plot_options} list contains parameters passed to \code{\link[graphics]{par}}:
#' \itemize{
#'   \item \code{cex}: Character expansion factor for text (default: 1.0)
#'   \item \code{mar}: Plot margins (default: c(1.5, 0.25, 1.2, 0.25))
#'   \item \code{oma}: Outer margin parameters (default: c(1.7, 2.5, 0, 7))
#'   \item \code{cex.axis}: Relative size of axis tick mark labels (default: 0.8)
#'   \item \code{mgp}: Margin line positions for axis title, labels, and line (default: c(3, 0.5, 0))
#' }
#'
#' @return Invisibly returns the original plotting parameters (for restoration)
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
#' # Compute profile likelihood
#' result = computeLikelihoodProfiles(
#'   params_current = true_params,
#'   negLogLikelihood = quadratic_likelihood,
#'   bounds = bounds,
#'   profile_options = list(
#'     grid_method = "adaptive",
#'     grid_points = 15,
#'     max_grid_range_multiplier = 1.0
#'   ),
#'   optimizer = "optim",
#'   verbose = FALSE,
#'   method = "L-BFGS-B"
#' )
#'
#' # Basic plotting
#' plotLikelihoodProfiles(result)
#'
#' # With options and true values
#' plotLikelihoodProfiles(result,
#'   options = list(max_cols = 2, log_scale = TRUE, show_ci = FALSE),
#'   plot_options = list(cex = 0.9, cex.axis = 0.7),
#'   true_values = true_params)
#'
#' @export
plotLikelihoodProfiles = function(profile_result, options = list(), plot_options = list(), true_values = NULL) {

  # Validate input
  if (!is.list(profile_result) || !all(c("profiles", "confidence_intervals", "options") %in% names(profile_result))) {
    stop("profile_result must be a list from computeLikelihoodProfiles()")
  }

  # Set default options
  default_options = list(
    max_cols = 3,
    show_ci = TRUE,
    show_optimal = TRUE,
    show_threshold = TRUE,
    log_scale = FALSE,
    lwd = 1.5,
    line_x = 0.5,
    line_y = 1.5,
    cex.lab = NULL
  )

  # Set default plot options (par parameters)
  default_plot_options = list(
    cex = 1.0,
    mar = c(1.5, 0.25, 1.2, 0.25),
    oma = c(1.7, 2.5, 0, 7),
    cex.axis = 0.8,
    mgp = c(3, 0.5, 0)
  )

  # Merge with user options
  general_options = utils::modifyList(default_options, options)
  par_options = utils::modifyList(default_plot_options, plot_options)

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

  # Setup grid layout with outer margins
  n_cols = min(general_options$max_cols, n_params)
  n_rows = ceiling(n_params / n_cols)

  # Set up plotting area with outer margins
  par_args = c(list(mfrow = c(n_rows, n_cols)), par_options)
  old_par = do.call(graphics::par, par_args)
  on.exit(graphics::par(old_par))

  # Plot each parameter
  for (i in 1:n_params) {
    # Calculate grid position
    col = ((i - 1) %% n_cols) + 1
    is_leftmost_col = (col == 1)

    true_val = if (!is.null(true_values)) true_values[i] else NULL

    plotSingleProfile(
      profile_data = profiles[[i]],
      param_name = param_names[i],
      param_index = i,
      confidence_interval = confidence_intervals[i, ],
      profile_options = profile_opts,
      general_options = general_options,
      par_options = par_options,
      true_value = true_val,
      is_leftmost_col = is_leftmost_col
    )
  }

  # Add common labels
  y_label = if (general_options$log_scale) "Log-Likelihood Difference" else "Likelihood Ratio"
  graphics::mtext("Parameter Value", side = 1, line = general_options$line_x, outer = TRUE, adj = 0.5, cex = general_options$cex.lab)
  graphics::mtext(y_label, side = 2, line = general_options$line_y, outer = TRUE, adj = 0.5, cex = general_options$cex.lab)

  # Add common legend on the right side
  addProfileLegend(general_options, par_options, true_values)

}

#' Add Common Profile Legend
#'
#' Internal function to add a common legend on the right side of the plot.
#'
#' @param general_options General plotting options
#' @param par_options Parameters for par()
#' @param true_values Optional true parameter values
#'
#' @keywords internal
addProfileLegend = function(general_options, par_options, true_values) {

  # Build legend items
  legend_items = c("Profile")
  legend_colors = c("blue")
  legend_lty = c(1)
  legend_lwd = c(general_options$lwd)

  if (general_options$show_optimal) {
    legend_items = c(legend_items, "Optimized")
    legend_colors = c(legend_colors, "green")
    legend_lty = c(legend_lty, 2)
    legend_lwd = c(legend_lwd, 2)
  }

  if (!is.null(true_values)) {
    legend_items = c(legend_items, "True Value")
    legend_colors = c(legend_colors, "black")
    legend_lty = c(legend_lty, 3)
    legend_lwd = c(legend_lwd, 2)
  }

  if (general_options$show_threshold) {
    legend_items = c(legend_items, "Threshold")
    legend_colors = c(legend_colors, "red")
    legend_lty = c(legend_lty, 2)
    legend_lwd = c(legend_lwd, 1)
  }

  if (general_options$show_ci) {
    legend_items = c(legend_items, "CI Bounds")
    legend_colors = c(legend_colors, "magenta")
    legend_lty = c(legend_lty, 2)
    legend_lwd = c(legend_lwd, 1)
  }

  # Reset par() settings and create a new plot area for legend
  # Similar to plotDistrib approach
  old_par_legend = graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(graphics::par(old_par_legend))

  # Create invisible plot that covers the entire figure
  graphics::plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')

  # Add legend in the right side
  graphics::legend("right",
         legend = legend_items,
         col = legend_colors,
         lty = legend_lty,
         lwd = legend_lwd,
         cex = par_options$cex * 0.9,
         bg = "white",
         box.col = "black")
}#' Plot Single Parameter Profile
#'
#' Internal function to plot a single parameter's profile likelihood.
#'
#' @param profile_data Profile data for single parameter
#' @param param_name Parameter name for plot title
#' @param param_index Parameter index
#' @param confidence_interval Confidence interval bounds [lower, upper]
#' @param profile_options Profile computation options
#' @param general_options General plotting options
#' @param par_options Parameters for par()
#' @param true_value Optional true parameter value
#' @param is_leftmost_col Logical indicating if this is leftmost column
#'
#' @keywords internal
plotSingleProfile = function(profile_data, param_name, param_index, confidence_interval,
                             profile_options, general_options, par_options, true_value = NULL, is_leftmost_col = TRUE) {

  # Check if profile computation failed
  if (profile_data$failed) {
    # Create empty plot with error message
    graphics::plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", yaxt = "n",
         main = param_name, cex.main = par_options$cex)
    graphics::text(0.5, 0.5, "Profile computation\nfailed",
         cex = par_options$cex, col = "red", adj = 0.5)
    return(invisible())
  }

  # Extract data
  grid_values = profile_data$grid_values
  profile_costs = profile_data$profile_costs
  optimal_cost = profile_data$optimal_cost
  optimized_value = profile_data$optimal_param_value


  # Check for valid data
  if (length(grid_values) < 2) {
    graphics::plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", yaxt = "n",
         main = param_name, cex.main = par_options$cex)
    graphics::text(0.5, 0.5, "Insufficient\nvalid data",
         cex = par_options$cex, col = "red", adj = 0.5)
    return(invisible())
  }

  # Filter out invalid values
  valid_indices = is.finite(profile_costs)
  if (sum(valid_indices) < 2) {
    graphics::plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", yaxt = "n",
         main = param_name, cex.main = par_options$cex)
    graphics::text(0.5, 0.5, "No valid\ndata points",
         cex = par_options$cex, col = "red", adj = 0.5)
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

  if (general_options$log_scale) {
    # Log-likelihood difference
    y_values = -cost_diff
    threshold_value = -profile_options$ll_ratio_threshold / 2
  } else {
    # Relative likelihood
    y_values = exp(-cost_diff)
    threshold_value = exp(-profile_options$ll_ratio_threshold / 2)
  }

  # Set up plot limits
  x_range = range(sorted_grid)
  x_margin = diff(x_range) * 0.05
  xlim = c(x_range[1] - x_margin, x_range[2] + x_margin)

  # Y-axis limits: common [0,1] for natural space, individual for log-space
  if (general_options$log_scale) {
    y_range = range(y_values, na.rm = TRUE)
    if (general_options$show_threshold) {
      y_range = range(c(y_range, threshold_value), na.rm = TRUE)
    }
    y_margin = diff(y_range) * 0.1
    ylim = c(y_range[1] - y_margin, max(y_range[2] + y_margin, 0.1))
    yaxt_setting = "s"  # Always show y-axis in log space
  } else {
    ylim = c(0, 1)  # Common limits for natural space
    yaxt_setting = ifelse(is_leftmost_col, "s", "n")
  }

  # Create base plot
  graphics::plot(NULL,
       xlim = xlim, ylim = ylim, xlab = "", ylab = "",
       yaxt = yaxt_setting, main = param_name, cex.main = par_options$cex)

  # Add grid
  graphics::grid(col = "lightgray", lty = 2)

  # Plot profile likelihood curve
  graphics::lines(sorted_grid, y_values, lwd = general_options$lwd, col = "blue")

  # Add optimal parameter value
  if (general_options$show_optimal) {
    graphics::abline(v = optimized_value, col = "green", lwd = 2, lty = 2)
  }

  # Add true value if provided
  if (!is.null(true_value)) {
    graphics::abline(v = true_value, col = "black", lwd = 2, lty = 3)
  }

  # Add likelihood threshold
  if (general_options$show_threshold) {
    graphics::abline(h = threshold_value, col = "red", lwd = 1, lty = 2)
  }

  # Add confidence interval bounds
  if (general_options$show_ci && !any(is.na(confidence_interval))) {
    graphics::abline(v = confidence_interval[1], col = "magenta", lwd = 1, lty = 2)
    graphics::abline(v = confidence_interval[2], col = "magenta", lwd = 1, lty = 2)
  }
}
