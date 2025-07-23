#' likelyprofiler: Profile Likelihood Confidence Intervals
#'
#' Compute profile likelihood confidence intervals for parameters using 
#' grid-based optimization with customizable methods and parallel computation support.
#' Get started quickly by checking vignette \code{vignette("gaussian_mixture_profiling", package = "likelyprofiler")}.
#'
#' Profile likelihood is a robust method for uncertainty quantification that 
#' provides confidence intervals for individual parameters by fixing each parameter 
#' at different values and optimizing over the remaining parameters. This approach 
#' is particularly valuable when asymptotic approximations may be unreliable or 
#' when the likelihood surface is non-quadratic.
#'
#' The main function is \code{\link{computeLikelihoodProfiles}}, which computes 
#' profile likelihood confidence intervals for all parameters of a given 
#' likelihood function. The package supports multiple optimization algorithms 
#' including base R's \code{\link[stats]{optim}}, \code{\link[DEoptim]{DEoptim}}, and custom user-defined optimizers.
#'
#' @section Main Functions:
#' * \code{\link{computeLikelihoodProfiles}}: Compute profile likelihood confidence 
#'   intervals for all parameters using grid-based optimization. Supports parallel 
#'   computation and multiple optimization algorithms.
#' * \code{\link{plotLikelihoodProfiles}}: Plot the computed profile likelihood 
#'   results to visualize confidence intervals and likelihood surfaces.
#'
#' @section Key Features:
#' * **Multiple optimizers**: Supports \code{\link[stats]{optim}}, \code{\link[DEoptim]{DEoptim}}, and custom optimization functions.
#' * **Grid methods**: Uniform and adaptive grid spacing for parameter exploration
#' * **Parallel computation**: Optional cluster-based parallel execution for faster computation
#' * **Robust validation**: Comprehensive input validation and error handling
#' * **Flexible configuration**: Customizable grid parameters, optimization options, and confidence levels
#'
#' @section Basic Usage:
#' The typical workflow involves:
#' 1. Define your negative log-likelihood function
#' 2. Set parameter bounds and current parameter values
#' 3. Call \code{\link{computeLikelihoodProfiles}} with your configuration
#' 4. Visualize results with \code{\link{plotLikelihoodProfiles}}
#'
#' For a complete example, see: \code{vignette("gaussian_mixture_profiling", package = "likelyprofiler")}
#'
#' @examples
#' # Define a simple quadratic likelihood function
#' quadratic_likelihood = function(params) {
#'   x = params[1]
#'   y = params[2]
#'   (x - 2)^2 + 2 * (y - 1)^2
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
#'   profile_options = list(grid_points = 15),
#'   verbose = FALSE
#' )
#'
#' # View results
#' print(result$summary)
#' print(result$confidence_intervals)
#'
#' @aliases likelyprofiler
#' @name likelyprofiler-package
"_PACKAGE"
