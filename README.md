# likelyprofiler

Profile likelihood confidence intervals for parameter uncertainty quantification.

This package provides robust methods for computing profile likelihood confidence intervals using grid-based optimization. Profile likelihood is particularly valuable when asymptotic approximations may be unreliable or when the likelihood surface is non-quadratic. The package supports multiple optimization algorithms including base R's `optim`, `DEoptim`, and custom user-defined optimizers, with optional parallel computation for faster execution.

## Version
The package was tested using R version 4.5.0 (2025-04-11).

## Installation

The following must be installed prior to installing the package:
* [R](https://www.r-project.org/)
* [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or [Xcode](https://mac.install.guide/commandlinetools/) on Mac.
* `remotes` package. Install by running `install.packages("remotes")` in R. Using larger package development `devtools` will also work.
* [git](https://git-scm.com/)

Then run the following command in R:
```R
remotes::install_git(url = "https://gitlab.com/csb.ethz/likelyprofiler", build_vignettes = TRUE, dependencies = TRUE)
```
Alternatively, you can set `build_vignettes = FALSE` if you do not wish to use the vignettes for documentation.

And attach the package with

```R
library(likelyprofiler)
```

Type `?likelyprofiler` for the package documentation, including list of functions.

## Examples

Get started with the main vignette:
```R
vignette("gaussian_mixture_profiling", package = "likelyprofiler")
```

Basic usage example:
```R
# Define a simple quadratic likelihood function
quadratic_likelihood <- function(params) {
  x <- params[1]
  y <- params[2]
  (x - 2)^2 + 2 * (y - 1)^2
}

# Set up the problem
true_params <- c(2.0, 1.0)
names(true_params) <- c("param_x", "param_y")
bounds <- list(lower = c(0, 0), upper = c(4, 3))

# Compute profile likelihood
result <- computeLikelihoodProfiles(
  params_current = true_params,
  negLogLikelihood = quadratic_likelihood,
  bounds = bounds,
  profile_options = list(grid_points = 15),
  verbose = FALSE
)

# View results
print(result$summary)
print(result$confidence_intervals)

# Plot results
plotLikelihoodProfiles(result)
```
