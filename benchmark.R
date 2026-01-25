# SMCO Benchmark Script
# Version: 1.1.0
# Release Date: January 25, 2026
#
# Strategic Monte Carlo Optimization (SMCO) Algorithm
# By: Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, and Guodong Zhang
#
# Citation: X. Chen, Z. Chen, W.Y. Gao, X. Yan, & G. Zhang, Optimization via the strategic law
#   of large numbers, Proc. Natl. Acad. Sci. U.S.A. 123 (4) e2519845123,
#   https://doi.org/10.1073/pnas.2519845123 (2026).
#
# This script benchmarks SMCO on standard test functions.
# Users can configure the test function, dimension, and SMCO parameters below.

#####################################
# USER CONFIGURATION
#####################################

# Test function selection
# Available functions (see testfuncs.R for details):
#   Scalable (any dimension): "Rastrigin", "Ackley", "Griewank", "Rosenbrock",
#                             "Michalewicz", "DixonPrice", "Zakharov", "Qing",
#                             "SqNorm", "-SqNorm", "ANN"
#   Fixed 2D only: "Dropwave", "Shubert", "McCormick", "Bukin6", "Damavandi",
#                  "Mishra06", "EggHolder", "Six-Hump Camel"

name_config <- "Rastrigin"  # Test function name
dim_config <- 10            # Dimension (ignored for fixed 2D functions)

# Benchmark settings
n_replications <- 10        # Number of Monte Carlo replications
random_seed <- 123          # Random seed for reproducibility

# SMCO configuration
SMCO_options <- list(
  n_starts = 32,            # Number of starting points
  iter_max = 200,           # Maximum iterations per start
  bounds_buffer = 0.05,     # Bounds extension factor
  buffer_rand = TRUE,       # Randomize bounds extension
  tol_conv = 1e-8,          # Convergence tolerance
  refine_search = TRUE,     # Enable refinement phase
  refine_ratio = 0.5,       # Fraction of iterations for refinement
  partial_option = "center",# "center" (two-sided) or "forward" (one-sided)
  use_runmax = TRUE,        # Track running maximum
  use_parallel = FALSE,     # Enable parallel execution
  verbose = FALSE           # Print progress
)

# ANN-specific configuration (only used if name_config = "ANN")
ANN_options <- list(
  input_dim = 3,            # Number of input features

  hidden_dim = 5,           # Number of hidden neurons
  n = 500,                  # Number of training samples
  use_truth = FALSE,        # Use noiseless target (TRUE) or noisy (FALSE)
  noise_sigma = 1           # Noise standard deviation
)

#####################################
# END USER CONFIGURATION
#####################################

cat("==============================================\n")
cat("SMCO Benchmark\n")
cat("==============================================\n\n")

# Load required packages
if (!requireNamespace("qrng", quietly = TRUE)) {
  stop("Package 'qrng' is required. Install with: install.packages('qrng')")
}
library(qrng)

# Source SMCO algorithm and test functions
source("SMCO.R")
source("testfuncs.R")

# Helper function: generate uniform starting points
generate_unif_starts <- function(n_starts, bounds_lower, bounds_upper) {
  d <- length(bounds_lower)
  result <- matrix(NA, nrow = n_starts, ncol = d)
  for (j in 1:d) {
    result[, j] <- bounds_lower[j] + runif(n_starts) * (bounds_upper[j] - bounds_lower[j])
  }
  return(result)
}

# Get test function configuration
if (name_config == "ANN") {
  test_config <- assign_config(name = name_config, d = dim_config, opt_config = ANN_options)
} else {
  test_config <- assign_config(name = name_config, d = dim_config, opt_config = NULL)
}

bounds_lower <- test_config$bounds_lower
bounds_upper <- test_config$bounds_upper
test_func <- test_config$f
test_func_minus <- function(x) (-test_func(x))

# Print configuration
cat("Configuration:\n")
cat("  Test function:", name_config, "\n")
cat("  Dimension:", length(bounds_lower), "\n")
cat("  Bounds: [", bounds_lower[1], ",", bounds_upper[1], "]^d\n")
cat("  Replications:", n_replications, "\n")
cat("  n_starts:", SMCO_options$n_starts, "\n")
cat("  iter_max:", SMCO_options$iter_max, "\n")
cat("  True minimum:", ifelse(is.na(test_config$f_min), "unknown", test_config$f_min), "\n")
cat("  True maximum:", ifelse(is.na(test_config$f_max), "unknown", test_config$f_max), "\n\n")

# Run benchmark for minimization
run_benchmark <- function(to_maximize) {
  direction <- if (to_maximize) "MAXIMIZATION" else "MINIMIZATION"
  cat("Running", direction, "benchmark...\n")

  if (to_maximize) {
    test_func_opt <- test_func
    true_opt <- test_config$f_max
  } else {
    test_func_opt <- test_func_minus
    true_opt <- test_config$f_min
  }

  time_results <- numeric(n_replications)
  fopt_results <- numeric(n_replications)

  set.seed(random_seed)

  for (i in 1:n_replications) {
    if (i %% max(1, floor(n_replications / 5)) == 1 || i == n_replications) {
      cat("  Replication", i, "of", n_replications, "\n")
    }

    start_points <- generate_unif_starts(SMCO_options$n_starts, bounds_lower, bounds_upper)

    start_time <- Sys.time()
    results <- SMCO_multi(f = test_func_opt,
                          bounds_lower, bounds_upper,
                          start_points = start_points,
                          opt_control = SMCO_options)
    end_time <- Sys.time()

    time_results[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    fopt <- results$best_result$f_optimal
    fopt_results[i] <- ifelse(to_maximize, fopt, -fopt)
  }

  # Calculate statistics
  time_avg <- mean(time_results)
  time_sd <- sd(time_results)
  fopt_mean <- mean(fopt_results)
  fopt_sd <- sd(fopt_results)

  if (!is.na(true_opt)) {
    AE <- abs(fopt_results - true_opt)
    rMSE <- sqrt(mean(AE^2))
    AE_median <- median(AE)
    AE_95 <- quantile(AE, 0.95)
  } else {
    rMSE <- NA
    AE_median <- NA
    AE_95 <- NA
  }

  return(list(
    direction = direction,
    true_opt = true_opt,
    time_avg = time_avg,
    time_sd = time_sd,
    fopt_mean = fopt_mean,
    fopt_sd = fopt_sd,
    rMSE = rMSE,
    AE_median = AE_median,
    AE_95 = AE_95
  ))
}

# Run benchmarks
results_min <- run_benchmark(to_maximize = FALSE)
results_max <- run_benchmark(to_maximize = TRUE)

# Print results
cat("\n==============================================\n")
cat("RESULTS\n")
cat("==============================================\n\n")

print_result <- function(res) {
  cat(res$direction, "\n")
  cat("  True optimum:     ", ifelse(is.na(res$true_opt), "unknown", round(res$true_opt, 6)), "\n")
  cat("  Found (mean/sd):  ", round(res$fopt_mean, 6), " / ", round(res$fopt_sd, 6), "\n")
  if (!is.na(res$rMSE)) {
    cat("  rMSE:             ", round(res$rMSE, 6), "\n")
    cat("  AE median:        ", round(res$AE_median, 6), "\n")
    cat("  AE 95th pctl:     ", round(res$AE_95, 6), "\n")
  }
  cat("  Time (mean/sd):   ", round(res$time_avg, 3), "s / ", round(res$time_sd, 3), "s\n\n")
}

print_result(results_min)
print_result(results_max)

cat("Benchmark complete.\n")
