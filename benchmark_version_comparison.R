# Benchmark Comparison: SMCO v1.0.0 vs v1.1.0
# Version: 1.1.0
# Release Date: January 25, 2026
#
# Strategic Monte Carlo Optimization (SMCO) Algorithm
# By: Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, and Guodong Zhang
#
# This script compares performance between v1.0.0 and v1.1.0 implementations.
# Results are used to document performance improvements in README.md.

cat("==============================================\n")
cat("SMCO Version Comparison: v1.0.0 vs v1.1.0\n")
cat("==============================================\n\n")

#####################################
# Load both versions into separate environments
#####################################

# Load v1.0.0 into isolated environment
env_v100 <- new.env()
source("v1.0.0/SMCO.R", local = env_v100)

# Load v1.1.0 into isolated environment
env_v110 <- new.env()
source("SMCO.R", local = env_v110)

# Load test functions
source("testfuncs.R")

#####################################
# Configuration
#####################################

# Test dimensions
dimensions <- c(2, 10, 50)

# Variants to test (mapped to opt_control settings)
variants <- list(
  SMCO = list(refine_search = FALSE, iter_boost = 0),
  SMCO_R = list(refine_search = TRUE, iter_boost = 0, refine_ratio = 0.5),
  SMCO_BR = list(refine_search = TRUE, iter_boost = 1000, refine_ratio = 0.5)
)

# Common settings
common_settings <- list(
  n_starts = 10,
  iter_max = 200,
  bounds_buffer = 0.05,
  tol_conv = 1e-8,
  use_runmax = TRUE,
  use_parallel = FALSE,
  seed = 123
)

# Number of replications for timing
n_reps <- 3

#####################################
# Test function: Rastrigin (minimization via negation)
#####################################

run_comparison <- function(dim, variant_name, variant_settings, n_reps = 3) {
  # Set up test function (Rastrigin - minimize)
  config <- config_rast(dim)
  f_min <- function(x) -config$f(x)  # Negate for minimization (SMCO maximizes)

  bounds_lower <- config$bounds_lower
  bounds_upper <- config$bounds_upper

  # Generate shared starting points using v1.1.0's generator
  set.seed(123)
  start_points <- env_v110$generate_sobol_points(common_settings$n_starts,
                                                  bounds_lower, bounds_upper,
                                                  seed = 123)

  # Build opt_control for each version
  opt_control_v100 <- c(common_settings, variant_settings)
  opt_control_v100$iter_nstart <- common_settings$n_starts

  opt_control_v110 <- c(common_settings, variant_settings)

  # Run v1.0.0 with timing
  times_v100 <- numeric(n_reps)
  results_v100 <- list()
  for (i in 1:n_reps) {
    t_start <- Sys.time()
    results_v100[[i]] <- env_v100$SMCO_multi(f_min, bounds_lower, bounds_upper,
                                              start_points = start_points,
                                              opt_control = opt_control_v100)
    t_end <- Sys.time()
    times_v100[i] <- as.numeric(difftime(t_end, t_start, units = "secs"))
  }

  # Run v1.1.0 with timing
  times_v110 <- numeric(n_reps)
  results_v110 <- list()
  for (i in 1:n_reps) {
    t_start <- Sys.time()
    results_v110[[i]] <- env_v110$SMCO_multi(f_min, bounds_lower, bounds_upper,
                                              start_points = start_points,
                                              opt_control = opt_control_v110)
    t_end <- Sys.time()
    times_v110[i] <- as.numeric(difftime(t_end, t_start, units = "secs"))
  }

  # Compute statistics
  mean_time_v100 <- mean(times_v100)
  mean_time_v110 <- mean(times_v110)
  speedup <- (mean_time_v100 - mean_time_v110) / mean_time_v100 * 100

  # Get best f_optimal (negated back to original scale for reporting)
  f_opt_v100 <- -results_v100[[1]]$best_result$f_optimal
  f_opt_v110 <- -results_v110[[1]]$best_result$f_optimal

  list(
    dim = dim,
    variant = variant_name,
    time_v100 = mean_time_v100,
    time_v110 = mean_time_v110,
    speedup_pct = speedup,
    f_opt_v100 = f_opt_v100,
    f_opt_v110 = f_opt_v110,
    f_opt_diff = f_opt_v110 - f_opt_v100
  )
}

#####################################
# Run all comparisons
#####################################

results <- list()
idx <- 1

for (dim in dimensions) {
  cat(sprintf("\n--- Testing dimension %d ---\n", dim))

  for (variant_name in names(variants)) {
    cat(sprintf("  Variant: %s ... ", variant_name))

    result <- run_comparison(dim, variant_name, variants[[variant_name]], n_reps)
    results[[idx]] <- result
    idx <- idx + 1

    cat(sprintf("Speedup: %.1f%%\n", result$speedup_pct))
  }
}

#####################################
# Print summary table
#####################################

cat("\n\n==============================================\n")
cat("SUMMARY: Performance Comparison\n")
cat("==============================================\n\n")

cat(sprintf("%-10s %-10s %12s %12s %10s %14s %14s\n",
            "Dimension", "Variant", "v1.0.0 (s)", "v1.1.0 (s)", "Speedup",
            "f_opt v1.0.0", "f_opt v1.1.0"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (r in results) {
  cat(sprintf("%-10d %-10s %12.3f %12.3f %9.1f%% %14.4f %14.4f\n",
              r$dim, r$variant, r$time_v100, r$time_v110, r$speedup_pct,
              r$f_opt_v100, r$f_opt_v110))
}

#####################################
# Compute overall statistics
#####################################

speedups <- sapply(results, function(r) r$speedup_pct)
mean_speedup <- mean(speedups)
max_speedup <- max(speedups)

# Speedup by dimension
cat("\n\nSpeedup by Dimension:\n")
for (dim in dimensions) {
  dim_results <- results[sapply(results, function(r) r$dim == dim)]
  dim_speedups <- sapply(dim_results, function(r) r$speedup_pct)
  cat(sprintf("  %d-D: %.1f%% average speedup\n", dim, mean(dim_speedups)))
}

cat(sprintf("\nOverall: %.1f%% average speedup (max %.1f%%)\n", mean_speedup, max_speedup))

#####################################
# Key optimizations in v1.1.0
#####################################

cat("\n\n==============================================\n")
cat("KEY OPTIMIZATIONS IN v1.1.0\n")
cat("==============================================\n")
cat("
1. Pre-computed h_step matrix: Avoids repeated division in main loop
2. In-place vector modification: Reduces memory allocation in compute_partial_signs
3. Metadata-based runmax tracking: Single vector copy at end instead of per-iteration
4. Pre-computed fixed bounds: When buffer_rand=FALSE, bounds computed once
5. Vectorized Sobol scaling: Uses R's column-wise recycling
6. Input validation: New validate_smco_inputs() for better error messages
7. Adaptive defaults: n_starts defaults to max(5, sqrt(dim))
")

cat("\n==============================================\n")
cat("Benchmark completed.\n")
cat("==============================================\n")
