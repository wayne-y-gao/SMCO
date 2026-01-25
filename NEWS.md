# NEWS

## SMCO v1.1.0 (January 25, 2026)

### Performance Optimizations

This release includes significant performance improvements, particularly for high-dimensional optimization problems.

#### Changes in `compute_partial_signs()`
- **In-place vector modification**: Instead of creating copies of the input vector for each perturbation, the function now modifies the vector in-place and restores original values after each dimension evaluation
- **Metadata-based runmax tracking**: Instead of copying the full d-dimensional vector on each improvement, the function now stores only the dimension index and perturbed value, reconstructing the best point once at the end

#### Changes in `SMCO_single()`
- **Pre-computed h_step lookup table**: For `iter_max < 1000`, all step sizes are pre-computed in a matrix before the main loop, eliminating repeated division operations
- **Cached fixed bounds**: When `buffer_rand = FALSE`, the extended bounds are computed once before the loop instead of every iteration
- **Pre-computed convergence threshold**: The convergence check threshold is computed once before the loop

#### Changes in `generate_sobol_points()`
- **Vectorized scaling**: Replaced the dimension-wise loop with a vectorized operation using R's column-wise recycling

### Performance Benchmark Results

Benchmark on 50-dimensional Rastrigin function (20 replications, 32 starting points, 300 iterations):

#### Maximization

| Algorithm | v1.0.0 | v1.1.0 | Speedup |
|-----------|--------|--------|---------|
| SMCO      | 3.108s | 2.794s | **10.1% faster** |
| SMCO_R    | 3.106s | 2.831s | **8.9% faster** |
| SMCO_BR   | 3.620s | 2.864s | **20.9% faster** |

#### Minimization

| Algorithm | v1.0.0 | v1.1.0 | Speedup |
|-----------|--------|--------|---------|
| SMCO      | 5.272s | 2.983s | **43.4% faster** |
| SMCO_R    | 5.474s | 3.015s | **44.9% faster** |
| SMCO_BR   | 5.140s | 3.060s | **40.5% faster** |

**Note:** Performance gains scale with problem dimension due to reduced memory allocation overhead from eliminated vector copies.

### Correctness Verification

All accuracy metrics (rMSE, AE50, AE95, AE99) are identical between v1.0.0 and v1.1.0, confirming that the optimizations do not affect numerical results.

### Repository Streamlining

- **Merged `testfunc_NN1L.R` into `testfuncs.R`**: All test functions (deterministic and neural network) are now in a single file
- **Added `benchmark.R`**: New user-configurable benchmark script for testing SMCO
- **Removed comparison algorithms**: `GD.R`, `SignGD.R`, `SPSA.R`, and `RunComparison.R` are no longer included (available in v1.0.0 archive)

---

## SMCO v1.0.0 (July 16, 2025)

Initial public release of the SMCO algorithm.

### Features
- Core SMCO algorithm with multi-start optimization
- Refinement search phase for local optimization
- Boosted search with configurable starting index
- Support for parallel execution
- Comprehensive benchmarking framework
- Deterministic and neural network-based test functions
- Comparison implementations: Gradient Descent, Sign Gradient Descent, SPSA
