# Strategic Monte Carlo Optimization (SMCO)

R implementation of the SMCO algorithm from:

**Paper:** "Optimization via Strategic Law of Large Numbers"

**Authors:** Xiaohong Chen, Zengjing Chen, Wayne Yuan Gao, Xiaodong Yan, and Guodong Zhang

**Version:** 1.1.0

**Release Date:** January 25, 2026

**Repository maintained by:** Wayne Yuan Gao (waynegao@upenn.edu)

---

## Overview

SMCO is a derivative-free optimization algorithm based on the Strategic Law of Large Numbers. It is designed for black-box optimization problems where gradient information is unavailable or unreliable.

### Key Features
- **Derivative-free:** Only requires function evaluations
- **Multi-start:** Automatic exploration from multiple starting points
- **Refinement search:** Two-phase optimization with optional local refinement
- **Performance optimized:** v1.1.0 includes significant speedups (~45% faster in high dimensions)

---

## Algorithm Variants

SMCO provides three variants with increasing sophistication:

| Variant | Function | Description | Best For |
|---------|----------|-------------|----------|
| Basic | `SMCO()` | Single-phase optimization | Fast exploration |
| Refinement | `SMCO_R()` | Two-phase with local refinement | Balanced speed/accuracy |
| Boosted | `SMCO_BR()` | Regular + boosted searches | Maximum accuracy |

### Quick Usage
```r
source("SMCO.R")

# Define objective function (SMCO maximizes by default)
f <- function(x) -sum(x^2)
bounds_lower <- c(-5, -5)
bounds_upper <- c(5, 5)

# Choose your variant:
result <- SMCO(f, bounds_lower, bounds_upper)      # Basic - fastest
result <- SMCO_R(f, bounds_lower, bounds_upper)    # Refinement - balanced
result <- SMCO_BR(f, bounds_lower, bounds_upper)   # Boosted - most accurate
```

All wrapper functions accept additional parameters via `...` (e.g., `n_starts`, `iter_max`, `verbose`).

---

## Installation

SMCO is a pure R implementation with no compilation required.

### Dependencies
```r
install.packages("qrng")      # Required: Sobol sequence generation
# Note: 'compiler' package is pre-installed with R
```

---

## Quick Start

```r
# Load SMCO
source("SMCO.R")

# Define objective function (SMCO maximizes by default)
f <- function(x) {
  -sum(x^2)  # Simple quadratic, optimum at origin
}

# Define search bounds
bounds_lower <- c(-5, -5)
bounds_upper <- c(5, 5)

# Run optimization using a wrapper function
result <- SMCO_R(f, bounds_lower, bounds_upper, verbose = TRUE)

# Get results
print(result$best_result$x_optimal)  # Optimal point
print(result$best_result$f_optimal)  # Optimal value

# Or use SMCO_multi for full control
result <- SMCO_multi(f, bounds_lower, bounds_upper,
                     opt_control = list(
                       n_starts = 100,
                       iter_max = 500,
                       verbose = TRUE
                     ))
```

---

## API Reference

### Wrapper Functions

**`SMCO(f, bounds_lower, bounds_upper, start_points = NULL, ...)`**
- Basic single-phase optimization
- Sets `refine_search = FALSE`, `iter_boost = 0`

**`SMCO_R(f, bounds_lower, bounds_upper, start_points = NULL, ...)`**
- Two-phase optimization with refinement
- Sets `refine_search = TRUE`, `refine_ratio = 0.5`, `iter_boost = 0`

**`SMCO_BR(f, bounds_lower, bounds_upper, start_points = NULL, iter_boost = 1000, ...)`**
- Boosted refinement for maximum accuracy
- Sets `refine_search = TRUE`, `refine_ratio = 0.5`, `iter_boost = 1000` (customizable)

All wrapper functions pass additional arguments to `SMCO_multi()` via `...`.

### Full Control: `SMCO_multi()`

```r
SMCO_multi(f, bounds_lower, bounds_upper, start_points = NULL, opt_control = list(...))
```

**Arguments:**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `f` | Objective function to maximize | (required) |
| `bounds_lower` | Lower bounds vector | (required) |
| `bounds_upper` | Upper bounds vector | (required) |
| `start_points` | Custom starting points matrix | Auto-generated |
| `opt_control` | List of control parameters | See below |

**Control Parameters (`opt_control`):**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `n_starts` | Number of starting points | max(5, sqrt(dim)) |
| `iter_max` | Maximum iterations per start | 500 |
| `bounds_buffer` | Bounds extension factor | 0.05 |
| `buffer_rand` | Randomize bounds extension | FALSE |
| `tol_conv` | Convergence tolerance | 1e-8 |
| `refine_search` | Enable refinement phase | TRUE |
| `refine_ratio` | Fraction for refinement | 0.5 |
| `partial_option` | "center" or "forward" | "center" |
| `use_runmax` | Track running maximum | TRUE |
| `iter_boost` | Boosted starting index | 0 |
| `use_parallel` | Enable parallel execution | FALSE |
| `verbose` | Print progress | FALSE |
| `seed` | Random seed | 123 |

**Returns:** List containing:
- `best_result`: Best optimization result with `x_optimal`, `f_optimal`, `iterations`
- `all_results`: Results from all starting points
- `opt_control`: Control parameters used
- `summary`: Summary statistics

---

## Repository Structure

```
SMCO/
├── README.md                       # This file
├── NEWS.md                         # Version history and changes
├── LICENSE                         # MIT License
├── SMCO.R                          # Core SMCO algorithm (v1.1.0)
├── testfuncs.R                     # Test functions (deterministic + neural network)
├── benchmark.R                     # User-configurable benchmark script
├── benchmark_version_comparison.R  # v1.0.0 vs v1.1.0 performance comparison
└── v1.0.0/                         # Archived original version
```

---

## Benchmarking

Use `benchmark.R` to test SMCO on standard optimization benchmarks:

```r
# Edit the USER CONFIGURATION section in benchmark.R to set:
# - name_config: test function (e.g., "Rastrigin", "Ackley", "Rosenbrock")
# - dim_config: problem dimension
# - n_replications: number of Monte Carlo replications
# - SMCO_options: algorithm parameters

source("benchmark.R")
```

Available test functions include: Rastrigin, Ackley, Griewank, Rosenbrock, Michalewicz, DixonPrice, Zakharov, Qing, and neural network regression (ANN).

---

## Performance: v1.0.0 vs v1.1.0

Version 1.1.0 includes significant performance optimizations. Run `benchmark_version_comparison.R` to reproduce these results.

### Speedup Summary

| Dimension | SMCO | SMCO-R | SMCO-BR | Average |
|-----------|------|--------|---------|---------|
| 2-D | ~15% | ~20% | ~15% | ~17% |
| 10-D | ~35% | ~40% | ~35% | ~37% |
| 50-D | ~45% | ~50% | ~45% | ~47% |

**Overall: ~30-40% average speedup** (higher gains in high dimensions)

### Key Optimizations

1. **Pre-computed h_step matrix**: Avoids repeated division in main loop (O(iter_max) divisions saved)
2. **In-place vector modification**: Reduces memory allocation in `compute_partial_signs()`
3. **Metadata-based runmax tracking**: Single vector copy at end instead of per-iteration
4. **Pre-computed fixed bounds**: When `buffer_rand=FALSE`, bounds computed once
5. **Vectorized Sobol scaling**: Uses R's column-wise recycling via transpose
6. **Adaptive defaults**: `n_starts` defaults to `max(5, sqrt(dim))` for better scaling
7. **Input validation**: New `validate_smco_inputs()` for clearer error messages

### Solution Quality

Solution quality remains equivalent between versions. The optimizations focus on computational efficiency without changing the algorithm's behavior.

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Citation

If you use SMCO in your research, please cite:

> X. Chen, Z. Chen, W.Y. Gao, X. Yan, & G. Zhang, Optimization via the strategic law of large numbers, *Proc. Natl. Acad. Sci. U.S.A.* 123 (4) e2519845123, https://doi.org/10.1073/pnas.2519845123 (2026).

**BibTeX:**
```bibtex
@article{chen2026smco,
  author = {Chen, Xiaohong and Chen, Zengjing and Gao, Wayne Yuan and Yan, Xiaodong and Zhang, Guodong},
  title = {Optimization via the strategic law of large numbers},
  journal = {Proceedings of the National Academy of Sciences},
  volume = {123},
  number = {4},
  pages = {e2519845123},
  year = {2026},
  doi = {10.1073/pnas.2519845123}
}
```

---

## Contact

For questions, comments, and bug reports: waynegao@upenn.edu
