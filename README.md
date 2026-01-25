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

# Run optimization
result <- SMCO_multi(f, bounds_lower, bounds_upper,
                     opt_control = list(
                       n_starts = 100,
                       iter_max = 200,
                       verbose = TRUE
                     ))

# Get results
print(result$best_result$x_optimal)  # Optimal point
print(result$best_result$f_optimal)  # Optimal value
```

---

## API Reference

### Main Function: `SMCO_multi()`

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
| `n_starts` | Number of starting points | 100 |
| `iter_max` | Maximum iterations per start | 200 |
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
├── README.md           # This file
├── NEWS.md             # Version history and changes
├── LICENSE             # MIT License
├── SMCO.R              # Core SMCO algorithm (v1.1.0)
├── testfuncs.R         # Test functions (deterministic + neural network)
├── benchmark.R         # User-configurable benchmark script
└── v1.0.0/             # Archived original version
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
