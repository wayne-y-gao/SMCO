# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Strategic Monte Carlo Optimization (SMCO) - an R implementation of a novel optimization algorithm based on the paper "Optimization via Strategic Law of Large Numbers" by Chen, Chen, Gao, Yan, and Zhang.

**Citation:**
> Chen X, Chen X, Gao WY, Yan Z, Zhang M (2025). Optimization via strategic law of large numbers. *Proceedings of the National Academy of Sciences*, 122(4), e2419436122. https://doi.org/10.1073/pnas.2419436122

This is a pure R research codebase with no build system. All scripts are self-contained.

## Repository Structure

```
SMCO/
├── SMCO.R              # Core algorithm (v1.1.0)
├── testfuncs.R         # Test functions
├── benchmark.R         # Benchmarking script
├── NEWS.md             # Version history
├── README.md           # Project documentation
├── CLAUDE.md           # This file (project guidance)
├── LICENSE             # License file
└── v1.0.0/             # Archived original implementation
    ├── SMCO.R          # Core algorithm
    ├── RunComparison.R # Benchmarking framework
    ├── testfuncs.R     # Test functions
    ├── testfunc_NN1L.R # Neural network test functions
    ├── GD.R            # Gradient descent
    ├── SignGD.R        # Sign gradient descent
    ├── SPSA.R          # SPSA optimizer
    ├── LICENSE         # License file
    └── README.md       # Original readme
```

## Development Guidelines

- **`v1.0.0/` is a read-only archive** - DO NOT modify files in this folder
- **New development occurs at root level** - Create new files and versions here
- **Reference `v1.0.0/` for understanding** - Use the archived code as reference for algorithm logic

## Running the Code (Current v1.1.0)

### Using SMCO as a standalone optimizer
```r
source("SMCO.R")

result <- SMCO_multi(f, bounds_lower, bounds_upper, start_points,
                     opt_control = list(
                       n_starts = 100,        # number of starting points
                       iter_max = 200,        # max iterations per start
                       bounds_buffer = 0.05,  # bounds extension factor
                       refine_search = TRUE,  # enable refinement phase
                       refine_ratio = 0.5,    # fraction of iterations for refinement
                       partial_option = "center",  # "center" or "forward"
                       use_runmax = TRUE,     # track running maximum
                       use_parallel = FALSE
                     ))
```

### Running benchmarks
```r
source("benchmark.R")
# Runs SMCO on standard test functions and outputs results
```

## Running the Code (v1.0.0 Archive)

The following instructions apply to the archived v1.0.0 implementation.

### Using SMCO v1.0.0
```r
source("v1.0.0/SMCO.R")

result <- SMCO_multi(f, bounds_lower, bounds_upper, start_points,
                     opt_control = list(
                       n_starts = 100,
                       iter_max = 200,
                       bounds_buffer = 0.05,
                       refine_search = TRUE,
                       refine_ratio = 0.5,
                       partial_option = "center",
                       use_runmax = TRUE,
                       use_parallel = FALSE
                     ))
```

### Replicating paper experiments (v1.0.0)
1. Navigate to the `v1.0.0/` directory
2. Edit the configuration section at the end of `RunComparison.R`
3. Run: `source("v1.0.0/RunComparison.R")`

## Architecture

### Current v1.1.0 (SMCO.R at root)
Optimized implementation with performance improvements:
- `compute_partial_signs()` - computes discrete partial derivative signs
- `SMCO_single()` - base algorithm iterations
- `SMCO_single_refine()` - adds refinement search phase
- `SMCO_single_boost()` - adds boosted search with higher starting index
- `SMCO_multi()` - multi-start wrapper, returns best result

Key v1.1.0 optimizations:
- Removed `compiler::cmpfun()` dependency (modern R JIT handles this)
- Simplified Sobol sequence generation
- Improved default parameters

### v1.0.0 Archive
Original implementation with comparison framework including alternative optimizers (GD, SignGD, SPSA) and comprehensive benchmarking against packages like GenSA, DEoptim, GA, pso, etc.

## Key Dependencies

Required R packages:
- `qrng` - Sobol sequence generation for starting points

## Key Concepts

- **Bounds buffer**: SMCO extends the search domain slightly beyond specified bounds (controlled by `bounds_buffer`)
- **Refinement search**: Two-phase optimization - initial exploration then local refinement
- **Running maximum**: Tracks best solution found during gradient sign computations
- **Partial options**: "center" uses two-sided perturbations, "forward" uses one-sided

## Contact

Repository maintained by Wayne Yuan Gao (waynegao@upenn.edu)
