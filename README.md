# GLODS-SI

**Global and Local Optimization using Direct Search — A Scale-Invariant Approach**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)]()

> **Note**: Implementation files will be added upon paper acceptance. The repository structure and documentation are provided for reference.

---

GLODS-SI is a scale-invariant reformulation of the [GLODS](https://link.springer.com/article/10.1007/s10898-014-0224-9) framework for derivative-free, bound-constrained global optimization.

## Key Innovation

GLODS-SI addresses the sensitivity of GLODS to heterogeneous variable scales through a **two-space formulation**:

- **Geometric operations** (polling, distances, merging) are performed in **normalized coordinates** [0,1]^n
- **Objective evaluation** is carried out in the **original variables** via an inverse affine map

This ensures that step sizes and merging decisions have a consistent, dimensionless interpretation independent of physical units and variable ranges.

## Features

- Scale-invariant geometry for robust performance under heterogeneous variable scales
- Multistart strategy with intelligent merging to avoid redundant convergence
- Preserves all convergence guarantees of the original GLODS framework
- Validated on 63 benchmark problems under 7 scaling regimes (contrast up to 10^8)

## Repository Structure

```
GLODS_SI/
├── README.md               # This file
├── LICENSE.md              # MIT License
├── CITATION.cff            # Citation metadata
├── src/                    # Core algorithm implementation
│   ├── glods_si.m          # Main GLODS-SI algorithm
│   ├── add.m               # Merging mechanism
│   └── ...
├── scaling/                # Scaling transformation utilities
│   ├── generate_scale_factors.m
│   ├── generate_scaled_bounds.m
│   └── create_scaled_wrapper.m
├── experiments/            # Scripts to reproduce paper results
│   ├── run_experiments.m
│   └── generate_figures.m
├── figures/                # Figures from the paper
└── docs/                   # Additional documentation
```

## Quick Start

### Basic Usage

```matlab
% Add GLODS-SI to path
addpath('src');

% Define your problem
f = @(x) sum(x.^2);    % Objective function
lb = [-5; -5];         % Lower bounds
ub = [5; 5];           % Upper bounds

% Run GLODS-SI
[x_opt, f_opt, output] = glods_si(f, lb, ub);

fprintf('Optimal value: %.6f\n', f_opt);
```

### With Heterogeneous Scaling

```matlab
% Add scaling utilities
addpath('scaling');

% Original problem
f_orig = @rosenbrock;
lb_orig = [-5; -5];
ub_orig = [5; 5];

% Generate scale factors (e.g., extreme scaling with kappa = 1e8)
n = 2;
s = generate_scale_factors(n, 'extreme');

% Create scaled bounds
[lb_scaled, ub_scaled] = generate_scaled_bounds(lb_orig, ub_orig, s);

% Create wrapper that preserves landscape structure
f_scaled = create_scaled_wrapper(f_orig, lb_orig, ub_orig, lb_scaled, ub_scaled);

% Run GLODS-SI (automatically handles normalized geometry)
[x_opt, f_opt] = glods_si(f_scaled, lb_scaled, ub_scaled);
```

## Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha0`  | 0.1     | Initial step size in normalized space |
| `r0`      | 0.2     | Initial comparison radius (typically 2 * alpha0) |
| `beta`    | 0.5     | Step-size contraction factor |
| `gamma`   | 2.0     | Step-size expansion factor |
| `tol`     | 1e-8    | Step-size stopping tolerance |
| `maxevals`| 20000   | Maximum function evaluations |

## Benchmark Suite

The benchmark problems and pre-computed samples used in the paper are available in a separate repository:

👉 [DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite)

This includes:
- 63 bound-constrained test problems
- 7 scaling strategies (baseline to kappa = 10^8)
- 441 self-contained problem instances
- Pre-computed initial samples in HDF5/MAT format

## Numerical Results

GLODS-SI maintains high success rates (87–94%) across all scaling regimes, while GLODS degrades significantly under scale heterogeneity:

| Scaling Strategy | Contrast (kappa) | GLODS | GLODS-SI |
|------------------|------------------|-------|----------|
| Baseline         | 1                | 84.1% | 87.3%    |
| Progressive      | 10^6             | 76.2% | 88.9%    |
| Extreme          | 10^8             | 65.1% | 93.7%    |
| Halton oscillatory | 10^6           | 65.1% | 92.1%    |

## Citation

If you use GLODS-SI in your research, please cite:

```bibtex
@article{madeira2026glods_si,
  title={{GLODS-SI}: Global and Local Optimization using Direct Search 
         -- A Scale-Invariant Approach},
  author={Madeira, Jos{\'e} F. A.},
  journal={Journal of Global Optimization},
  year={2026},
  note={Submitted}
}
```

For the original GLODS framework:

```bibtex
@article{custodio2015glods,
  title={{GLODS}: Global and Local Optimization using Direct Search},
  author={Cust{\'o}dio, A. L. and Madeira, J. F. A.},
  journal={Journal of Global Optimization},
  volume={62},
  pages={1--28},
  year={2015},
  doi={10.1007/s10898-014-0224-9}
}
```

## Related Work

- **GLODS** (original): http://ferrari.dmat.fct.unl.pt/personal/alcustodio/glods
- **DFO_Benchmark_Suite**: https://github.com/aguilarmadeira/DFO_Benchmark_Suite

## License

This project is licensed under the MIT License — see [LICENSE.md](LICENSE.md) for details.

## Author

**José F. A. Madeira**  
IDMEC, Instituto Superior Técnico, Universidade de Lisboa  
ISEL, Instituto Politécnico de Lisboa  
Email: aguilarmadeira@tecnico.ulisboa.pt

## Acknowledgments

This work was supported by Fundação para a Ciência e a Tecnologia (FCT) through LAETA (project [DOI: 10.54499/UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).
