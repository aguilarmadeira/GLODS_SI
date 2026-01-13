# src — Core Algorithm Implementation

This directory contains the MATLAB implementation of the GLODS-SI algorithm.

## Main Files

| File | Description |
|------|-------------|
| `glods_si.m` | Main GLODS-SI algorithm entry point |
| `add.m` | Merging mechanism (add procedure) |
| `poll.m` | Polling step with positive basis |
| `search.m` | Search step (2^n-Centers strategy) |
| `normalize.m` | Coordinate normalization utilities |

## Algorithm Overview

GLODS-SI operates in a two-space formulation:

1. **Normalized space** [0,1]^n: All geometric operations (polling, distance computations, merging decisions)
2. **Original space** [l,u]: Objective function evaluations via inverse affine map

### Key Components

**Variable Normalization**
```matlab
% Forward: x -> y (original to normalized)
y = (x - lb) ./ (ub - lb);

% Inverse: y -> x (normalized to original)
x = lb + y .* (ub - lb);
```

**Merging Mechanism**

The `add` procedure maintains a list of sampled points with:
- Position (in normalized coordinates)
- Objective value (evaluated in original coordinates)
- Activity status (active/inactive)
- Step size and comparison radius

Points are merged when within each other's comparison radius, with dominance determined by objective values and sufficient decrease conditions.

## Usage

```matlab
% Add to path
addpath('src');

% Basic call
[x_opt, f_opt, output] = glods_si(f, lb, ub);

% With options
opts.alpha0 = 0.1;      % Initial step size
opts.maxevals = 10000;  % Max evaluations
opts.tol = 1e-8;        % Stopping tolerance

[x_opt, f_opt, output] = glods_si(f, lb, ub, opts);
```

## References

- Custódio, A.L., Madeira, J.F.A. (2015). GLODS: Global and Local Optimization using Direct Search. *J. Global Optim.* 62, 1–28.
