# scaling — Scaling Transformation Utilities

This directory contains utilities for applying deterministic scaling transformations to optimization problems.

## Files

| File | Description |
|------|-------------|
| `generate_scale_factors.m` | Generate per-variable scale factors for a given strategy |
| `generate_scaled_bounds.m` | Compute scaled bounds from original bounds and scale factors |
| `create_scaled_wrapper.m` | Create objective function wrapper that preserves landscape structure |

## Scaling Strategies

Seven deterministic scaling strategies are implemented:

| Strategy | Contrast (kappa) | Description |
|----------|------------------|-------------|
| `baseline` | 1 | No scaling (identity) |
| `progressive` | 10^6 | Cyclic geometric pattern |
| `extreme` | 10^8 | Binary partition (half small, half large) |
| `sobol_oscillatory` | 10^6 | Sobol-sequence driven continuous pattern |
| `sobol_digit_oscillatory` | 10^8 | Sobol digit-based binary assignment |
| `halton_oscillatory` | 10^6 | Halton-sequence driven pattern |
| `spatial_thermal` | ~9 x 10^4 | Multiphysics-inspired (spatial + thermal) |

## Usage

### Generate Scale Factors

```matlab
n = 10;  % Problem dimension

% Baseline (no scaling)
s = generate_scale_factors(n, 'baseline');

% Extreme scaling (kappa = 1e8)
s = generate_scale_factors(n, 'extreme');

% Progressive scaling (kappa = 1e6)
s = generate_scale_factors(n, 'progressive');
```

### Apply Scaling to Bounds

```matlab
lb_orig = -5 * ones(n, 1);
ub_orig = 5 * ones(n, 1);

s = generate_scale_factors(n, 'extreme');
[lb_scaled, ub_scaled] = generate_scaled_bounds(lb_orig, ub_orig, s);
```

### Create Problem Wrapper

```matlab
% Original problem
f_orig = @rosenbrock;
lb_orig = [-5; -5];
ub_orig = [5; 5];

% Scaling
s = generate_scale_factors(2, 'extreme');
[lb_s, ub_s] = generate_scaled_bounds(lb_orig, ub_orig, s);

% Create wrapper (preserves landscape structure)
f_wrapper = create_scaled_wrapper(f_orig, lb_orig, ub_orig, lb_s, ub_s);

% Now use f_wrapper with scaled bounds
x_scaled = rand(2,1) .* (ub_s - lb_s) + lb_s;
f_val = f_wrapper(x_scaled);  % Evaluates at corresponding original point
```

## Coordinate Mapping

The wrapper preserves problem structure through relative coordinate mapping:

```
t = (x_scaled - lb_scaled) ./ (ub_scaled - lb_scaled)   % Relative position [0,1]
x_orig = lb_orig + t .* (ub_orig - lb_orig)             % Map to original space
f = f_original(x_orig)                                   % Evaluate
```

This ensures:
- Global and local minima are preserved in relative coordinates
- Landscape structure (multimodality, basins) remains unchanged
- Only the geometric scaling of the search space is modified

## References

See Appendix B of the GLODS-SI paper for full mathematical details.
