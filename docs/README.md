# docs — Documentation

Additional documentation for GLODS-SI.

## Contents

| Document | Description |
|----------|-------------|
| `algorithm_details.md` | Detailed algorithm description |
| `convergence_analysis.md` | Summary of convergence guarantees |
| `api_reference.md` | Function API reference |

## Quick Reference

### Main Function

```matlab
[x_opt, f_opt, output] = glods_si(f, lb, ub, options)
```

**Inputs:**
- `f` — Objective function handle `f(x)` returning scalar
- `lb` — Lower bounds (n x 1 vector)
- `ub` — Upper bounds (n x 1 vector)
- `options` — (Optional) Structure with algorithm parameters

**Outputs:**
- `x_opt` — Best solution found
- `f_opt` — Objective value at x_opt
- `output` — Structure with optimization history

### Options Structure

```matlab
options.alpha0    = 0.1;    % Initial step size (normalized space)
options.r0        = 0.2;    % Initial comparison radius
options.beta      = 0.5;    % Step contraction factor
options.gamma     = 2.0;    % Step expansion factor
options.tol       = 1e-8;   % Stopping tolerance
options.maxevals  = 20000;  % Maximum evaluations
options.ninit     = 30;     % Number of initial points
options.verbose   = false;  % Display progress
```

## Related Resources

- **Paper**: Madeira (2026), "GLODS-SI: Global and Local Optimization using Direct Search — A Scale-Invariant Approach"
- **Original GLODS**: [http://ferrari.dmat.fct.unl.pt/personal/alcustodio/glods](http://ferrari.dmat.fct.unl.pt/personal/alcustodio/glods)
- **Benchmark Suite**: [https://github.com/aguilarmadeira/DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite)
