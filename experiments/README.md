# experiments — Reproducibility Scripts

This directory contains scripts to reproduce the numerical experiments from the GLODS-SI paper.

## Files

| File | Description |
|------|-------------|
| `run_experiments.m` | Main script to run all experiments |
| `run_single_strategy.m` | Run experiments for a single scaling strategy |
| `generate_figures.m` | Generate all figures from the paper |
| `compute_data_profiles.m` | Compute data/performance profiles |
| `compare_glods_vs_glods_si.m` | Head-to-head comparison |

## Prerequisites

1. **GLODS-SI code**: Ensure `src/` is in your MATLAB path
2. **Benchmark problems**: Download from [DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite)
3. **Original GLODS** (for comparison): Available at http://ferrari.dmat.fct.unl.pt/personal/alcustodio/glods

## Running Experiments

### Full Experiment Suite

```matlab
% Add paths
addpath('../src');
addpath('../scaling');

% Run all experiments (may take several hours)
run_experiments;
```

### Single Scaling Strategy

```matlab
% Run only baseline experiments
run_single_strategy('baseline');

% Run extreme scaling experiments
run_single_strategy('extreme');
```

### Generate Figures

```matlab
% After running experiments, generate paper figures
generate_figures;
% Figures saved to ../figures/
```

## Experimental Setup

As described in the paper:

- **Test set**: 63 bound-constrained problems (n = 2 to 10)
- **Scaling strategies**: 7 (baseline + 6 heterogeneous)
- **Stopping criteria**: alpha < 1e-8 OR 20,000 evaluations
- **Initialization**: 30 Sobol points per run
- **Performance measure**: Data profiles with tau = 1e-5

## Output Files

Results are saved to:
- `results/` — Raw experimental data (MAT files)
- `../figures/` — Generated figures (PDF/PNG)

## Expected Results

The experiments should reproduce Table 3 from the paper:

| Strategy | Contrast | GLODS (%) | GLODS-SI (%) |
|----------|----------|-----------|--------------|
| baseline | 1:1 | 84.1 | 87.3 |
| progressive | 10^6 | 76.2 | 88.9 |
| extreme | 10^8 | 65.1 | 93.7 |
| sobol_oscillatory | 10^6 | 71.4 | 90.5 |
| sobol_digit_oscillatory | 10^8 | 61.9 | 90.5 |
| halton_oscillatory | 10^6 | 65.1 | 92.1 |
| spatial_thermal | ~9e4 | 68.3 | 93.7 |

## Notes

- Experiments are deterministic (Sobol initialization)
- Total runtime: approximately 4-8 hours on a modern desktop
- Memory usage: < 4 GB RAM
