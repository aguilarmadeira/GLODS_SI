# Numerical results

This folder contains the numerical results that accompany the paper:

> J. F. A. Madeira,
> *GLODS-SI: Scale-Invariant Global-Local Direct Search for Engineering
> Design Optimization*,
> Journal of Computational Design and Engineering, 2026.
> Manuscript ID JCDE-2026-065.

The benchmark suite is the 63 bound-constrained test problems described
in Appendix A of the paper, evaluated under the eight scaling strategies
of Appendix B. The companion repository
[DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite)
provides the corresponding 504 self-contained MATLAB wrappers
(63 instances × 8 strategies).

---

## Folder structure

Three head-to-head and joint comparisons are provided. Each folder
contains, for the eight scaling strategies:

```
results/
├── GLODS_vs_GLODSSI/                 (GLODS    vs GLODS-SI)
├── NOMAD_vs_GLODSSI/                 (NOMAD    vs GLODS-SI)
└── GLODS_vs_NOMAD_vs_GLODSSI/        (GLODS, NOMAD, GLODS-SI jointly)
```

Each subfolder contains:

| File pattern                       | Count | Description                            |
|------------------------------------|------:|----------------------------------------|
| `Profile_<comparison>_<strategy>_BEST.pdf` |   8 | Data profile (PDF), one per strategy   |
| `Summary_<comparison>_<strategy>_BEST.txt` |   8 | Per-strategy success rates at selected evaluation budgets |
| `latex_table_<comparison>.tex`            |   1 | LaTeX tables aggregating all strategies |

---

## Mapping to the paper

| Paper element | Source files                                                                  |
|---------------|--------------------------------------------------------------------------------|
| Figure 3 (baseline)        | `GLODS_vs_NOMAD_vs_GLODSSI/Profile_3way_baseline_BEST.pdf`        |
| Figure 4 (extreme)         | `GLODS_vs_NOMAD_vs_GLODSSI/Profile_3way_extreme_BEST.pdf`         |
| Figure 5 (Halton oscillatory) | `GLODS_vs_NOMAD_vs_GLODSSI/Profile_3way_halton_oscillatory_BEST.pdf` |
| Figure 6 (spatial-thermal) | `GLODS_vs_NOMAD_vs_GLODSSI/Profile_3way_spatial_thermal_BEST.pdf` |
| Table 3, columns GLODS / GLODS-SI | `GLODS_vs_GLODSSI/latex_table_GLODS_vs_GLODSSI.tex`        |
| Table 3, column NOMAD             | `NOMAD_vs_GLODSSI/latex_table_NOMAD_vs_GLODSSI.tex`        |

Profiles for the four scaling strategies not shown in the paper
(`moderate`, `progressive`, `sobol_oscillatory`, `sobol_digit_oscillatory`)
are also provided for completeness in the same `GLODS_vs_NOMAD_vs_GLODSSI/`
folder.

---

## Experimental setup

All runs share the following configuration:

| Setting                  | Value                                            |
|--------------------------|--------------------------------------------------|
| Test problems            | 63 (Appendix A of the paper)                     |
| Scaling strategies       | 8 (Appendix B of the paper)                      |
| Evaluation budget        | 20,000 function evaluations                      |
| Initial sample           | 30 Sobol points                                  |
| Tolerance                | τ = 10⁻⁵                                         |
| GLODS / GLODS-SI         | `parameters_glods_si.m` (in repository root)     |
| NOMAD                    | NOMAD 4 default settings                         |

Runs are deterministic: the Sobol initialization and the search/poll
mechanisms used here use no random number generator. Each problem was
executed once per algorithm.

---

## Reading the summary files

Each `Summary_*.txt` file follows the same layout:

```
GLODS vs GLODS-SI | baseline | BEST | tau=1e-05 | so=5 ss=0
Valid: 63/63
Criterion: H<=f_min+tau*(f_0-f_min) (joint min across all algs)

Evals    | GLODS (original)     | GLODS-SI (scale-invariant) | SI-GLODS
-------------------------------------------------------------------------
100      |   17/63( 27.0%) |   17/63( 27.0%) |     +0
500      |   35/63( 55.6%) |   36/63( 57.1%) |     +1
...
final    |   55/63( 87.3%) |   57/63( 90.5%) |     +2
```

The `final` line corresponds to the success rate at the full evaluation
budget, which is the value reported in Table 3 of the paper.

---

## Reproducibility

To re-run any experiment, use the source code in the repository root
together with the test wrappers from
[DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite).
Each scaling strategy folder in that suite (e.g.,
`prob_var_continuous/scaled/baseline_1/`) contains 63 self-contained
wrappers that can be passed directly to `glods_si`:

```matlab
addpath('prob_var_continuous/scaled/baseline_1');
addpath('GLODS_SI');

% Example: run GLODS-SI on the Ackley 10D wrapper, baseline scaling
[lb, ub] = ackley_10D(10);
[profile, P, f, alfa, r, fevals] = glods_si(@ackley_10D, [], [], lb, ub);
```

The default parameters in `parameters_glods_si.m` (Sobol initialization,
30 initial points, α₀ = 0.1, r₀ = 0.2, τ = 10⁻⁵, 20,000 evaluations)
correspond to the configuration used to produce the results in this
folder.
