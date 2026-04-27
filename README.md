# GLODS-SI

**Scale-Invariant Global-Local Direct Search for Engineering Design Optimization**

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-pending-orange.svg)](https://doi.org/)

GLODS-SI is a derivative-free direct-search method for bound-constrained
global optimization problems whose decision variables span widely
different physical scales. It extends the original GLODS framework
(Custódio & Madeira, 2015) by performing all geometric operations —
polling, distance computation, point merging — in normalized coordinates
`y ∈ [0, 1]^n`, while objective evaluations remain in the original
variable space `x ∈ [ℓ, u]`. This *two-space formulation* makes the
method's behaviour invariant to the units chosen for each variable and
robust to heterogeneous variable scales.

This repository contains the reference MATLAB implementation that
accompanies the paper:

> J. F. A. Madeira,
> *GLODS-SI: Scale-Invariant Global-Local Direct Search for Engineering
> Design Optimization*,
> Journal of Computational Design and Engineering, 2026.
> Manuscript ID **JCDE-2026-065**.

The numerical results that produced Tables 3 and Figures 3–6 of the
paper are provided as final artifacts in [`results/`](results/).

---

## Repository contents

```
GLODS_SI/
├── README.md                  (this file)
├── LICENSE                    (GNU LGPL v3)
├── CITATION.cff               (citation metadata)
├── glods_si.m                 (main solver)
├── parameters_glods_si.m      (default parameters; matches paper Section 4)
├── driver_glods_si.m          (example driver)
├── aluffi_pentini_2D.m        (illustrative test problem, 2D)
└── results/                   (numerical results that produced the paper's
                               tables and figures; see results/README.md)
```

---

## Installation

Requirements: MATLAB R2018b or later (no toolboxes required).

1. Clone or download this repository.
2. Add the repository folder to your MATLAB path:
   ```matlab
   addpath('path/to/GLODS_SI');
   ```

That's it.

---

## Quick start

Run the bundled demo on the 2D Aluffi-Pentini problem:

```matlab
cd GLODS_SI
driver_glods_si
```

The driver runs GLODS-SI with the default parameters listed in
`parameters_glods_si.m` and reports:

- objective value reached;
- best point in work-space coordinates;
- total number of function evaluations;
- a convergence profile (best-so-far value vs. function evaluations).

Aluffi-Pentini is a 2D problem with a single global minimum at
`x* = (-1.0465, 0)` and `f* = -0.3523`, plus two local minima.
GLODS-SI typically reaches the global minimum within a few hundred
function evaluations.

---

## Using GLODS-SI on your own problem

GLODS-SI expects:

- a function handle `@my_objective` that accepts a column vector
  `x ∈ R^n` and returns a scalar `f`;
- bound vectors `lb`, `ub` of size `n × 1`.

The solver is then called as:

```matlab
[profile, Plist, flist, alfa, radius, fevals] = ...
    glods_si(@my_objective, [], [], lb, ub);
```

Outputs:

| Output     | Meaning                                                    |
|------------|------------------------------------------------------------|
| `profile`  | Best-so-far objective value vs. function evaluations       |
| `Plist`    | Final list of points returned by the solver                |
| `flist`    | Objective values at the points in `Plist`                  |
| `alfa`     | Final step sizes (one per point in `Plist`)                |
| `radius`   | Final comparison radii                                     |
| `fevals`   | Total number of function evaluations performed             |

For benchmarking against the test suite used in the paper, see the
companion repository
[DFO_Benchmark_Suite](https://github.com/aguilarmadeira/DFO_Benchmark_Suite),
which provides 504 self-contained wrappers (63 instances × 8 scaling
strategies) ready to be passed as the first argument to `glods_si`.

---

## Default parameters

The default values in `parameters_glods_si.m` match the convention used
in the experimental section of the paper (Section 4):

| Parameter           | Value     | Meaning                                 |
|---------------------|-----------|-----------------------------------------|
| `alfa_ini`          | 0.1       | Initial step size in normalized space   |
| `radius_ini`        | 0.2       | Initial comparison radius (= 2·`alfa_ini`)|
| `tol_stop`          | 1e-5      | Stopping tolerance                      |
| `max_fevals`        | 20000     | Maximum number of function evaluations  |
| `nPini`             | 30        | Number of initial sample points         |
| `list`              | 6         | Initialization: Sobol sequences         |
| `suf_decrease`      | 0         | Integer-lattice sufficient decrease     |
| `cache`             | 1         | Cache previously evaluated points       |

Two changes relative to the original GLODS defaults are worth noting:

- `alfa_ini` and `radius_ini` are *dimensionless* (fractions of each
  variable range), since geometric operations occur in `[0, 1]^n`.
- `radius_ini = 2 · alfa_ini = 0.2` is the convention used to generate
  Table 3 of the paper. It satisfies `r_0 ≥ d_max · alfa_ini` (with
  `d_max = 1` for the positive basis `[I -I]`), preventing nearby
  initial points from being merged immediately and favouring
  exploration of the global structure.

---

## Numerical results

The [`results/`](results/) folder contains the numerical results that
produced the tables and figures of the accompanying paper, organized
into three sub-folders:

### Main comparison (3-way)

[`results/GLODS_vs_NOMAD_vs_GLODSSI/`](results/GLODS_vs_NOMAD_vs_GLODSSI/)
contains the joint 3-algorithm data profiles, including the figures
used in the manuscript:

- **Figure 3** — baseline scaling (κ = 1)
- **Figure 4** — extreme scaling (κ = 10⁸)
- **Figure 5** — Halton oscillatory scaling (κ = 10⁶)
- **Figure 6** — spatial–thermal scaling (κ ≈ 9 × 10⁴)

Profiles for the four scaling strategies not shown in the paper are
also provided, together with the auto-generated LaTeX summary table.

### Pairwise comparisons

- [`results/GLODS_vs_GLODSSI/`](results/GLODS_vs_GLODSSI/) — source of
  the GLODS and GLODS-SI columns of Table 3.
- [`results/NOMAD_vs_GLODSSI/`](results/NOMAD_vs_GLODSSI/) — source of
  the NOMAD column of Table 3.

Each sub-folder contains, for the eight scaling strategies considered
in the paper: PDF data profiles, ASCII summary tables, and the
auto-generated LaTeX success-rate tables. See
[`results/README.md`](results/README.md) for the full description.

These files are final artifacts and can be inspected directly without
rerunning the experiments.

---

## License

This software is distributed under the **GNU Lesser General Public
License version 3** (LGPL-3.0-or-later). This license is inherited
from the original GLODS framework (Custódio & Madeira, 2015), on which
GLODS-SI is based.

See [`LICENSE`](LICENSE) for the full license text.

---

## Citation

If you use GLODS-SI in academic work, please cite the accompanying paper:

```bibtex
@article{Madeira2026GLODSSI,
  author  = {Madeira, J. F. A.},
  title   = {{GLODS-SI}: Scale-Invariant {Global--Local} Direct Search
             for Engineering Design Optimization},
  journal = {Journal of Computational Design and Engineering},
  year    = {2026},
  note    = {Manuscript ID JCDE-2026-065}
}
```

A `CITATION.cff` file is also provided for tools that consume that
metadata format (GitHub, Zenodo, etc.).

### Underlying framework

> Custódio, A. L., Madeira, J. F. A. (2015).
> *GLODS: Global and Local Optimization using Direct Search.*
> Journal of Global Optimization, 62, 1–28.
> [doi:10.1007/s10898-014-0224-9](https://doi.org/10.1007/s10898-014-0224-9)

---

## Acknowledgments

This work was supported by *Fundação para a Ciência e a Tecnologia*
(FCT) through LAETA (project
[UID/50022/2025](https://doi.org/10.54499/UID/50022/2025)).

---

## Contact

**José F. Aguilar Madeira**
IDMEC, Instituto Superior Técnico, Universidade de Lisboa
ISEL, Instituto Politécnico de Lisboa
Email: aguilarmadeira@tecnico.ulisboa.pt
ORCID: [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)
