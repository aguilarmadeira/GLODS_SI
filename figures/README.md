# figures — Paper Figures

This directory contains the figures from the GLODS-SI paper.

## Figure List

| Figure | File | Description |
|--------|------|-------------|
| Fig. 1 | `scale_bias_illustration.pdf` | Scale bias in GLODS merging mechanism |
| Fig. 2 | `two_space_formulation.pdf` | GLODS vs GLODS-SI: two-space mechanism |
| Fig. 3 | `baseline_profile.pdf` | Data profile: baseline scaling (kappa = 1) |
| Fig. 4 | `progressive_profile.pdf` | Data profile: progressive scaling (kappa = 10^6) |
| Fig. 5 | `extreme_profile.pdf` | Data profile: extreme scaling (kappa = 10^8) |
| Fig. 6 | `sobol_oscillatory_profile.pdf` | Data profile: Sobol oscillatory (kappa = 10^6) |
| Fig. 7 | `sobol_digit_profile.pdf` | Data profile: Sobol digit-driven (kappa = 10^8) |
| Fig. 8 | `halton_profile.pdf` | Data profile: Halton oscillatory (kappa = 10^6) |
| Fig. 9 | `spatial_thermal_profile.pdf` | Data profile: spatial-thermal (kappa ~ 9e4) |

## Regenerating Figures

To regenerate these figures from experimental data:

```matlab
cd('../experiments');
generate_figures;
```

Figures will be saved in both PDF and PNG formats.

## Figure Format

- Vector format: PDF (for paper submission)
- Raster format: PNG at 300 DPI (for presentations)
- Color scheme: Blue (GLODS), Red dashed (GLODS-SI)

## Usage in Publications

These figures are provided under the MIT license. When using in publications, please cite the GLODS-SI paper.
