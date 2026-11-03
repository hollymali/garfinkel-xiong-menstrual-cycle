# Reference-amplitude model variant

`MODEL_reference_amplitude.py` is a side-by-side alternative to `MODEL.py`.
No existing Python file is changed. The endocrine model, hormone conditions,
six-second Euler solver, warm-up, and `run()` return format are identical to
the current model; only the eight thermal weights below are reweighted.

## Why this variant exists

The current shared thermal-clock drive produces daily oscillations that are
too large in both constant-hormone conditions, especially pregnancy/high E2.
The variant targets constant-condition amplitudes near one-half of cycling
while treating the established cycling bimodality as a hard shape constraint.

| Condition | Reference amplitude | Current model | Reference-amplitude variant |
|---|---:|---:|---:|
| Cycling E2 / non-pregnant | 2.449 °C | 2.024 °C | 2.012 °C |
| Constant low E2 / aged | 1.388 °C | 1.766 °C | 0.971 °C |
| Constant high E2 / pregnancy | 0.824 °C | 1.629 °C | 0.944 °C |

The references are the mean of the five individual 24-hour peak-to-trough
ranges. The model values are the range of the retained 120-hour output; its
constant-condition cycles repeat, so this is also their daily range.
The variant's constant-low/cycling and constant-high/cycling amplitude ratios
are 0.483 and 0.469, respectively.

## Reweighted parameters

| Parameter | `MODEL.py` | Variant |
|---|---:|---:|
| Conservation clock weight | 2.27 | 2.365 |
| N conservation weight | 0.60 | 0.979 |
| Progesterone conservation weight | 0.0274 | 0.0379 |
| Loss clock weight | 0.78 | 0.580 |
| E2 loss weight | 0.0000869 | 0.000566 |
| Basal heat rate | 49.44 | 48.11 |
| Conservation gain | 9.76 | 8.97 |
| Hormonal-loss gain | 15.29 | 24.98 |

This is not a uniform amplitude scaling. The modestly larger N, progesterone,
and E2 weights distinguish the two constant hormone states, while the loss
clock weight and the two downstream heat-balance gains set their smaller
clock-driven ranges. The cycling coefficients were constrained against the
entire original standardized cycling trace, not only its peak and trough. Its
exact-solver correlation with the original cycling waveform is 0.978, and the
two-lobe pattern remains visible in every modeled day.

The final 5% increase from the fitted N weight (0.932 to 0.979) is entirely
downstream of the endocrine cascade. It gives the first cycling lobe a small
additional lift without changing GnRH, LH, or E2 at any time point.

## Overlay validation

Using exactly the alignment rules in `plot_model.py`, the side-by-side variant
changes comparison RMSE as follows:

| Comparison | Current model | Variant |
|---|---:|---:|
| Cycling vs non-pregnant `mousebodytemp` | 0.694 °C | 0.754 °C |
| Constant low E2 vs aged 2011 (clock anchored) | 0.863 °C | 0.609 °C |
| Pregnancy vs pregnant-control `mousebodytemp` | 0.697 °C | 0.422 °C |

The constant conditions are necessarily smoother and more regular than the
digitized animals because the model gives them constant hormones plus one
idealized clock. This coefficient-only reweighting improves both constant
comparisons while preserving cycling bimodality. The cycling RMSE increases
slightly because the half-amplitude ratio and shape constraint were prioritized
over an unconstrained fit to that trace.

## Run without touching the current outputs

```powershell
python plot_reference_amplitude_model.py
```

This writes the side-by-side files prefixed with `reference_amplitude_` plus
the two standalone 300-DPI manuscript panels:

- `manuscript_cycling_e2_waveform_aligned.png`
- `manuscript_pregnancy_high_e2_waveform_aligned.png`

The dependency-free fit that generated the unrounded coefficients is
`fit_temperature_amplitudes.js`.
