# Meeting 7 intermediate checkpoint — 2026-09-06

This checkpoint preserves the exact model version that generated
`temperature_timeseries.png` before revising the interpretation of the SCN,
progesterone, and activity terms.

## Implemented thermal constants

- Resting heat: `0.32 kcal/h`
- Whole-body conductance: `0.02425 kcal/(h·°C)` (provisional mass scaling for
  an approximately 25-g female; update when the matched mass is available)
- Body heat capacity: `0.018379 kcal/°C`
- Ambient temperature: `22 °C`
- Maximum E2-associated thermogenic fraction: `0.11`
- E2 thermogenesis Hill function: `K = 9.0 pg/mL`, `n = 1.41`
- E2 heat-conservation Hill function: `K = 7.49 pg/mL`, `n = 1.83`
- Maximum E2-associated whole-body conductance reduction: `0.0275`
- SCN temperature-equivalent scale: `0.3 °C`
- Activity temperature-equivalent scale: `1.9 °C`
- Progesterone temperature-equivalent scale: `0.3 °C`
- Pregnant E2: `50.88 pg/mL`
- Pregnant SCN multiplier: `0.51`
- Pregnant activity multiplier: `0.7322`
- Pregnant progesterone multiplier: `1.0`

## Reproduced outputs

| Condition | Mean (°C) | Minimum (°C) | Maximum (°C) | Range (°C) |
|---|---:|---:|---:|---:|
| Cycling E2 / nonpregnant | 36.92916 | 36.11504 | 38.53490 | 2.41987 |
| Constant-high E2 only | 37.24184 | 36.65352 | 38.74444 | 2.09093 |
| Pregnancy comparison | 37.48759 | 37.10031 | 38.58707 | 1.48677 |

The pregnancy comparison raises the trough by `0.98527 °C` and the peak by
`0.05217 °C` relative to the cycling-E2 trace.

## Questions deliberately left open for the next revision

1. The factors `G × 0.3 °C` and `G × 1.9 °C` convert reported or assigned
   temperature-equivalent effects into heat-flow terms using `Q = GΔT`. This
   is dimensionally consistent, but it is a phenomenological conversion—not
   evidence that SCN, progesterone, or activity act biologically through
   conductance. Consider expressing these as shifts in a defended-temperature
   target instead.
2. The pregnancy SCN multiplier of `0.51` was introduced by applying a reported
   49% reduction in inducible BAT capacity to the SCN-associated thermal term.
   The current SCN term represents a broader circadian temperature influence,
   so that BAT result does not directly justify halving the full SCN term. The
   user has flagged that the next revision should remove this multiplication
   (`m_SCN = 1` during pregnancy) and then recalibrate any remaining fitted
   pregnancy parameter transparently.
3. This checkpoint intentionally retains `m_SCN = 0.51` so the saved code and
   saved figure remain mutually reproducible.

## Files in this checkpoint

- `MODEL.py`: implemented endocrine and temperature equations
- `PLOT.py`: time-series plotting code
- `temperature_timeseries.png`: generated 72-hour result
- `README.md`: literature rationale and coefficient audit

