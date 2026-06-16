# Rasgon et al. (2003) Menstrual-Cycle Model

A reference reimplementation of the HPG-axis (hypothalamic–pituitary–gonadal) menstrual-cycle
model from:

> Rasgon NL, Pumphrey L, Prolo P, Elman S, Negrao AB, Licinio J, Garfinkel A.
> **"Emergent Oscillations in Mathematical Model of the Human Menstrual Cycle."**
> *CNS Spectrums.* 2003;8(11):805–814. https://doi.org/10.1017/S1092852900019246

The model tracks four hormones — estradiol (E2), follicle-stimulating hormone (FSH),
luteinizing hormone (LH), and progesterone (P4) — across two phase-specific ODE systems
(Early Follicular Phase and Mid-Luteal Phase), each driven by an explicit GnRH pulse
generator. It reproduces the paper's central claim: the cycle's irregular oscillations
emerge from steroid-feedback dynamics rather than from an external clock.

## Contents

| File | Description |
|------|-------------|
| `menstrual_model.py` | The model: EFP/MLP ODE right-hand sides, Hill feedback functions, GnRH pulse generator, and integrators (adaptive `solve_ivp` RK45 plus a fixed-step classic RK4 matching the original Berkeley Madonna scheme). |
| `plot_model.py` | Time-series diagnostics driver. Imports the model and renders 5×2 panels (GnRH/FSH/LH/E2/P4, full run + zoom). |
| `fig_efp_faithful.png` | Output of the faithful Early Follicular Phase transcription. |
| `emergent-oscillations-in-mathematical-model-of-the-human-menstrual-cycle.pdf` | The source paper. |

## Usage

Requires `numpy`, `scipy`, and `matplotlib`.

```bash
# Smoke test (prints per-hormone min/max/range for both phases)
python menstrual_model.py

# Generate the time-series figures
python plot_model.py
```

## Modernization choices

Two changes from the original 2003 Berkeley Madonna implementation:

1. **Integration** uses `scipy.integrate.solve_ivp` (adaptive RK45). A fixed-step RK4
   path (`run_rk4`) reproduces the original DT = 0.25 min scheme for cross-checking.
2. **GnRH pulses** are reimplemented as a step-size-independent forcing function: an
   internal phase variable advances at rate `1/Period(state)` and a narrow Gaussian
   bump of peak height `Amplitude(state)` fires at each integer phase — replacing
   Berkeley Madonna's `PULSE()` builtin, whose spike height was tied to the step size.

All equations and parameters were transcribed from Appendices A–D of the paper. The
model is qualitative: hormone magnitudes are not meaningful in absolute terms — only the
dynamics and shapes are.

## License / attribution

Code is a derivative reimplementation for study and verification. The bundled paper PDF
is © the original authors / *CNS Spectrums* (Cambridge University Press) and is included
for reference only.
