# Extracted Model — Rasgon et al. (2003)

**Source:** Rasgon NL, Pumphrey L, Prolo P, Elman S, Negrao AB, Licinio J, Garfinkel A.
"Emergent Oscillations in Mathematical Model of the Human Menstrual Cycle."
*CNS Spectrums.* 2003;8(11):805–814. https://doi.org/10.1017/S1092852900019246

This file captures everything quantitative in the paper: the two phase-specific ODE
systems (Early Follicular Phase, Mid-Luteal Phase), their parameters, the Hill/sigmoid
feedback functions, and the GnRH pulse generator. Transcribed from Appendices A–D by
high-resolution rendering of the original scan (the plain-text PDF layer is garbled).

State variables (4 hormones), with the paper's units:
- `E2`  — estradiol (estrogen), pmol/L
- `FSH` — follicle-stimulating hormone, IU/L
- `LH`  — luteinizing hormone, IU/L
- `P4`  — progesterone, nmol/L
- `GnRH`— gonadotropin-releasing hormone, driven by an explicit pulse generator (input, not an ODE)

Integration in the original: 4th-order Runge–Kutta, fixed step **0.25 min**, run for
**8,000 min (5.5 days)**.

---

## Appendix A — EARLY FOLLICULAR PHASE (EFP)

### Differential equations
```
dE2/dt  = .01 * (FSH->E2 positive feedback) - (kE2  * E2)
dFSH/dt = 2e-7 * [ (E2->FSH negative feedback) + (GnRH->FSH positive feedback) ] - (kFSH * FSH)
dLH/dt  = (.1 * Multiplier * GnRH) - (kLH * LH)
dP4/dt  = (.045 * LH) - (kP4 * P4)
```

### Clearance constants
```
kE2  = .011
kFSH = .009
kLH  = .002
kP4  = 9.2e-5
```

### State-dependent control variables (piecewise on E2)
```
Multiplier = if E2 < 128 -> 2
             elif E2 > 800 -> .5
             else -> -.002 * E2 + 2.26

Period     = if E2 < 250 -> 60
             elif E2 > 400 -> 100
             else -> .267 * E2 - 6.75          # minutes between GnRH pulses

Amplitude  = if E2 < 250 -> 5
             elif E2 > 400 -> 3
             else -> -.013 * E2 + 8.2          # GnRH pulse amplitude
```

### Feedback (Hill) functions — EFP  *(all confirmed against Appendix C plots)*
```
FSH->E2 positive feedback  = 150 * (.15*FSH)^7  / (1 + (.15*FSH)^7)  + 250   # rises 250 -> 400, mid FSH~7
GnRH->FSH positive feedback = 5  * (.2*GnRH)^15 / (1 + (.2*GnRH)^15) + 5     # rises 5 -> 10, mid GnRH=5
E2->FSH negative feedback  = 9 / (1 + (.0025*E2)^15) + 1                     # falls 10 -> 1, mid E2=400  [RECONSTRUCTED]
```

> **NOTE 1 — RESOLVED via Appendix C middle plot.** The printed
> `E2->FSH negative feedback = (1 + (.0025*E2)^15)/(5*(.2*GnRH)^15)` is corrupted: the
> `5*(.2*GnRH)^15` bled in from the `GnRH->FSH` fraction below, and the numerator/offset
> were lost. The Appendix C plot (x: E2 250–600, y: 0–10) shows a **downward** Hill: ≈10
> for E2<350, ≈1 for E2>500, midpoint ≈E2=400 — i.e. `9/(1+(.0025*E2)^15) + 1`
> (k=.0025 ⇒ midpoint E2=400, matching the print; range 1–10 read from the plot).
> The GnRH term's printed denominator exponent `^3` is a typo for `^15` (the plot confirms
> matching exponents, half-rise at GnRH=5).

---

## Appendix A — MID-LUTEAL PHASE (MLP)

### Differential equations
```
dE2/dt  = .005  * [ (FSH->E2 positive feedback) - (kE2 * E2) ]
dFSH/dt = .05   * [ (5 * GnRH) + (E2->FSH negative feedback) - (kFSH * FSH) ]
dLH/dt  = .0014 * [ (20 * GnRH) + (E2->LH negative feedback) + (P4->LH negative feedback) - (kLH * LH) ]
dP4/dt  = .08   * [ (3 * LH) - (kP4 * P4) ]
```

### GnRH pulse generator (MLP "self-priming": each pulse triggers two more)
```
GnRH = PULSE(amplitude, 1, period)
     + PULSE(amplitude, 3, period + 20)
     + PULSE(amplitude, 2, period + 30)
```
(One primary pulse plus two subsequent pulses — the GnRH self-priming effect.)

> **NOTE 4 — PULSE semantics (Berkeley Madonna).** The `PULSE(...)` syntax + fixed-step
> RK4 at DT=0.25 min identify the original tool as **Berkeley Madonna**, whose builtin is
> `PULSE(volume, firstpulse, interval)`: it returns `volume/DT` during the single timestep
> at t = firstpulse, firstpulse+interval, firstpulse+2*interval, … and 0 otherwise — so
> each pulse is an **impulse whose time-integral (area) = volume = amplitude**. Mapping the
> MLP generator:
> ```
> train 1: amplitude, first @ t=1, every  period      min
> train 2: amplitude, first @ t=3, every (period+20)  min
> train 3: amplitude, first @ t=2, every (period+30)  min
> ```
> The three slightly different intervals beat against each other, producing the drifting
> doublet/triplet bursts the paper describes (FSH/LH going "in and out of synchrony"). The
> EFP uses a single such train (amplitude, period). **Modernization choice:** since
> `volume/DT` is step-size dependent, reimplement GnRH as an explicit forcing function —
> a sum of narrow normalized kernels (area = amplitude) at the pulse times, or discrete
> impulse events between solver segments — rather than copying the DT-scaled spike.

### Clearance constants
```
kE2  = .9
kFSH = .4
kLH  = 1.9
kP4  = .0835
```

### State-dependent control variables (piecewise on E2+P4)
```
Period    = if (E2+P4) < 282.55 -> 120
            elif (E2+P4) > 293.3 -> 240
            else -> 11.111 * (E2+P4) - 3018.88889

Amplitude = if (E2+P4) < 275 -> 10
            elif (E2+P4) > 310 -> 1
            else -> -.2571 * (E2+P4) + 80.7143
```

### Feedback (Hill) functions — MLP  *(checked against Appendix D plots)*
```
E2->FSH negative feedback = (9.5 + (.36*E2)^30)    / (1 + (.36*E2)^30)    + 0.5   # falls 10 -> 1.5, mid E2~2.78
E2->LH  negative feedback = (5   + (.36*E2)^30)    / (1 + (.36*E2)^30)    + 5     # falls 10 -> 6,   mid E2~2.78
FSH->E2 positive feedback = .75 * (.2*FSH)^8       / (1 + (.2*FSH)^8)     + 2.5   # rises 2.5 -> 3.25, mid FSH=5  [* not + : see NOTE 2]
P4->LH  negative feedback = (5   + (.00347*P4)^175)/ (1 + (.00347*P4)^175) + 5    # falls 10 -> 6,   mid P4~288   [offset uncertain: see NOTE 3]
```

> **NOTE 2 — RESOLVED.** The naming is consistent once the form is understood:
> a **negative-feedback** term uses an *additive* numerator constant
> `(A + x^n)/(1 + x^n) + C` with A>1, which **decreases** from `A+C` (low x) to `1+C`
> (high x). A **positive-feedback** term uses a *multiplicative* leading constant
> `A * x^n/(1 + x^n) + C`, which **increases** from `C` to `A+C`. The Appendix D plots
> confirm every term's direction. The one transcription fix: `FSH->E2 positive feedback`
> is printed additively (`.75 + ...`) but its plot **rises** 2.5→3.25, so the leading
> operator must be multiplication (`.75 * ...`), matching the EFP `FSH->E2` form. Half-
> saturation: E2 terms at .36*E2=1 ⇒ E2≈2.78; P4 term at .00347*P4=1 ⇒ P4≈288. High
> exponents (8/30/175) make these near-hard switches.

> **NOTE 3 — P4->LH offset uncertain.** Three sources disagree on the additive constant:
> the high-res scan shows `+2.5`, the garbled text layer shows `-2.5`, and the Appendix D
> plot (range 10→6, identical to E2->LH) implies `+5`. Defaulting to **+5** to reproduce
> the published plot; revisit if MLP dynamics look off.

---

## Comparison dataset (for validation — NOT contained in this PDF)
Model output was compared to: Licinio J, Negrao AB, Mantzoros C, et al. "Synchronicity of
frequently sampled, 24-h concentrations of circulating leptin, luteinizing hormone, and
estradiol in healthy women." *PNAS.* 1998;95:2541–2546.
- 6 healthy women, mean age 25.5 ± 1.6 y, follicular phase
- Blood sampled **every 7 min for 1,442 min** → **207 samples/subject** (8 AM to 8 AM)
- LH by RIA: sensitivity 0.1 U/L; intra-assay CV 2.6%; inter-assay CV 5.4%
- Estradiol by RIA: sensitivity 8 pg/mL; intra-assay CV 7.0%; inter-assay CV 8.1%

## Stated limitations (verbatim sense)
- Model is **qualitative, not fully quantitative**; hormone numbers are not meaningful in
  absolute terms — only the dynamics/shapes are.
- Inhibin and neurohumoral elements are omitted.
- Ovariectomy experiment: removing ovarian steroids returned GnRH/gonadotropins to simple
  periodic behavior (aperiodicity is steroid-feedback-driven "deterministic chaos").

---

## Ambiguities — all resolved (see notes above)
1. **NOTE 1** ✅ EFP `E2->FSH negative feedback = 9/(1+(.0025*E2)^15)+1` (from Appendix C plot);
   GnRH exponent is `^15` both places (printed `^3` was a typo).
2. **NOTE 2** ✅ Feedback direction is encoded by additive (decreasing) vs multiplicative
   (increasing) numerator constants; MLP `FSH->E2` operator fixed `+`→`*`.
3. **NOTE 3** ⚠️ P4->LH offset defaulted to `+5` (plot-consistent); flagged, low risk.
4. **NOTE 4** ✅ `PULSE` = Berkeley Madonna impulse train (area = amplitude); to be
   reimplemented as an explicit, step-size-independent forcing function.
