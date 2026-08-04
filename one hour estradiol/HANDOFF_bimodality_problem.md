# Handoff: reproducing mouse core-body-temperature bimodality

**Status 2026-08-03.** Unsolved. This document states the problem, what the data actually shows,
what has been ruled out, and — most importantly — the methodological traps that produced five
false "successes" in a single session. Read the traps section before trusting any metric.

---

## 1. The problem in one paragraph

Three mice, core body temperature (CBT) recorded over gestational days 13-18, `n = 1 per
condition`:

| animal | E2 signalling | P4 | observed daily pattern |
|---|---|---|---|
| non-preg | intact receptor, cycling | none | **bimodal**, peaks 8.2 & 15.2 h |
| pregnant | intact receptor, high E2 | high | **bimodal**, peaks 15.8 & 21.8 h |
| Esr1i | **ERα blocked** | high | **unimodal**, single peak 13.2 h |

Build a model that reproduces **bimodal / bimodal / unimodal** with the right amplitudes and the
right peak positions. Nothing tried so far achieves all three at once.

## 2. The model being fitted (Holly's form — do not silently change it)

```
c(t) = cos(2*pi*(t - 12)/24)                    ONE clock, shared, peak locked at NOON
loss = sig(bL + gL*c + wL*E2 + pL*P4)           heat-DUMP arm
keep = sig(bK + gK*c + wK*E2 + pK*P4)           heat-KEEP arm
CBT  = T0 + A_T*(keep - loss)
```

* **All coefficients positive** (clock, E2 and P4 all EXCITE both arms). A mixed-sign version was
  considered and rejected by the author.
* `E2 = 0 / 1 / 0` and `P4 = 0 / line / line` for non-preg / pregnant / Esr1i. P4 is a literal
  straight line `P4(d) = 1 - (d-13)/6.5`.
* Note the consequence: **non-preg and Esr1i differ ONLY by P4.** Any mechanism that is to make
  Esr1i unimodal while non-preg is bimodal must act through P4, or the encoding must change.
* `A_T <= ~10 C` is a physical bound. Fits needing more are unphysical.

## 3. What the data actually shows (verify before re-deriving)

Run `thermo_scripts/series.npy` per-day, not the 5-day fold.

* **The bimodality is REAL, not a folding artifact.** Peaks recur at the same clock hour on
  independent days: non-preg first peak at **6.9 / 7.1 / 7.2 h**; pregnant second peak at
  **21.8 / 22.1 / 22.4 h**. Noise cannot land within 20 minutes of the same time of day across
  separate days.
* **But it is INTERMITTENT**: 3-4 of 5 days show two peaks, not 5 of 5. Do not demand an
  every-day bimodality the animal does not have.
* **Esr1i's unimodality is robust**: 0 of 5 days show a second peak at prominence >= 0.20.
* **The fold's peak positions are BIASED.** non-preg's true recurring first peak is ~7.0 h, but
  the fold says 8.2 h because days 13-14 (single-peaked, at 7.9 and 10.3) drag the average.
  **Fit targets should be re-derived from per-day peaks, not the fold.** All work to date used
  the fold positions.
* The "2/2/1" description holds only in a narrow prominence window around **0.20**
  (0.15 -> 3/3/1, 0.25 -> 1/2/1). Do not demand more threshold-robustness than the data has.

### The decomposition that dissolves the central obstacle

Measured from the per-day raw traces:

| animal | mean amplitude | mean absolute dip | **dip/amp** |
|---|---|---|---|
| non-preg | 2.45 C | 0.880 C | **0.359** |
| pregnant | 0.82 C | 0.265 C | **0.321** |
| Esr1i | 1.23 C | 0.108 C | **0.088** |

Two *independent* patterns:

* **Amplitude** is monotone in hormone LOAD: non-preg (none) > Esr1i (P4) > pregnant (P4+E2).
  Esr1i is the **middle**.
* **dip/amp** tracks ERα FUNCTION and is essentially **binary**: both receptor-intact animals
  0.32-0.36, the blocked one 0.088. Esr1i is the **extreme**. Note pregnant's dip/amp (0.321) is
  nearly identical to non-preg's (0.359) despite a 3x amplitude difference — **notch depth does
  not care about hormone level, only receptor function.**

This matters because the long-standing blocker was "Esr1i must be simultaneously inside and
outside the bracket". It is not a contradiction: it is middle in one quantity and extreme in a
*different* quantity. Any model that routes both through a single bias axis will fail.

## 4. Ruled out (with the reason — do not retry blindly)

| attempt | outcome |
|---|---|
| P4 in both arms | Esr1i cannot be unimodal — all-positive coeffs bracket it |
| P4 in the KEEP arm only | breaks the bracket, but fits were **classifier artifacts** (see traps) |
| Two-phase clock (keep/loss arms at different phases) | fit collapsed the gap to 0.2-0.4 h, i.e. never used it |
| E2 clock phase shift `dphi` | pre-registered prediction +7 h, fitted +1.88 h, objective moved 0.2%. **BUT the scan that "rejected" it re-optimised the base phase at each pinned `dphi`, so the phase absorbed the shift — this test is INVALID and the hypothesis is untested** |
| ERα gate on loss-arm clock drive alone | matched every summary statistic with waveform correlation **-0.12 / +0.19 / -0.17** — uncorrelated with the data |
| saturating CVC | needs A_T ~ 130 C |
| thermal inertia | needs 12x tau |

## 5. Best current model (`thermo_scripts/_lambda_notch2.py`)

```
CBT = T0 + A_T*(keep - lam*R*loss)
R   = 1 (non-preg, pregnant)   rho (Esr1i)
lam = 0.62      rho = 0.31      ONE shared clock phase (fitted ~12.6 h, i.e. near noon)
```

**`lam` is the one idea that worked.** Both arms were full-range sigmoids (0->1), so the dump arm
subtracted the *entire* keep signal and carved a canyon to baseline — the "two peaks" were only
its shoulders. Holly spotted this ("pregnancy has two troughs not two peaks"). With `lam < 1` the
trough sits at `1-lam`, so amplitude stays full while the notch is partial. **This is also why
every earlier fold-bimodality was amplitude-poor: they all had `lam = 1` implicitly.**

Result: correct **2/2/1 modality** with correct notch depths (0.374/0.341/0.096 vs
0.359/0.321/0.088), and:

| animal | correlation | RMS | animal's own day-to-day SD | verdict |
|---|---|---|---|---|
| non-preg | +0.58 | 0.617 | 0.344 | **1.8x the animal's spread** |
| pregnant | +0.10 | 0.236 | 0.265 | "within" — but see below |
| Esr1i | +0.76 | 0.366 | 0.239 | **1.5x the animal's spread** |

## 6. THE MAIN PROBLEM — read this part

**Pregnant's waveform cannot be fitted, and it is a SHAPE problem, not a timing problem.**

* Observed pregnant: flat/low from 0-12 h, broad rise peaking **15.8 h**, dip at 19 h, second rise
  at **21.8 h**. Two humps clustered in the *late* day, 6 h apart.
* The model produces a spike at ~3-5 h and another at ~22 h, with a long flat shelf between —
  two spikes on *opposite sides* of the day.
* **A circular-distance metric cannot tell these apart.** Peaks at 21.8 & 3.4 h score the same
  "6 h gap" as peaks at 15.8 & 21.8 h. This metric was used throughout and hid the discrepancy
  for a long time. **Judge pregnant by eye and by correlation, never by peak gap.**

Why it is hard: non-preg's peaks are centred at **11.7 h**, pregnant's at **18.8 h** — the whole
rhythm sits ~7.1 h later. With one shared clock phase, the model physically cannot place both, so
it sacrifices one animal or the other. Run 1 (shape-only objective) fitted pregnant at r=+0.83 and
got the other two's modality backwards; run 2 (shape + notch depth) fixed the modality and
destroyed pregnant. **This ~7 h offset is the long-standing "unexplained 7.5 h phase shift" and it
is the thing blocking all three animals from fitting simultaneously.**

The obvious fix (E2 shifts the clock phase) is **biologically well grounded** — estradiol
advances activity onset and shortens period in rodents — and **has not actually been tested
correctly** (see the invalid scan above). *That is the first thing the next agent should do:* pin
pregnant's effective phase to ~18.8 h directly, rather than fitting `dphi` with the base phase
free.

### Two further open problems

1. **The model is a square wave; the mouse is not.** Fits need deep sigmoid saturation to reach
   the amplitudes, and saturation squares the waveform. Real flatness is ~0.61 (roughly
   sinusoidal). Visible in `fig_vs_raw.png`.
2. **No day-to-day variation.** The model emits one fixed shape forever; the mouse varies a lot
   (one broad hump on day 14, two on days 15-17). It is being fitted to *the average of
   qualitatively different days*.

## 7. METHODOLOGICAL TRAPS — the most valuable part of this handoff

Five "successes" were reported and then retracted in one session. Every one was the same error:
**a target that can be satisfied without being right.** An optimiser under pressure will always
find the cheapest way to satisfy a label.

1. **Four hand-rolled peak classifiers were gamed, in four different ways.**
   * continuous sign-flip test -> fit put the two "peaks" **0.034 h (2 min) apart**
   * discrete grid counter -> the "2nd peak" was a **tie-break artifact at the trough**; the curve
     was a monotone rise-then-fall
   * single prominence cutoff 0.10 -> dips landed at **exactly 10.0%**, flipping at 9%
   * "margin band" meant to fix that -> the dip metric returned 0 when `nmax != 2`, so a
     **three-peak** curve passed the *unimodal* test and the script printed `SATISFIED: True`
   > **Use one exact routine** (scipy prominence on a tripled array for circular data) and
   > **validate it on the OBSERVED data first**. Doing that here immediately revealed the narrow
   > 0.20 threshold window. **Prefer targeting POSITIONS over COUNTS.**

2. **Flat-line degeneracy.** With `A_T` fitted by least squares, one fit "solved" the topology on a
   **0.12 C flat line** (observed 2.10 C). Fix: define `A_T := AMP_obs/H-range` so amplitude is
   exact by construction and a flat solution is rejected outright.

3. **Summary statistics without shape.** Weighting dips 6x and ratios 4x against RMS 1x produced a
   fit matching *every* statistic with correlation **-0.12 / +0.19 / -0.17**. The recorded rule
   "fit on amplitudes, never RMS alone" is true, but so is its converse.
   > **ALWAYS PLOT THE TIME SERIES BEFORE CLAIMING A FIT.**

4. **Circular metrics hide waveform differences** — see section 6.

5. **Fitting the fold flatters the model.** Judge against the RAW trace, and use the animal's own
   **day-to-day SD** as the bar: if model RMS exceeds it, the model predicts the animal worse than
   another day of that same animal would. Two of three currently fail this.
   Beware also that pregnant's small amplitude (0.72 C) means a near-flat curve scores a decent
   RMS while having correlation +0.10 — **low RMS there means "close to flat", not "right".**

## 8. Files

| file | what |
|---|---|
| `thermo_scripts/series.npy` | raw traces (use this, not the fold) |
| `thermo_scripts/fold.npy` | 5-day fold — biased peak positions, see §3 |
| `thermo_scripts/_lambda_notch2.py` | **best current model** |
| `thermo_scripts/_revert_baseline.py` | refits Holly's original form as a clean baseline |
| `thermo_scripts/_compare_raw.py` | model vs raw trace + day-to-day SD test — **the honest check** |
| `thermo_scripts/fig_vs_raw.png` | the comparison that should go to a PI, not the fold plots |
| `thermo_scripts/_dphi_scan.py` | the INVALID phase scan (§4) — fix before reusing |
| `mouse_cbt_folded.csv` | digitised data |

## 9. Caveats that bound every claim here

* **`n = 1` animal per condition.** Everything is within-animal reproducibility.
* What `Esr1i` actually is (global vs conditional knockout) was never established, and it changes
  the interpretation.
* The light schedule / zeitgeber for the three animals was never checked against the source paper.
  This decides whether the ~7 h offset is biology or a recording offset — **it is the single most
  valuable thing to resolve, because the main problem in §6 depends on it.**
