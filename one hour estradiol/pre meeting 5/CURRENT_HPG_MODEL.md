# Current HPG model: start here

[Exact plotted equations, all21 constants and corresponding raw LH/E2 traces](../meeting11/DENSER_REGIME_EQUATIONS_AND_TRACES.md).
The24h/72h/14day figures use the unchanged tested nominal center at1.5s,
with806,400 full-resolution CSV samples and matching source/audit metadata.
Presentation does not select a new fully qualified regime or retune constants.

**Latest verification: constants-only and active E2→N confirmed; robustness is partial.**
[Fresh uncertainty/source/feedback assessment](../meeting11/DENSER_REGIME_ROBUSTNESS_AND_FEEDBACK.md).
The physical function body is identical to the frozen reference; only18
existing defaults differ. Fixed N reverses day/night raw medians12.39/3.63
to5.17/9.00 at3s/14days. Fresh128 simultaneous±2.5% draws give121/128
local E2 ordering passes,112/128 combined E2/LH/SCN ordering passes and0/128
full passes across both warmups. Absolute amplitudes vary; the full return
is mixed-sign despite continuously inhibitory E2→N. The same-draw finer
repeat gives117/128 local E2 ordering and108/128 combined ordering passes,
with0/128 full passes. E2 ordering passes113/128 on both grids;8 pass3s
but fail1.5s,4 reverse. Finer day pooled medians5.724–15.593 show amplitude
sensitivity. All detailed limitations are recorded in the verification report.

**Latest: more actual daytime E2 pulses, with the approved amplitude scale.**
[Denser raw E2/LH comparison and full audit](../meeting11/MORE_DAYTIME_E2_PULSES.md).
The3s diagnostic has8 eligible daytime peaks/day versus5 previously,
local medians12.891/3.628. The14-day diagnostic has113 versus73 day peaks,
but still fails daytime/control preservation and nighttime concentration
checks. Its feedback is active and mixed-sign, with a restoring finite
midpoint mean-response probe. No qualified replacement or fresh±2.5%
validation is claimed; current search continues with unchanged physical rates.
The finer14-day comparison confirms109 versus72 eligible daytime peaks
(+51.4%),7–8/day, medians12.390/3.628 and finer daily day medians10.781–15.106.
Night maximum31.151 and control retention.38818 remain failures. Full fine
audits and raw first/weakest-day12 comparisons are in the latest report.
The continuation completed2,349 nominal settings in15 families, with no
both-window ordinal joint pass. A+48min SCN-phase alternative passes
14-day concentration containment (night max22.936) but changes the raw
day median to8.187 and still fails control retention.39624. Both variants
remain diagnostics; their full raw/cycle records and differing feedback
evidence are explicit in the report.

**Latest: Licinio individual-pulse contrast remains primary; no qualified replacement.**
[Latest continuation, raw E2/LH plots and full failures](../meeting11/LICINIO_PRIORITY_CONTINUATION.md).
The new14-day diagnostic has pooled raw medians13.412/2.022, but day9 falls
to1.963/1.613 and worst daytime/control retention.18062. Do not promote the
pooled medians or first-day plot. No new±2.5% qualification is claimed.
Original physical equations remain unchanged; older entries below are provenance.

**Latest tolerance: ±2.5% is sufficient; no regime is selected.**
[Current search and rejected nominal waveform](../meeting11/SMALL_RADIUS_ACTIVE_FEEDBACK_SEARCH.md).
The strongest-looking new trial fails longer-window control retention/local
ordering, has a nighttime escape, and gives0/64 fresh joint±2.5% draws
across both warmups. This applies to that candidate, not the model class.
The physical structure and all corrected amplitude/LH/feedback requirements
remain. Older recommendations and percentages below are superseded.

**Latest requirement: retain LH suppression and an active negative return.**
The user requires nighttime-high LH peaks with SCN contribution and is seeking
feedback-mediated irregularity. E2→N inhibition alone is insufficient when
LH suppression makes the full amplitude return reinforcing. No regime is
currently qualified for all corrected behavior/feedback/robustness conditions.
[Correction, measured feedback signs and timing ablation](../meeting11/NIGHT_HIGH_LH_CORRECTION.md).
The active-negative diagnostic has a restoring midpoint response and greater
phase-timing variability with N feedback, but fails nighttime E2 bounds and
the accepted daytime/control loss limit. It is not a replacement calibration.

**Correction: no-suppression recommendation withdrawn.** The user requires
LH peaks highest at night, with SCN contribution. That regime instead shifts
SCN/LH peaks into daytime, so its99.2% statistic omits a required behavior.
It is not a qualified preferred regime. See
[the nighttime-LH correction](../meeting11/NIGHT_HIGH_LH_CORRECTION.md).

**Current choice: remove the decreasing LH suppression Hill; retain E2→N.**
[New assessment, exact inline equations and raw E2/LH/N traces](../meeting11/WITHOUT_LH_SUPPRESSION_ASSESSMENT.md).
The user prioritizes robust day/night differences over a particular mechanism.
The rounded center in `meeting11/without_LH_suppression_trials/MODEL_selected.py`
has local day/night medians13.712/.362pg/mL and actual awake/sleep LH
intervals100.825/139.05min. Fresh simultaneous±10% changes to all19 remaining
constants retain every local mean/median ordering and sleep slowing256/256
at BOTH96/192hwarmups (3s,3days each). Including both approximate reference
concentration spans gives254/256=99.2%, not100% joint robustness.

Day medians still span7.35–18.73, individual day excursions can be2.92;
nominal E2 peak34.06 does not fit the reference48 peak. Day troughs17.44
remain higher than night11.85. Half-lives4/35/22.5 are qualified fitted
alternatives. No whole-box or full reference-waveform qualification is claimed.
The old uninhibited comparator becomes identical, so its strict inequality
is inapplicable; no substitute control is invented. Historical results below
keep their original status and must not override this latest user direction.

**Latest: retain an effective E2→N return path.**
[Active-feedback search and measured limits](../meeting11/MEANINGFUL_FEEDBACK_SEARCH.md).
A constant-only trial gives 128/128 fresh joint ±5% sampled passes at both
warmups (all 21 constants, three-second steps), including E2 comparisons,
reference spans and longer sleep LH. N varies29%; fixing it at its mean
changes the nominal day median10.60→7.18. But absolute day medians span
4.36–16.93; nominal peak37.68 does not fit the full48 reference peak.
At ±10%, E2 passes225/256 (87.9%); including sleep LH slowing187/256 (73.0%).
The full amplitude loop is mostly reinforcing because LH→E2 slope is mostly
negative. Thus no fully qualified ±10% or conventional stabilizing negative-
feedback replacement is established. Feedback-free trials below are not adopted.

**Latest simplification test:**
[E2→N removal and linear LH activation](../meeting11/MODEL_SIMPLIFICATION_ASSESSMENT.md).
Removing E2 feedback leaves nominal local medians effectively unchanged
(9.74/.55), with fresh 19-constant ±10% pass rates 94.1%/96.1% at 15%/20%
allowed day loss. Later-window failures remain. The direct linear-LH trial
has good inhibited pulses but retains only about 12.5% of its control day
amplitude, so it fails the preservation goal. The original reference source
and its previously reported percentages are preserved; the new sources are
explicit trials, not a fully qualified replacement.

[Exact equations and structural/assumption audit of the robustness trial](../meeting11/RANGE_TRIAL_EQUATIONS_AND_AUDIT.md).
No constants changed. Source/center identity and saved traces check exactly;
the audit flags phase-driven GnRH, weak N feedback, low-LH operating-point
sensitivity, the constant-inflow concentration scale and fitted half-lives.

**User-approved daytime loss is now 15–20% versus matched control.**
[Saved ±10% reassessment](../meeting11/DAYTIME_LOSS_TOLERANCE_ASSESSMENT.md):
trial every-day joint pass rates are 245/256 (95.7%) allowing 15% loss and
252/256 (98.4%) allowing 20%. Both concentration spans and nighttime conditions
remain required. Later-window ordering failures and weak absolute day pulses
remain; this is not a full-box or full-horizon guarantee.

**Latest robustness search: no full ±10% regime established.**
[`meeting11/RANGE_ROBUSTNESS_SEARCH.md`](../meeting11/RANGE_ROBUSTNESS_SEARCH.md)
records 7,317 paired search evaluations. A new range-focused trial improves
fresh every-day relative comparisons from 147/256 to 215/256 on identical
seed-1510 directions; all 256 initial trial draws retain nighttime attenuation
and both approximate reference concentration spans. However, daytime raw
medians span 2.74–18.05 pg/mL and a longer-warmup draw reverses median ordering,
including at 3/1.5-second steps. The trial is not a qualified replacement.
The nominal range correction below keeps its limited status.

**Latest nominal candidate: corrected Subject 2 concentration scale.**
[`meeting11/SUBJECT2_RANGE_PRESERVED_BEHAVIOR.md`](../meeting11/SUBJECT2_RANGE_PRESERVED_BEHAVIOR.md)
and [`meeting11/MODEL_subject2_range_candidate.py`](../meeting11/MODEL_subject2_range_candidate.py)
supersede the oversized E2 profiles below for the user's calibration request.
Nominal E2 ranges are 14.69–48.24 day and 14.50–22.67 night pg/mL, versus the
approximate raster reference 9.26–48.16 / 9.16–28.23. Local every-peak medians
are 13.14/4.90, means 12.03/4.56; matched control medians 8.89/11.19. Trough
medians are 17.86/14.87. Actual LH intervals remain 100/135 min awake/sleep.
Only existing constants change; the physical function body remains identical
and fully inline. Half-lives stay at fitted alternatives 4/35/22.5 min.

Nominal comparisons and concentration scale survive finer steps, longer warmup
and a 14-day record. Fresh 256 simultaneous changes across all 21 constants:
±2.5% pooled/every-day comparisons hold 256/256 and 254/256; ±10% holds 173/256
and 139/256. Failures are retained, with no silent loss or range tolerance.
**No fully qualified ±10% regime.** The report includes phone-readable reference/
model traces on a shared E2 axis, all values, missing/tiny pulse counts, frequency,
source verification and numerical checks. Lowest reference values near 9, exact
timing/frequency and the unobserved high control levels remain unfitted.

**±10% assessment:**
[`meeting11/TEN_PERCENT_ROBUSTNESS_ASSESSMENT.md`](../meeting11/TEN_PERCENT_ROBUSTNESS_ASSESSMENT.md).
The current 100/135-minute nominal waveform is not robust to simultaneous
±10% changes across all 21 constants: fresh pooled/every-day comparisons hold
210/256 and 189/256; daytime local medians span 4.38–76.16 pg/mL. Selected
daytime-preservation and night/control failures persist with finer steps and
longer warmup. A 100/120-minute trial improves pooled/every-day counts to
256/256 and 254/256 on a different independent seed, but has daytime medians
5.40–121.61 and a substantially changed waveform (nominal peak 337.82,
daytime trough maximum 142.02). It is not promoted as a replacement. Broader
Hill/kinetic/N-feedback selection found no qualified full-waveform ±10% regime;
this is not an impossibility claim. The nominal candidate below remains intact.

**Preceding comparative candidate, oversized E2 scale:**
[`meeting11/ROBUST_SLEEP_REGIME.md`](../meeting11/ROBUST_SLEEP_REGIME.md) and
[`meeting11/MODEL_robust_sleep_candidate.py`](../meeting11/MODEL_robust_sleep_candidate.py).
Nominal local E2 median day/night excursions are 24.31/5.15 pg/mL, with actual
LH intervals about 100 min awake / 135 during sleep. Only existing constants
change; rates stay fully inline and half-lives remain 4/35/22.5 min.

Fresh paired seed 2026101474, 256 simultaneous changes across all 21 constants:
at ±2.5%, all relative comparisons hold on every day in 256/256 versus previous
239/256. At ±5%, new pooled/every-day counts are 255/256 and 252/256, versus
210/256 and 175/256 before. One pooled draw loses 4% of its daytime median
against control; worst individual-day loss is 17%. These remain failures.
Daytime medians still span 15.34–41.78 at ±2.5% and 9.82–56.80 at ±5%.
**This improves relative margins, but does not establish a robust fixed large
daytime span or a fully qualified ±5% regime.** No unagreed floor/tolerance is
introduced. Longer sleep intervals survive all sampled draws; missing E2
peaks are zero. Nominal behavior survives finer steps, longer warmup and a
14-day record. The selected failing draw still fails at finer steps but passes
after longer warmup. Small morning pulses, elevated midday troughs and high
unfitted control E2 remain material limits. Full raw plots, individual values,
paired comparisons, numerical checks and source verification are in the report.

**Previous longer-sleep-interval option:**
[`meeting11/LONGER_SLEEP_INTERVALS.md`](../meeting11/LONGER_SLEEP_INTERVALS.md) and
[`meeting11/MODEL_longer_sleep_intervals_candidate.py`](../meeting11/MODEL_longer_sleep_intervals_candidate.py).
Actual LH intervals are about 100 min awake / 135 during sleep (previous 100/118).
Nominal every-peak E2 medians 25.02/5.20 day/night retain the preceding scale and
the matched-control comparison. Only four existing sleep/SCN/N constants change;
physical rates remain inline and half-lives stay 4/35/22.5 min. Daytime peak and
trough extremes rise, and early pulses become smaller. Nominal comparisons
survive finer steps and longer warmup. Fresh ±2.5% pooled comparisons hold 63/64,
with one cycle night/control median failure; ±5% holds 55/64. No fully qualified
radius is claimed for this option. The preceding fit remains intact below.

**Preceding larger-span candidate:**
[`meeting11/LARGER_SPAN_LOWER_TROUGHS.md`](../meeting11/LARGER_SPAN_LOWER_TROUGHS.md)
and [`meeting11/MODEL_larger_span_feedback_candidate.py`](../meeting11/MODEL_larger_span_feedback_candidate.py).
The user requested larger excursions and daytime troughs closer to nighttime.
Every-peak median day/night excursions are 25.62/5.21 pg/mL (cycle medians
31.52/5.07), versus matched-control 16.19/9.36. Distinct trough medians are
14.24/7.99, improving the earlier 25.25/8.99; midday troughs still reach 23.17.
Only existing constants were retuned; physical rates remain fully inline.
Clearance alternatives change to GnRH/LH/E2 4/35/22.5 min, with qualified human
sources in the report. E2 frequency slows to 0.611/0.528 per hour. Control mean
E2 is high (331.3/341.8) and is not fitted to an observed record.

All 64 fresh sampled ±2.5% combinations retain pooled mean/median comparisons
for both local measurements; each-day cycle comparisons also hold. Nominal and
three vulnerable draws retain pooled inequalities with finer steps and longer
warmup. At ±5%, only 56/64 retain every pooled comparison: nighttime/control
attenuation has failures. No full ±5% regime, absolute acceptance threshold or
whole-box guarantee is claimed. Raw waveform, all individual values, missing/
small-pulse counts, frequency, source verification and limitations are recorded.
The preceding `larger_span` shape trial is not qualified: its longer-warmup
daytime preservation failed. Historical range-based claims remain withdrawn.

**Previous Subject 2-scale candidate:** [`meeting11/TRUE_LOCAL_PULSE_CANDIDATE.md`](../meeting11/TRUE_LOCAL_PULSE_CANDIDATE.md)
and [`meeting11/MODEL_balanced_hills_pulse_candidate.py`](../meeting11/MODEL_balanced_hills_pulse_candidate.py).
Corrected per-cycle median day/night excursions are 13.39/3.97 pg/mL versus
7.49/8.98 in its fully coupled control. This preserves nominal daytime pulse size
on roughly the Subject 2 scale. All 64 sampled ±2.5% combinations preserve the
relative mean/median inequalities for both local metrics; day medians span
9.48–16.31. ±5% still has control/retention failures. No full ±5% pulse-robust
claim or new absolute acceptance threshold is established. Source and actual
phone-readable traces, full individual records, numerical checks and material
waveform/control limitations are in that report. The earlier results below
remain historical; their withdrawn range-based claims are not reinstated.

For a new chat, read [`NEW_CHAT_HANDOFF.md`](NEW_CHAT_HANDOFF.md) first.

**Updated pulse criterion (October 14): the presentation recommendation is withdrawn.**
Read [`meeting11/TRUE_E2_PULSE_AMPLITUDE_ASSESSMENT.md`](../meeting11/TRUE_E2_PULSE_AMPLITUDE_ASSESSMENT.md)
before using the historical results below. Amplitude is assessed from each E2
peak and its immediately adjacent troughs, with the smaller rise/fall excursion
reported per LH cycle and missing E2 peaks counted as zero. Median day/night
excursions are 0.73/1.90 pg/mL for the former presentation regime and 0.40/1.24
for the rejected higher-production alternative. Both have smaller typical daytime
excursions. This persists with finer steps and longer warmup. The earlier ±6%
and 60/64 at ±15% counts do not establish robustness for the actual pulse target.
Model constants and equations have not been retuned in this correction.

## Historical results using LH-cycle ranges

Broader robustness (October 14):
[`meeting11/BROADER_ROBUSTNESS_RESULTS.md`](../meeting11/BROADER_ROBUSTNESS_RESULTS.md).
The recommended model retains every original amplitude goal in all 64 tested
joint ±6% combinations, with all six vulnerable-case numerical checks passing.
The weakest day range is 9.08 pg/mL, close to its 9.04 retention bound. Each-day
night < day ordering survives all 64 tested ±10% combinations, but only 57/64
retain every original amplitude goal at that radius.
The separate `meeting11/MODEL_broader_amplitude_experiment.py` preserves daytime
retention, the gap, and each-day ordering in all 64 fresh ±15% combinations.
Only 60/64 also retain the earlier absolute-night bound and matched-control
attenuation. Mean daytime E2 rises to 54.17 pg/mL and distinct daytime E2 peaks
become less frequent. The experiment is not a replacement for the recommended
presentation regime. These are sampled results, not whole-box guarantees.

Recommended presentation regime (October 14, after the user accepted constant
retuning of existing N feedback): `meeting11/MODEL_presentation_regime.py`.
Read `meeting11/PRESENTATION_PARAMETER_REGIME.md` for the full nominal table,
equations, primary-source qualifications, and presentation wording. GnRH/LH/E2
half-lives are 4/57/27.5 min; E2 27.5 rounds White et al.'s published IV estimate.
All other coefficients match the rounded day-preserving candidate below. Nominal
mean LH-cycle E2 ranges are 11.11 day / 2.71 night pg/mL, gap 8.40. A fresh seed
tests 64 independent simultaneous +/-5% coefficient/clearance combinations:
all 64 retain the specified day amplitude and night/day gap, with day range
9.39–12.32, night 1.37–4.97, minimum gap 6.86. The weakest four cases pass all
eight finer-step/longer-warmup checks. The full supplemental screen passes 62/64;
two fail only its chosen resolved daytime peak-frequency cutoff. This is sampled
local robustness, not a guarantee throughout the box or a complete human/rat
waveform validation. N relaxation and suppression remain fitted hypothesis
parameters; mean E2 is 35.64 / 9.86 with a pronounced daily envelope. No states
or explicit delays were added. Both earlier hypothesis sources and the canonical
comparison remain unchanged. Historical notes below refer to those earlier fits.

Day-preserving candidate (October 14):
`meeting11/DAY_PRESERVING_ROBUSTNESS_RESULTS.md` and the separate self-contained
`meeting11/MODEL_day_preserving_candidate.py`. Mean cycle amplitudes are
11.79 day / 2.89 night pg/mL, versus 10.04 / 8.57 in the active prior fit.
All 32 independent joint +/-5% draws, including all three clearance half-lives,
preserve the stated daytime-retention and amplitude-gap goals; the weakest
cases also pass finer-step and longer-warmup checks. Three of 32 +/-10% draws
fail daytime retention, and wider LH/GnRH clearance alternatives still fail.
This is bounded sampled local robustness. It requires fitted N relaxation
7.69 h, LH half-life 57 min, retuned existing coefficients, and a pronounced
daily E2 envelope. Mean E2 becomes 38.77 / 10.87 day/night, and LH inhibition
acts during day too. No new state or helper function is added. The active prior
fit and canonical comparison remain unchanged; do not silently substitute this
candidate or call it a full reference-waveform fit.

Joint changes (October 14): `meeting11/JOINT_PARAMETER_ROBUSTNESS.md`.
With 16 endocrine coefficients varied independently together within 1%, 5%,
and 10%, the three-day average night < day amplitude ordering survives
24/32, 16/32, and 10/32 sampled combinations at six-second resolution.
Coordinated changes can compensate: source +1% (13 to 13.13) and inhibition
exponent -1% (3 to 2.97) restore cycle ordering on all nine days across the
standard/finer-step/longer-warmup checks. Local excursions do not pass every
check. This pair is an experiment and was not promoted; the active defaults
remain unchanged. Sample fractions are not physiological probabilities.

Other-parameter robustness (October 14): `meeting11/OTHER_PARAMETER_ROBUSTNESS.md`.
Across 105 one-at-a-time changes in 18 endocrine coefficients with half-lives
fixed, night E2 remains 36.7--62.4% below its matched uninhibited model control.
The stronger night < day ordering is fragile: no tested symmetric +/-5% band
preserves it on every day. The promising +/-1% N-gain and minimum-SCN bands
also fail added three-second numerical checks. The active coefficients and
saved current fit were not changed.

Literature and equation audit: `meeting11/HALF_LIFE_AND_EQUATION_AUDIT.md`.
GnRH 4 min and E2 30 min have qualified human-study support; LH 20 min
approximates a reported fast disappearance phase, not the longer endogenous
clearance estimates in cycling women. The current night < day amplitude
ordering is sensitive to clearance: changing only LH to 21 min or E2 to
27.45 min reverses it. The inherited two fifth-power central Hills also create
steep effective feedback. These limitations do not change the saved fit.

The current working hypothesis is
`meeting11/MODEL_smooth_lh_suppression.py`. Read
`meeting11/CONSTANT_INFLOW_RESCUE_RESULTS.md` for its inline equation, checked
amplitudes, controls, and time series. The user is trying to preserve nighttime
attenuation with ordinary clearance, constant input, and modest Hill functions.
The current equation is
dE2/dt = 13 + 65*LH^4/(1+LH^4)*3^3/(3^3+LH^3) - 2*ln(2)*E2.
The constant source is 13 pg/mL/h, production scale 65, activation exponent 4,
inhibition midpoint 3 IU/L and exponent 3, and E2 half-life 30 minutes. Both LH
Hill exponents are restricted to 1--4. There is no capacity state or recovery
equation; the only functions are the two simulation entry points with inline rates.

This recovers the qualitative attenuation and ordering: standard E2 full-cycle
ranges are 10.04 day / 8.57 night pg/mL. Night range is below day on all nine
retained days across 6-second, 3-second, and longer-warmup checks, and 54--59%
below a matched uninhibited model control. Day cycle range retains 100--106% of
the original comparison baseline. Local excursions are 5.84 / 5.72 pg/mL, with
local daytime size retaining only 62--64% of the original. The stricter 90%
retention criterion by both amplitude measures is not met. The inhibitory
midpoint is reached during day too; do not call this exclusively night-triggered.
It is a qualitative fitted hypothesis, not a complete reference-waveform fit.

Sleep remains a binary 23:00--07:00 input with period multiplier 1+0.75*Sleep.
GnRH, N feedback, LH clearance, SCN, and thermal coefficients are unchanged.
The amplitude recovery required retuning the existing LH activation and inhibition
parameters and the E2 half-life alongside constant inflow. Constant input alone
was not sufficient. All reported acceptance checks use the standalone scalar
model with matched controls; the active defaults exactly match the saved trajectory.
`meeting11/constant_inflow_rescue_trials/` contains the screens, validation,
trajectories, fully coupled and fixed-LH controls, PNG/SVG plots, and a CSV.

The previous no-constant-input model is `MODEL_regular_clearance_history.py`,
with historical results in `INSTANTANEOUS_LH_HILL_RESULTS.md`. The capacity-state
version is `MODEL_smooth_capacity_history.py`, with historical results in
`SMOOTH_LH_SUPPRESSION_RESULTS.md`. The canonical `MODEL.py` remains the original
comparison baseline with unchanged SHA256.

The development history below documents the original baseline and preceding trials.

Before modifying it, read:

1. `meeting11/README.md` for the equations, parameters, and current behavior.
2. `meeting11/MODEL_DEVELOPMENT_LOG.md` for the agreed constraints, failed
   approaches, numerical comparisons, and the unresolved E2-amplitude issue.
3. `meeting11/test_model.py` for the behaviors currently enforced.

The original baseline's unresolved result was that its instantaneous, monotone
LH-to-E2 Hill inflow can produce lower nighttime E2 levels but not smaller
nighttime peak-to-trough E2 excursions while nighttime LH pulses remain
larger and less frequent. Parameter sweeps, a globally shorter LH half-life,
an approximately five-hour nighttime LH interval, and a fixed 30-minute
LH-to-E2 delay have already been tested. Do not repeat those trials without
first reviewing the development log.

The continuing investigation is recorded in
`meeting11/JOINT_CALIBRATION_FINDINGS.md` and
`meeting11/LICINIO_CONSISTENCY_AUDIT.md`. The user requires additions to be
justified and consistent across Licinio's measurements. A numerically improved
shared-SCN ovarian prototype has an unsupported human nocturnal phase and is
not a validated replacement. The canonical model remains `meeting11/MODEL.py`.

The earlier user-authorized hypothesis experiment adds an instantaneous
decreasing LH Hill multiplier to E2 inflow. See
`meeting11/LH_HILL_INHIBITION_RESULTS.md` and
`meeting11/MODEL_lh_hill_inhibition_trial.py`. Only the inhibitory parameters
were swept. The checked trial default (K_I = 2.5 IU/L, m = 32, residual = 0.375)
reduces nighttime full-LH-cycle E2 variation by about 51–57 percent, retaining
about 75–79 percent of daytime variation across checks. The small nighttime-
below-daytime amplitude ordering is not robust to step size and warmup. The
plain Hill tending to zero achieves about 37–39 percent nighttime attenuation
with about 73–77 percent daytime retention and reverses the desired E2 mean
ordering. These are hypothesis trials; the canonical model remains the baseline.

An earlier user clarification required night E2 amplitude below day, rather
than merely attenuation relative to the old nighttime trace. That trial's
defaults are K_I = 2.07 IU/L, m = 37.5, residual = 0.4125. See
`meeting11/LH_HILL_AMPLITUDE_ORDERING_RESULTS.md`. With only these inhibitory
parameters changed, local E2 excursions are 3.33 day / 2.32 night pg/mL and
whole-LH-cycle E2 ranges are 5.51 day / 5.34 night pg/mL at six-second resolution.
Both amplitude orderings hold at three-second steps and after 192 h warmup;
whole-cycle night ranges are lower on each of the nine retained days across
these checks. This is a modest full-cycle reversal with reduced daytime range,
not a complete Subject 2 fit. No additional state or mechanism was introduced.
`meeting11/MODEL_lh_hill_inhibition_trial.py` retains this earlier ordered setting.

An earlier user clarification required preserving daytime E2 amplitude and proposed
desensitization only above an LH threshold absent from daytime LH. The
then-preferred hypothesis trial was `meeting11/MODEL_threshold_desensitization_trial.py`;
see `meeting11/THRESHOLD_DAYTIME_PRESERVATION_RESULTS.md`. Its fixed threshold is
5.5 IU/L and no daytime LH crosses it in any of the three numerical checks.
One intracellular capacity state permits suppression after threshold crossings
and recovery afterward, with fitted times of 10 and 80 minutes. All pre-existing
HPG coefficients stay fixed. At six-second resolution, E2 full-cycle ranges are
9.80 day / 6.63 night pg/mL (original baseline: 9.69 / 16.77), and local excursions
are 8.79 / 3.66 pg/mL. Day full-cycle amplitude retains 101–104 percent of baseline
and local day excursions retain 94–99 percent across checks. Night amplitude is
smaller on every retained day. Only night triggers desensitization, though some
recovering capacity carries into early day. Final acceptance uses the standalone
model because bulk floating-point evaluation can shift hybrid pulse phase.
This superseded the reduced-daytime-amplitude instantaneous trial. The user then
requested the smooth simplification described at the top of this file.
The canonical `MODEL.py` remains the comparison baseline.

