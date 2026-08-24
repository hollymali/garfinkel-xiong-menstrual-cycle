# Final endocrine–temperature model

This is the finalized dynamic model used for the mouse-age body-temperature
comparison. It keeps the circadian-modulated endocrine model, a simple
conservation-versus-loss heat balance, and a separate temperature-activated
self-loss pathway. There is no temperature set point and no relaxation-time
parameter.

## Endocrine model

The endocrine circadian multiplier is

\[
C_{E2}(t)=\cos\left(2\pi[m(t)-840]/1440\right),
\qquad
M_{E2}(t)=1+0.7\frac{1+C_{E2}(t)}{2},
\]

where time is in minutes. Thus pulse amplitude is multiplied by a value from
1.0 to 1.7 and the endocrine clock peaks at 14:00.

\[
\frac{dN}{dt}=\frac{0.25}{1+(E_2/20.1)^5}-0.25N,
\]

\[
A_N(t)=6.8M_{E2}(t)\frac{N^4}{0.28^4+N^4},
\qquad
P_N=63+51\frac{N^2}{0.28^2+N^2}\ \mathrm{min},
\]

\[
\frac{dGnRH}{dt}=-\frac{\ln2}{4}GnRH+A_N\delta_{pulse},
\]

\[
\frac{dLH}{dt}=0.5GnRH-\frac{\ln2}{20}LH,
\qquad
\frac{dE_2}{dt}=0.1145LH-\frac{\ln2}{30}E_2.
\]

## Thermal signals

E2 and progesterone enter through one linear coefficient each. Neither hormone
has its own Hill function. The outer conservation and loss pathways remain
bounded by the same simple fourth-order Hill readout.

The shared thermal clock is

\[
C_T(t)=\frac{1+\cos\left(2\pi[h(t)-18.16]/24\right)}{2}.
\]

Its center was shifted from 11.71 to 18.16 hours to align the modeled daily
pattern with the mouse-age body-temperature traces. The upstream endocrine
clock remains unchanged at 14:00.

With \(H_4(x)=x^4/(1+x^4)\),

\[
\mathrm{Conservation}
=H_4\left(2.27C_T+0.60N+0.0274P\right),
\]

\[
\mathrm{HormonalLoss}
=H_4\left(0.78C_T+0.0000869E_2\right).
\]

Here P is measured in ng/mL and E2 in pg/mL. The E2 coefficient is small only
because raw E2 concentrations are numerically much larger than the unit-range
clock input.

## Final temperature equation

\[
\boxed{
\frac{dT}{dt}
=49.44
+9.76\,\mathrm{Conservation}
-15.29\,\mathrm{HormonalLoss}
-125\frac{T^5}{39.08^5+T^5}
}
\]

with \(T(0)=37.0\ ^\circ\mathrm C\). Rates are in degrees Celsius per hour
after scaling the heat balance by whole-body heat capacity.

Only simple coefficients were retuned after removing the two hormone-specific
Hill transforms. The thermal-clock center, outer pathway exponent, self-loss
exponent, self-loss maximum, and self-loss half-saturation value were left
unchanged from the accepted mouse-age fit.

## Hormone conditions

- Cycling: dynamic \(N(t)\) and \(E_2(t)\) from the endocrine model;
  \(P=0.2\) ng/mL.
- Constant low E2: \(E_2=4.865\) pg/mL, the exact minimum of the retained
  cycling trace; \(P=0.2\) ng/mL; \(N\) is held at its E2-dependent steady
  state.
- Pregnancy: \(E_2=215\) pg/mL, \(P=25\) ng/mL; \(N\) is held at its
  E2-dependent steady state.

## Final behavior

- Cycling: peak approximately 37.94 °C; trough approximately 35.90 °C.
- Constant low E2: peak approximately 37.95 °C; trough approximately 36.18 °C.
- Pregnancy: peak approximately 37.96 °C; trough approximately 36.33 °C.
- Peaks remain aligned, while the cycling-E2 trough remains the deepest.
- Removing the hormone-specific Hill transforms changes cycling RMSE by less
  than 0.001 °C and constant-low RMSE by less than 0.002 °C.

## Biological basis of the temperature-dependent loss term

The self-loss Hill term is a bounded phenomenological approximation of the
combined heat-dissipation response, not a claim that a single receptor follows
fifth-order binding kinetics. Experimental support for a temperature-activated,
nonlinear, and capacity-limited loss response includes:

- Cabanac and Massonnet measured sweating and vasomotor heat loss as functions
  of core temperature in humans
  ([J Physiol, 1977](https://pubmed.ncbi.nlm.nih.gov/856985/)).
- Fujii et al. measured core-temperature thresholds and post-threshold slopes
  for cutaneous vasodilation and sweating
  ([J Appl Physiol, 2014](https://pubmed.ncbi.nlm.nih.gov/25257867/)).
- Minson et al. found that local heating produced a finite cutaneous vascular
  conductance plateau, with age-dependent attenuation
  ([J Appl Physiol, 2002](https://pubmed.ncbi.nlm.nih.gov/12381749/)).
- Abbott and Saper showed in mice that activating warm-responsive median
  preoptic circuitry produces cutaneous vasodilation and hypothermia
  ([J Physiol, 2017](https://pubmed.ncbi.nlm.nih.gov/28786483/)).
- Krajewski-Hall et al. showed that activating NK3R pathways in female mice
  raises tail temperature while lowering core temperature, and that severe
  warmth produces nearly constant vasodilation
  ([Temperature, 2017](https://pubmed.ncbi.nlm.nih.gov/29687044/)).
- Skop et al. quantified mouse whole-body heat conductance and showed that the
  tail is only one component of total heat dissipation
  ([Am J Physiol Endocrinol Metab, 2020](https://pubmed.ncbi.nlm.nih.gov/32691633/)).

These studies justify the direction and bounded shape of the term. They do not
independently identify the fitted exponent \(n=5\) or half-saturation
\(K_T=39.08\ ^\circ\mathrm C\); those are effective model parameters.

## Files

- `MODEL.py`: endocrine simulation.
- `dynamic_temp_model.py`: canonical final temperature model and exports.
- `plot_temp_timeseries.py`: plots the three final model conditions.
- `fit_age_temperature.js`: reproducible constrained fit to the digitized
  3-month and 24-month female mouse traces.
- `mice_age_temperature_digitized.csv`: digitized age-comparison data.
- `final_age_temperature_comparison.png`: final mouse-data comparison.

Run `python dynamic_temp_model.py`, followed by
`python plot_temp_timeseries.py`, to regenerate the final outputs.
