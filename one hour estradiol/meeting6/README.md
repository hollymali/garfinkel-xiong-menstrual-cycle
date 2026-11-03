# Meeting 6 final endocrine-temperature model

This folder contains the final model, the two final figures, and the source
images/data used for comparison. The solver is deliberately compact and all
plotting is isolated in `plot_model.py`.

## Equations

### 1. Endocrine clock

All model time is measured in hours. The endocrine clock peaks at 14:00:

$$
X_{E2}(t)
=
\frac{1+\cos\left(
2\pi\frac{(t \bmod 24)-14}{24}
\right)}{2}
$$

$$
M_{E2}(t)
=
1+0.7X_{E2}(t)
$$

Thus, the unit clock $X_{E2}$ varies from $0$ to $1$, while the multiplier
$M_{E2}$ returned by `endocrine_clock()` varies from $1.0$ to $1.7$.

Biologically, the $N$-dependent factor represents the basal arcuate KNDy
pulse generator, while $X_{E2}$ represents VIP-associated circadian
modulation. The baseline $1$ in $M_{E2}=1+0.7X_{E2}$ therefore allows GnRH
pulses to continue when the modeled circadian modulation is minimal. This is
consistent with pulsatile LH persisting after SCN lesions
([Arendash and Gallo, 1979](https://pubmed.ncbi.nlm.nih.gov/379679/);
[Soper and Weick, 1980](https://pubmed.ncbi.nlm.nih.gov/7349965/)) and with
direct evidence that arcuate KNDy neurons generate GnRH/LH pulses
([Nagae et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7865162/)).
$X_{E2}=0$ should be interpreted as minimum circadian modulation, not proven
complete absence of VIP. The gain $0.7$ is a phenomenological model parameter,
not an experimentally established biological constant.

### 2. E2 feedback and pulse generation

$$
\frac{dN}{dt}
=
15\left(\frac{20.1^5}{20.1^5+E_2^5}\right)
-15N
$$

The production and turnover coefficient for $N$ is $15\ \mathrm{h}^{-1}$.

Pulse amplitude:

$$
A_N(t)
=
6.8\,M_{E2}(t)
\left(\frac{N^4}{0.28^4+N^4}\right)
$$

Pulse period:

$$
P_N(N)
=
1.05
+0.85\left(\frac{N^2}{0.28^2+N^2}\right)
\quad\text{hours}
$$

At each pulse time $t_p$, GnRH receives an instantaneous jump:

$$
\mathrm{GnRH}(t_p^+)
=
\mathrm{GnRH}(t_p^-)+A_N(t_p)
$$

### 3. Hormone cascade between pulses

The half-lives and their corresponding first-order decay constants are

$$
\begin{aligned}
t_{1/2,\mathrm{GnRH}} &= \frac{1}{15}\ \mathrm{h},
& k_{\mathrm{GnRH}}
&= \frac{\ln(2)}{1/15}
=15\ln(2)\ \mathrm{h}^{-1}, \\
t_{1/2,\mathrm{LH}} &= \frac{1}{3}\ \mathrm{h},
& k_{\mathrm{LH}}
&= \frac{\ln(2)}{1/3}
=3\ln(2)\ \mathrm{h}^{-1}, \\
t_{1/2,E_2} &= \frac{1}{2}\ \mathrm{h},
& k_{E_2}
&= \frac{\ln(2)}{1/2}
=2\ln(2)\ \mathrm{h}^{-1}.
\end{aligned}
$$

Using those decay constants, the cascade is

$$
\frac{d\,\mathrm{GnRH}}{dt}
=
-k_{\mathrm{GnRH}}\,\mathrm{GnRH}
$$

$$
\frac{d\,\mathrm{LH}}{dt}
=
30\,\mathrm{GnRH}
-k_{\mathrm{LH}}\,\mathrm{LH}
$$

$$
\frac{dE_2}{dt}
=
6.87\,\mathrm{LH}
-k_{E_2}\,E_2
$$

### 4. Thermal pathways

The thermal clock is a separate biological input from the endocrine pulse
multiplier. It peaks at $14.00$ h and ranges from $0$ to $1$:

$$
X_T(t)
=
\frac{1+\cos\left(
2\pi\frac{(t \bmod 24)-14.00}{24}
\right)}{2}
$$

Heat conservation:

$$
\mathrm{Conservation}(t)
=
\frac{
\left(2.365X_T(t)+0.932N(t)+0.0379P\right)^4
}{
1+\left(2.365_T(t)+0.932N(t)+0.0379P\right)^4
}
$$

Hormonal heat loss:

$$
\mathrm{HormonalLoss}(t)
=
\frac{
\left(0.58X_T(t)+0.000566E_2(t)\right)^4
}{
1+\left(0.58X_T(t)+0.000566E_2(t)\right)^4
}
$$

Here progesterone $P$ is in ng/mL and $E_2$ is in pg/mL.

### 5. Body-temperature equation

$$
\boxed{
\frac{dT}{dt}
=
48.81
+8.97\,\mathrm{Conservation}
-24.98\,\mathrm{HormonalLoss}
-125\left(\frac{T^5}{39.08^5+T^5}\right)
}
$$

The initial condition is

$$
T(0)=37.0\ ^\circ\mathrm{C}
$$

Temperature rates are in degrees Celsius per hour after whole-body
heat-capacity scaling. The last term is a bounded,
phenomenological temperature-activated heat-loss response; its fitted Hill
exponent is not interpreted as a receptor-binding stoichiometry.

## Conditions

- **Cycling:** dynamic $N(t)$ and $E_2(t)$, with $P=0.2$ ng/mL.
- **Constant low E2:** $E_2$ is the minimum of the retained cycling trace,
  $P=0.2$ ng/mL, and $N$ is fixed at its E2-dependent steady state.
- **Pregnancy:** $E_2=215$ pg/mL, $P=25$ ng/mL, and $N$ is fixed at its
  E2-dependent steady state.

## Solver and plots

`MODEL.py` uses one shared $1/600$-hour (six-second) grid and the same explicit
Euler style for both the endocrine cascade and temperature. Every calculated
time point passes directly into the thermal model and returned series.
The first 96 hours are discarded as warm-up. `run()` returns the final 120-hour
series directly to the plotting module; no intermediate CSV or JSON files are
needed.

## Temperature fit

The accepted comparison metrics are:

- Young/cycling mouse RMSE: 0.419 °C.
- Aged/constant-low-E2 mouse RMSE: 0.557 °C.
- Cycling temperature: peak 37.946 °C; trough 35.902 °C.
- Constant-low-E2 temperature: peak 37.949 °C; trough 36.182 °C.
- Pregnancy temperature: peak 37.959 °C; trough 36.329 °C.

Regenerate all final figures with:

```powershell
python plot_model.py
```

The plotting command imports `run()` from `MODEL.py` and writes:

- `final_model_timeseries.png`: the three final modeled conditions.
- `final_e2_timeseries.png`: a 24-hour cycling E2 trace from 08:00 to 08:00.
- `final_temperature_reference_comparison.png`: three five-day comparisons:
  cycling E2 against the non-pregnant `mousebodytemp` trace, constant low E2
  against the 24-month female 2011 trace, and pregnancy/high E2 against the
  pregnant-control `mousebodytemp` trace.
- `cycling_e2_2011_clock_anchored_comparison.png`: cycling E2 against the
  3-month female 2011 trace, with both five-day series beginning at 06:00/ZT0.
- `manuscript cycling e2 comparison.png`: the cycling-E2 model against the
  non-pregnant `mousebodytemp` trace, waveform-aligned over the five-day span.
- `cycling e2 2011 waveform aligned.png`: the cycling-E2 model against the
  3-month female 2011 trace, phase-aligned for a shape-only comparison.
- `constant low e2 2011 waveform aligned.png`: the constant-low-E2 model
  against the 24-month female 2011 trace, phase-aligned for shape comparison.
- `model e2 condition amplitude comparison.png`: the constant-low, cycling,
  and constant-high/pregnancy-hormone model traces overlaid on the same axes;
  the legend reports each five-day peak-to-trough amplitude.

The two `mousebodytemp` comparisons use correlation to select the best phase
shift within one 24-hour cycle. Those reported phase shifts are plotting
transformations for waveform comparison, not fitted biological delays. Their
x-axes are relative elapsed days; integer ticks are not called midnight
boundaries without the original acquisition timestamps or a stated protocol.

The constant-low-E2 comparison is not phase fitted. It uses the 2011 clock
anchor directly: each integer day begins at lights-on, 06:00 (ZT0), and the
shaded halves of the days are the 18:00-06:00 dark phases.

The separate cycling-E2/2011 comparison uses the same clock anchor and is also
not phase fitted. `run(output_start_clock_h=6.0)` extends the solver before
retaining its five-day output, so this view does not circularly wrap the model
or reuse the first day's values at the end.

The 2011 experiment itself did have a clock anchor: lights were on at 06:00
(ZT0) in a 12:12 light/dark cycle, with telemetry sampled every 5 minutes. The
source figure's black bars mark the 18:00-06:00 dark phases.

`[1] target.png` is the experimental temperature/E2 reference. The modeled E2
figure uses its panel-B alignment: 08:00 to 08:00 with matching time-of-day and
E2 axis bounds.

### Editing and rerunning the model

1. Change a parameter or equation in `MODEL.py` and save it.
2. Open a terminal in this `meeting6` folder.
3. Run `python plot_model.py`.

`plot_model.py` imports the saved `MODEL.py`, calls `run()`, and overwrites the
three final PNG files. Edit `MODEL.py` to change the biology; edit
`plot_model.py` only to change colors, labels, layout, or output filenames.

The main experimental controls near the top of `MODEL.py` are:

- `ENDOCRINE_CLOCK_PEAK_H` and `CLOCK_GAIN`: timing and strength of endocrine
  pulse modulation.
- `THERMAL_CLOCK_PEAK_H`: timing of the separate 0–1 thermal drive.
- The conservation and loss weights: how the clock and hormones enter each
  thermal pathway.
- `BASAL_HEAT_RATE`, `CONSERVATION_GAIN`, and `HORMONAL_LOSS_GAIN`: terms in
  the temperature-balance equation.
- `SELF_LOSS_MAX` and `SELF_LOSS_HALF_C`: strength and half-saturation of
  direct temperature-dependent heat loss.

Reference material is kept under `references/`; it is not consumed by the
solver. `thermoregulation_feedback_solid_arrows.png` is the final model
schematic.

## Biological references for self-regulated loss

- [Cabanac and Massonnet, J Physiol (1977)](https://pubmed.ncbi.nlm.nih.gov/856985/)
- [Fujii et al., J Appl Physiol (2014)](https://pubmed.ncbi.nlm.nih.gov/25257867/)
- [Minson et al., J Appl Physiol (2002)](https://pubmed.ncbi.nlm.nih.gov/12381749/)
- [Abbott and Saper, J Physiol (2017)](https://pubmed.ncbi.nlm.nih.gov/28786483/)
- [Krajewski-Hall et al., Temperature (2017)](https://pubmed.ncbi.nlm.nih.gov/29687044/)
- [Skop et al., Am J Physiol Endocrinol Metab (2020)](https://pubmed.ncbi.nlm.nih.gov/32691633/)
