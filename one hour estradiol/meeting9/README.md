# Meeting 9 minimal HPG model

This folder begins from the Meeting 8 endocrine framework while removing the
separate direct pituitary E2-to-LH inhibition term. The retained feedback loop
is

$$
E_2 \;\dashv\; N \longrightarrow \mathrm{GnRH}
\longrightarrow \mathrm{LH} \longrightarrow E_2 .
$$

The five dynamic quantities are aggregate KNDy/kisspeptin drive $N$, pulse
counter $\phi$, GnRH, LH, and E2. Time $t$ is in hours. The equations below
show the numerical defaults in `MODEL.py`; concentrations and state values are
understood in the units listed in the constants table.

## Complete HPG equations at the default settings

### Feedback and pulse-response functions

The E2-to-$N$ response is decreasing, whereas the two $N$ responses are
increasing. All three Hill exponents are 5:

$$
F_N(E_2)
= \frac{11^5}{11^5 + E_2^5},
$$

$$
H_P(N)
= \frac{N^5}{1^5 + N^5},
\qquad
H_A(N)
= \frac{N^5}{1^5 + N^5}.
$$

$H_P$ controls pulse period and $H_A$ controls pulse amplitude. They have the
same formula under the current defaults, but are written separately because
they play different roles in the model.

### Continuous dynamics between pulses

Between GnRH pulse events, the HPG state obeys

$$
\begin{aligned}
\frac{dN}{dt}
  &= 1\left[100\,\frac{11^5}{11^5 + E_2^5} - N\right], \\
\frac{d\,\mathrm{GnRH}}{dt}
  &= -10.397207708\,\mathrm{GnRH}, \\
\frac{d\,\mathrm{LH}}{dt}
  &= 1\,\mathrm{GnRH} - 2.079441542\,\mathrm{LH}, \\
\frac{dE_2}{dt}
  &= 15\,\mathrm{LH} - 1.890401402\,E_2.
\end{aligned}
$$

The three decay rates come directly from the default half-lives:

$$
\begin{aligned}
k_{\mathrm{GnRH}}
  &= \frac{\ln 2}{4/60\ \mathrm{h}}
   = 10.397207708\ \mathrm{h}^{-1}, \\
k_{\mathrm{LH}}
  &= \frac{\ln 2}{20/60\ \mathrm{h}}
   = 2.079441542\ \mathrm{h}^{-1}, \\
k_{E_2}
  &= \frac{\ln 2}{22/60\ \mathrm{h}}
   = 1.890401402\ \mathrm{h}^{-1}.
\end{aligned}
$$

There is no direct $E_2\to\mathrm{LH}$ multiplier in the Meeting 9 model.

### Pulse timing and instantaneous GnRH jump

The current $N$ value sets a pulse period between 75 and 51 minutes:

$$
P(N)
= 75 - (75-51)H_P(N)
= 75 - 24\,\frac{N^5}{1^5+N^5}
\quad \text{minutes}.
$$

The pulse counter accumulates cycles according to

$$
\frac{d\phi}{dt} = \frac{60}{P(N)}.
$$

A pulse occurs whenever $\phi$ crosses an integer. At each event time $t_j$,
GnRH receives the instantaneous jump

$$
\mathrm{GnRH}(t_j^+)
= \mathrm{GnRH}(t_j^-)
+ 100\,
\frac{\left[N(t_j^-)\right]^5}
     {1^5+\left[N(t_j^-)\right]^5}.
$$

$N$, LH, and E2 remain continuous across the event. Thus the pulse gain is
100, with no additive pulse-amplitude floor and no random or circadian pulse
term.

### Initial conditions and integration window

The simulation starts at the beginning of the discarded warm-up with

$$
\begin{aligned}
\phi(0) &= 0, &
\mathrm{GnRH}(0) &= 0, &
\mathrm{LH}(0) &= 5, &
E_2(0) &= 15, \\
N(0) &= 100\,\frac{11^5}{11^5+15^5}
     = 17.4974414.
\end{aligned}
$$

`MODEL.py` integrates these equations by explicit Euler with
$\Delta t=1/600$ hour (6 seconds), discards a 96-hour warm-up, and returns the
following 72 hours.

## Default constants

| Quantity | Model name | Default |
|---|---|---:|
| $N$ target ceiling | `N_GAIN` | 100 |
| $N$ turnover rate | fixed | 1 h$^{-1}$ |
| E2-to-$N$ midpoint | `N_E2_HALF_PG_ML` | 11 pg/mL |
| E2-to-$N$ Hill exponent | fixed | 5 |
| $N$-to-period midpoint | fixed | 1 $N$ unit |
| $N$-to-period Hill exponent | fixed | 5 |
| Minimum / maximum pulse period | fixed | 51 / 75 min |
| $N$-to-amplitude midpoint | fixed | 1 $N$ unit |
| $N$-to-amplitude Hill exponent | fixed | 5 |
| GnRH pulse gain | `GNRH_PULSE_MAX_JUMP` | 100 |
| GnRH-to-LH gain | `LH_GNRH_GAIN` | 1 |
| LH-to-E2 gain | `E2_GAIN` | 15 |
| GnRH half-life / decay rate | fixed | 4 min / 10.397207708 h$^{-1}$ |
| LH half-life / decay rate | fixed | 20 min / 2.079441542 h$^{-1}$ |
| E2 half-life / decay rate | fixed | 22 min / 1.890401402 h$^{-1}$ |

## Calibration and interpretation

The module-level calibration constants in `MODEL.py` are phenomenological:
the `N` production gain/target ceiling, E2-to-`N` midpoint, the GnRH pulse
gain, GnRH-to-LH gain, and LH-to-E2 gain. `N` turnover and both `N` response
midpoints are fixed at 1. The former Meeting 8 factors `102.4` and
`3.5` have been combined into one GnRH pulse-gain parameter because only their
product affects the equations. `N` gain and the two `N`-response midpoints are
partly scale-redundant, so they should not be interpreted as independently
identified biological quantities.

The interactive explorer reports the retained E2 range and robust
peak-height/peak-spacing variation. A setting passes the range constraint only
when its complete retained trajectory stays between 7 and 60 pg/mL.
It also plots the retained LH state against a 1--10 IU/L target band. Because
the LH state is now interpreted in IU/L, parameter sets should be accepted only
when both complete trajectories remain inside their respective target bands.

The complete retained trajectory is E2 11.39--45.84 pg/mL and LH 0.50--8.07
IU/L. The E2 trajectory is inside the Licinio target range of 7--60 pg/mL;
the LH minimum remains below the separate 1--10 IU/L target band. The defaults
are the chosen high-period irregular calibration rather than a claim that both
constraints have been exactly met.

The raw GnRH pulse, GnRH-to-LH, and LH-to-E2 gains are multiplicatively
coupled, so their individual effects are not identifiable from E2 alone. Each
gain remains directly adjustable in the explorer, without a responsiveness
term or helper function. Small changes can shift discrete pulse times; those
phase shifts then propagate through the nonlinear feedback loop. Consequently,
even a slider with a clear immediate effect need not change the 72-hour extrema
or irregularity monotonically.

## Python slider program

`parameter_sweep.py` is the non-HTML version of the explorer. It uses Python
and Matplotlib, reads its defaults and fixed constants from `MODEL.py`, and
contains five whole-number sliders for the non-literature-backed constants.
It shows live E2 and LH time series, both target bands, and separate E2 and LH
peak-to-peak interval plots. The interval-plot horizontal axis counts
successive pairs of peaks; the vertical axis is the number of minutes between
the two peaks in each pair.

From the `meeting9` folder, launch it with:

```powershell
python parameter_sweep.py
```

If Python reports that a plotting package is missing, install the two required
packages once:

```powershell
python -m pip install numpy matplotlib
```

Dragging is kept responsive with a 30-second-step moving preview. Releasing
the mouse automatically replaces that preview with the exact six-second-step
simulation used by `MODEL.py`; the title drops the words `moving preview` when
the exact redraw is ready. Only one point every three minutes is sent to the
plots, but every six-second value is included in the exact minimum and maximum
range checks. Use the **Reset defaults** button to return to the chosen
calibration above.

`parameter_explorer.html` remains available as a separate browser version.

## Complete core-temperature equations

Temperature is driven by the endocrine E2 trajectory but does not feed back
into the HPG subsystem. Let $T$ be core temperature, $T_a=21.25\ ^\circ
\mathrm{C}$ be ambient temperature, $B=0.32\ \mathrm{kcal/h}$ be OVX
background heat, $C_b$ be whole-body heat capacity, and $G$ be whole-body
thermal conductance.

### Heat capacity and conductance

The heat-capacity calculation uses the Zengin short-term WT-OVX body mass
$M=22.36\ \mathrm{g}$ and fat mass $M_f=3.98\ \mathrm{g}$:

$$
f_{\mathrm{fat}}=\frac{M_f}{M}=\frac{3.98}{22.36}=0.1779964,
$$

$$
c_b=0.9100-0.5231f_{\mathrm{fat}}
=0.8168901\ \frac{\mathrm{cal}}{\mathrm{g}\,{}^\circ\mathrm{C}},
$$

$$
C_b=\frac{c_bM}{1000}
=0.018265662\ \frac{\mathrm{kcal}}{{}^\circ\mathrm{C}}.
$$

The Abreu-Vieira conductance anchor is scaled from a $26.97$-g mouse to the
same $22.36$-g reference mass using approximate surface-area scaling:

$$
G=0.0255\left(\frac{22.36}{26.97}\right)^{2/3}
=0.022504399\ \frac{\mathrm{kcal}}{\mathrm{h}\,{}^\circ\mathrm{C}}.
$$

The corresponding thermal rate constant is

$$
\kappa=\frac{G}{C_b}=1.232060399\ \mathrm{h}^{-1}.
$$

### Shared SCN-E2 heat function

The thermal calibration uses a maximum shared E2-SCN thermogenic fraction of
`0.14`, Hill midpoint `K = 30 pg/mL`, Hill exponent `n = 5`, and SCN
E2-equivalent scale `k_SCN = 11 pg/mL`. With
$ZT(t)=t\bmod 24$, the nonnegative SCN-associated BAT drive, effective shared
input, Hill response, and heat input are

$$
S(t)=\frac{1-\cos\!\left(2\pi[ZT(t)-7]/24\right)}{2},
\qquad
X(t)=E_2(t)+11S(t),
$$

$$
H_{E_2+SCN}(t)=\frac{X(t)^5}{30^5+X(t)^5},
\qquad
Q_{E_2+SCN}(t)=0.14(0.32)H_{E_2+SCN}(t)
=0.0448H_{E_2+SCN}(t)\ \mathrm{kcal/h}.
$$

At `E2 = 39 pg/mL` and `SCN = 0`, this term contributes approximately
`0.11 B`, matching the intended Camporez-scale effect while leaving the gate
responsive to the rhythmic SCN-equivalent input at lower cycling E2 levels.

There is no separate additive SCN term in the Meeting 9 temperature equation;
the SCN signal acts only through $X(t)$ and the shared Hill response above.

### Activity input

Activity is represented by a five-hour rectangular dark-phase bout:

$$
A(t)=
\begin{cases}
1, & 12\le ZT(t)<17,\\
0, & \text{otherwise}.
\end{cases}
$$

Its temperature-equivalent scale is
$\Delta T_A=1.9\ ^\circ\mathrm{C}$. The multiplier $m_A$ is `1.0` in the
final Meeting 9 run, but remains an argument of the temperature-rate function.

### Temperature balance and Euler update

The complete heat-flow balance is

$$
C_b\frac{dT}{dt}
=B+0.14B\,H_{E_2+SCN}(t)
+Gm_A\Delta T_A A(t)
-G\left[T(t)-T_a\right].
$$

Equivalently, after division by $C_b$,

$$
\frac{dT}{dt}
=\frac{B}{C_b}
+\frac{0.14B}{C_b}H_{E_2+SCN}(t)
+\kappa m_A\Delta T_A A(t)
-\kappa\left[T(t)-T_a\right].
$$

Substituting the Meeting 9 defaults gives

$$
\frac{dT}{dt}
=17.519211732
+2.452689642\,H_{E_2+SCN}(t)
+2.340914758\,m_A A(t)
-1.232060399\left[T(t)-21.25\right]
\quad {}^\circ\mathrm{C}/\mathrm{h}.
$$

The simulation starts at $T(0)=37.0\ ^\circ\mathrm{C}$ at the beginning of
the 96-hour warm-up. On the shared six-second grid
$\Delta t=1/600\ \mathrm{h}$, `MODEL.py` applies explicit Euler:

$$
T_{i+1}=T_i+\Delta t\,
f_T\!\left(T_i,E_{2,i},S_i,A_i\right),
$$

where $f_T$ is the right-hand side of the temperature-rate equation above.
The final 72 hours are retained after the warm-up.

### Literature grounding for `k_SCN = 11 pg/mL`

`k_SCN` is a phenomenological E2-equivalent scale, not a measured biochemical
conversion between SCN activity and estradiol concentration. We calibrate it to
the approximately `0.38 °C` activity-independent day-night body-temperature
contrast reported in rats by [Scheer et al. (2005)](https://doi.org/10.1016/j.neuroscience.2004.12.012),
whose SCN-lesion experiments also support an SCN origin for the endogenous
circadian temperature rhythm.

For the normalized SCN waveform used here, its 12-hour light- and dark-phase
means are determined by the waveform geometry:

$$
\overline{S}_{\mathrm{light}}=0.1925,
\qquad
\overline{S}_{\mathrm{dark}}=0.8075.
$$

With thermal conductance `G = 0.0225044 kcal/(h °C)`, a `0.38 °C` contrast
requires a heat difference of `G(0.38) = 0.00855 kcal/h`. Because the maximum
shared E2-SCN heat term is `0.14(0.32) = 0.0448 kcal/h`, the required change in
the phase-average Hill output is approximately

$$
\Delta H_{\mathrm{required}}
=\frac{0.00855}{0.0448}
=0.1909.
$$

At a representative cycling value of `E2 = 30 pg/mL`, retaining
`k_SCN = 11 pg/mL` gives phase-average effective inputs of `32.12` and
`38.88 pg/mL`. The corresponding Hill outputs are `0.5845` and `0.7853`, so
`ΔH = 0.2008`. In the simple heat-balance conversion, this is a `0.400 °C`
contrast, close to the literature target of `0.38 °C`. We therefore retain the
rounded constant `k_SCN = 11 pg/mL`. The values `0.1925` and `0.8075` above are
mathematical averages of the chosen SCN waveform, not direct literature
measurements.

## Run

```python
import MODEL

result = MODEL.run()
e2 = result["E2_pg_ml"]
lh = result["LH_IU_L"]
temperature = result["temperature_c"]
```

The finalized `MODEL.run()` uses the module-level constants listed above.
