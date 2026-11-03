# Meeting 11 — final endocrine model

The rounded reference model has a **fixed one-hour ovarian response time**. `MODEL.py` defines the model and reference constants; `parameter_sweep.py` imports that same model.

`MODEL.py` keeps the Python equations and a small Windows speed helper for responsive slider updates. The helper runs only this model and returns LH/E2; the Python solver returns all states and also provides the fallback when the helper is unavailable.

Run `python parameter_sweep.py` to open the parameter controls and a companion E2 FFT power window. The E2 time series and FFT plots show the digitized Licinio Figure 2 day from `../meeting12/licinio.csv` in gray behind the live model. The reference day starts at 08:00 and repeats over the three-day comparison; the time-series view retains its 07:00 start. The spectrum updates after each parameter change and includes a full view plus zooms near 1, 2, 4–6, and 8–12 hours.

FFT processing matches `meeting12/compare.py`: both traces are sampled every seven minutes for 72 hours starting at 08:00 after the model warmup, linearly detrended, and transformed using `rfft`. Power is `abs(coefficients)**2 / N`, with the zero-frequency and final bins omitted. The gray reference is a repeated digitized day, not three independently recorded days; the plots show raw FFT power without normalization to each trace's maximum.

## Notation and clock

Time $t$ is in **hours from 07:00**, and $z(t)=t\bmod24$.

For the SCN equation, use **clock time in decimal hours**:

$$
c(t)=(t+7)\bmod24.
$$

Thus 07:00 is $c=7$, and 23:30 is $c=23.5$. The circadian phase $\theta_{\mathrm{clock}}$ is expressed on this same clock.

| Symbol | Meaning | Units |
|:--|:--|:--|
| $N$ | Central signal inhibited by E2 | Dimensionless |
| $G$ | GnRH state | Model units |
| $L$ | LH concentration | IU/L |
| $R$ | Ovarian responsiveness | Dimensionless |
| $E$ | Estradiol (E2) concentration | pg/mL |
| $\phi$ | Accumulated pulse phase | Cycles |

Sleep is 23:00–07:00, or model hours 16–24:

$$
S(t)=
\begin{cases}
1, & 16\le z(t)<24,\\
0, & 0\le z(t)<16.
\end{cases}
$$

## 1. Equations with symbolic constants

### Central E2 feedback

$$
\boxed{
\frac{dN}{dt}
=
\frac{K_N^h}{K_N^h+E^h}-N
}
$$

### Pulse timing and circadian amplitude

$$
P_0(N)=P_{\max}-(P_{\max}-P_{\min})\frac{N^h}{1+N^h}.
$$

$$
\boxed{
\frac{d\phi}{dt}
=
\frac{60}{P_0(N)\left[1+(m_s-1)S(t)\right]}
}
$$

The period coefficients are in minutes; 60 converts the phase rate to cycles per hour.

$$
M_{\mathrm{SCN}}(t)
=
\frac{M_{\max}+M_{\min}}{2}
+
\frac{M_{\max}-M_{\min}}{2}
\cos\!\left(\frac{2\pi[c(t)-\theta_{\mathrm{clock}}]}{24}\right).
$$

A GnRH pulse occurs at each integer phase crossing. Between pulses,

$$
\boxed{\frac{dG}{dt}=-k_GG}
$$

and each crossing adds

$$
\boxed{
\Delta G=\frac{N^h}{1+N^h}\,M_{\mathrm{SCN}}(t).
}
$$

Phase accumulates without resetting after a pulse.

### LH, ovarian responsiveness, and E2

$$
\boxed{\frac{dL}{dt}=G-k_LL}
$$

$$
\boxed{
\frac{dR}{dt}
=
\frac{1}{1+L^m}-R
}
$$

$$
\boxed{
\frac{dE}{dt}
=
b+VR\frac{L^a}{1+L^a}-k_EE
}
$$

LH stimulates E2 production through the increasing Hill term and lowers the target for ovarian responsiveness through the decreasing Hill term. Responsiveness approaches that target with a fixed one-hour response time.

For hormone half-lives $T_G,T_L,T_E$ in **minutes**, the hourly clearance rates are

$$
k_G=\frac{60\ln2}{T_G},
\qquad
k_L=\frac{60\ln2}{T_L},
\qquad
k_E=\frac{60\ln2}{T_E}.
$$

## 2. Equations with the final constants inserted

### Central feedback

$$
\boxed{
\frac{dN}{dt}=\frac{50^3}{50^3+E^3}-N
}
$$

### Pulse timing

$$
\boxed{
\frac{d\phi}{dt}
=
\frac{60}
{\left(75-10\frac{N^3}{1+N^3}\right)\left[1+0.57S(t)\right]}
}
$$

### GnRH decay and pulse size

Between pulses,

$$
\boxed{
\frac{dG}{dt}
=
-\frac{60\ln2}{4}\,G
=
-15\ln2\,G
}
$$

At each integer phase crossing,

$$
\boxed{
\Delta G
=
\frac{N^3}{1+N^3}
\left[
150+100\cos\!\left(\frac{2\pi[c(t)-23.5]}{24}\right)
\right].
}
$$

The SCN equation itself is

$$
\boxed{
M_{\mathrm{SCN}}(t)
=
150+100\cos\!\left(\frac{2\pi[c(t)-23.5]}{24}\right).
}
$$

The circadian pulse-amplitude factor ranges from **50 to 250** and peaks at **23:30**.

### LH

$$
\boxed{
\frac{dL}{dt}=G-\frac{60\ln2}{29}\,L
}
$$

### Ovarian responsiveness

$$
\boxed{
\frac{dR}{dt}=\frac{1}{1+L^4}-R
}
$$

### Estradiol

$$
\boxed{
\frac{dE}{dt}
=
20+250R\frac{L^4}{1+L^4}-\frac{60\ln2}{22.5}\,E
}
$$

The E2 clearance coefficient is exactly $\frac{8}{3}\ln2$ per hour.

## Parameters and one-at-a-time robustness

These are the visualization's approximate **one-parameter-at-a-time passing bands**: one parameter varies while all other parameters stay at their defaults. The percentages describe parameter changes, not a probability of passing.

| Visualization tab | Parameter | Default | Approximate passing range | Change from default |
|:--|:--|--:|:--|:--|
| E2 Constants | Basal E2 production $b$ (pg/mL/h) | 20 | 19.10625–20.80625 | -4.47% to +4.03% |
|  | LH-driven E2 production $V$ (pg/mL/h) | 250 | 245.78125–253.515625 | -1.69% to +1.41% |
| Hypothalamus | E2-to-N midpoint $K_N$ (pg/mL) | 50 | 49.46875–50.640625 | -1.06% to +1.28% |
|  | Minimum SCN drive $M_{\min}$ | 50 | 38.03125–55.046875 | -23.94% to +10.09% |
|  | Maximum SCN drive $M_{\max}$ | 250 | 243.046875–258.359375 | -2.78% to +3.34% |
| Literature grounded | Fast period coefficient $P_{\min}$ (min) | 65 | 60.8765625–68.7375 | -6.34% to +5.75% |
|  | Slow period coefficient $P_{\max}$ (min) | 75 | 73.2421875–76.40625 | -2.34% to +1.88% |
|  | Sleep period multiplier $m_s$ | 1.57 | 1.513578125–1.605815625 | -3.59% to +2.28% |
|  | SCN peak time | 23:30 | 22:02–00:08 (next day) | -8.88% to +3.88%* |
| Steepness | Central Hill exponent $h$ | 3 | 2.775–3.135 | -7.50% to +4.50% |
|  | LH activation exponent $a$ | 4 | 3.15375–5.59 | -21.16% to +39.75% |
|  | LH inhibition exponent $m$ | 4 | 3.55375–4.43625 | -11.16% to +10.91% |
| Half-lives | GnRH half-life $T_G$ (min) | 4 | 3.91–4.09125 | -2.25% to +2.28% |
|  | LH half-life $T_L$ (min) | 29 | 28.565–29.4621875 | -1.50% to +1.59% |
|  | E2 half-life $T_E$ (min) | 22.5 | 22.2328125–22.73203125 | -1.19% to +1.03% |

*SCN percentages use the solver's 16.5-hour phase reference; the clock-time window is equivalent to **87.86 minutes earlier through 38.36 minutes later** than 23:30. Phase percentages depend on the chosen time origin.

The saved bands were assessed over eight consecutive days after both 96-hour and 192-hour warmups at both 3-second and 1.5-second steps: **32 daily checks**. They describe sampled passing settings, rather than a guarantee for every intermediate value or for simultaneous parameter changes.

Although $P_{\min}=65$ minutes is the strong-$N$ asymptote, this reference has $0\le N\le1$, so its waking phase period is **70–75 minutes**. The sleep multiplier applies to that period.

## Initial conditions and numerical rule

$$
\begin{aligned}
G(0)&=0, & L(0)&=3, & E(0)&=15, & \phi(0)&=0,\\
N(0)&=\frac{50^3}{50^3+15^3}, &&
R(0)&=\frac{1}{1+3^4}=\frac{1}{82}.
\end{aligned}
$$

The solver uses forward Euler updates evaluated from the start of each step. The GnRH increment uses the start-of-step $N$ and circadian factor, multiplied by

$$
\left\lfloor\phi_{i+1}\right\rfloor-\left\lfloor\phi_i\right\rfloor.
$$

The explorer always shows the first three days after the normal 96-hour warmup, using 1.5-second steps. Its PASS/FAIL status identifies failing days, and the Failures tab lists each failed criterion for those three days. The green slider bands retain the existing eight-day, one-parameter-at-a-time results described above; they do not predict the outcome of simultaneous changes.

The SCN peak control accepts decimal clock hours or HH:MM: enter **23.5** or **23:30** for the reference peak. Numeric inputs, Lo/Hi endpoint buttons, the Robustness and Failures tabs, and both mode switches remain available.

## Run

From the project directory:

~~~powershell
python meeting11/MODEL.py
python meeting11/parameter_sweep.py
~~~

Running `parameter_sweep.py` opens the window directly, with no command-line options. The model needs NumPy; the explorer also needs Matplotlib. Open this README in a Markdown preview (Ctrl+Shift+V in VS Code) to render the equations.
