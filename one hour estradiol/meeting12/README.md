# Meeting 12 — ultradian parameter sweep and FFT comparison

Run `python meeting12/parameter_sweep.py` from the `one hour estradiol` directory. The explorer opens **one window** with E2/LH time series on the left, the full E2 spectrum and four zooms on the right, and parameter controls below. Both E2 and FFT panels place the gray digitized Licinio reference behind the live model.

## Current baseline

The defaults follow the corrected newer baseline in **Ultradian Rhythm Parameters** and the existing `model.py`. `sweep_model.py` provides the same equations with editable parameters and an optional Windows speed helper. Reset restores these values:

| Parameter | Default |
|---|---:|
| Basal E2 production | 6 |
| Responsive E2 production | 210.7 |
| LH activation exponent | 5 |
| LH desensitization exponent | 5 |
| E2-to-N midpoint (pg/mL) | 30.2 |
| Central Hill exponent | 2.36 |
| Fast period coefficient (min) | 67.5 |
| Slow period coefficient (min) | 76.73 |
| Sleep period multiplier | 1.55 |
| SCN minimum / maximum | 35 / 277 |
| SCN peak clock time | 23:48 (solver phase 16.8 h from 07:00) |
| GnRH / LH / E2 half-lives (min) | 4.09 / 29 / 22.53 |

The six unit coefficients remain fixed at 1, and ovarian responsiveness has a one-hour response time. Initial E2 is 19 pg/mL, matching `model.py`. Both the Python solver and Windows helper use start-of-step forward Euler updates at 1.5 seconds.

The chat's SCN minimum 25, basal production 5, desensitization exponent 8, and E2-to-N exponent 3.5 were separate experimental variants, not simultaneous baseline changes. The baseline retains the shared central exponent of 2.36 used by `model.py`. Isolating the E2-to-N-only variant would require a separate feedback exponent; the central slider currently changes both feedback and N's pulse maps.

## Interpretation and controls

- One parameter at a time resets other parameters to this baseline when a value is changed. Uncheck it to explore joint changes.
- The dotted slider marker indicates the current default. Fine scale shows a smaller exploration range; Full range allows larger changes. Lo/Hi select the visible slider limits.
- PASS/FAIL and the Failures tab evaluate the displayed three days after the 96-hour warmup.
- The Robustness tab explicitly labels the saved eight-day ranges as **historical**, with their earlier defaults. Those ranges have not been recomputed for this baseline and are not shown as current passing bands.

## Licinio and FFT calculation

`licinio.csv` contains a digitized day beginning at 08:00. It is repeated across the three-day comparison and aligned to the time-series clock, which starts at 07:00. These are repeated digitized measurements, not three independently recorded days.

The FFT follows `compare.py`: both model and reference are sampled every seven minutes over 72 hours from 08:00, linearly detrended, and transformed with `rfft`. Power is `abs(coefficients)**2 / N`; zero frequency and the final FFT bin are omitted. The plots retain raw power without normalization to each trace's maximum. Views cover the full spectrum and periods near 1, 2, 4–6, and 8–12 hours.

## Verification

```shell
python -B meeting12/check_parameter_sweep_fft.py
```

This checks agreement with `model.py`, Python/Windows solver parity, the SciPy FFT convention, single-window rendering, clock alignment, slider updates, reset, tabs, endpoint controls, and plot limits. It saves the combined preview to `parameter_sweep.png`. NumPy and Matplotlib are required by the explorer; SciPy is also used by the independent verification and `compare.py`.
