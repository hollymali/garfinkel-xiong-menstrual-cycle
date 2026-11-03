# Estradiol model project

The active project contains the meeting folders and the general analysis files at this level.

- **`meeting12/`**: current parameter sweep, updated ultradian baseline, digitized Licinio reference, and FFT comparisons.
- **`meeting11/`**: earlier endocrine model, equations, and saved plots.
- **`meeting_5/` through `meeting10/`**: earlier meeting work and supporting figures/data.
- **Top-level files**: general models, analyses, figures, and reference material. Older handoff notes describe earlier stages of the work; start with `meeting12/README.md` for the current model.

## Run the interactive comparison

From this directory:

```shell
python meeting12/parameter_sweep.py
```

This opens one window containing E2/LH time series, parameter controls, and live E2 FFT plots. Gray lines show the digitized Licinio day repeated across three days. Both spectra use seven-minute samples, linear detrending, and raw FFT power. See `meeting12/README.md` for details.

NumPy and Matplotlib are required. The independent FFT verification and `meeting12/compare.py` also use SciPy.

```shell
python -m pip install numpy matplotlib scipy
python -B meeting12/check_parameter_sweep_fft.py
```

The verification checks solver agreement, the current baseline, FFT convention, clock alignment, slider updates, and plot limits, then refreshes `meeting12/parameter_sweep.png`.

Local temporary, archive, audit, and older thermoregulation development folders are excluded from Git. They were moved out of the active project during the November 3, 2026 cleanup.
