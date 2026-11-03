# Estradiol model project

The active project contains the meeting folders and the general analysis files at this level.

- **`meeting11/`**: current endocrine model, interactive parameter sweep, equations, and saved plots.
- **`meeting12/`**: digitized Licinio reference day and the fixed-model FFT comparison.
- **`meeting_5/` through `meeting10/`**: earlier meeting work and supporting figures/data.
- **Top-level files**: general models, analyses, figures, and reference material. Older handoff notes describe earlier stages of the work; start with `meeting11/README.md` for the current model.

## Run the interactive comparison

From this directory:

```shell
python meeting11/parameter_sweep.py
```

This opens the parameter controls and a live E2 FFT power window. Gray lines show the digitized Licinio day repeated across three days. Both spectra use seven-minute samples, linear detrending, and raw FFT power. See `meeting11/README.md` for details.

NumPy and Matplotlib are required. The independent FFT verification and `meeting12/compare.py` also use SciPy.

```shell
python -m pip install numpy matplotlib scipy
python -B meeting11/check_parameter_sweep_fft.py
```

The verification checks the FFT convention, clock alignment, slider updates, and plot limits, then refreshes `meeting11/parameter_sweep_licinio.png` and `meeting11/parameter_sweep_fft.png`.

Local temporary, archive, audit, and older thermoregulation development folders are excluded from Git. They were moved out of the active project during the November 3, 2026 cleanup.
