# Pre meeting 5 — general analysis files

This folder contains the general models, analyses, figures, and historical handoff notes previously stored at the project root. The meeting folders and literature PDFs remain one level up.

- **`../meeting12/`**: current parameter sweep, updated ultradian baseline, digitized Licinio reference, and FFT comparisons.
- **`../meeting11/`**: earlier endocrine model, equations, and saved plots.
- **`../meeting_5/` through `../meeting10/`**: earlier meeting work and supporting figures/data.
- **`../`**: literature PDFs and meeting folders. Older handoff notes in this folder describe earlier stages of the work; start with `../meeting12/README.md` for the current model.

## Run the interactive comparison

From this directory:

```shell
python ../meeting12/parameter_sweep.py
```

This opens one window containing E2/LH time series, parameter controls, and live E2 FFT plots. Gray lines show the digitized Licinio day repeated across three days. Both spectra use seven-minute samples, linear detrending, and raw FFT power. See `../meeting12/README.md` for details.

NumPy and Matplotlib are required. The independent FFT verification and `../meeting12/compare.py` also use SciPy.

```shell
python -m pip install numpy matplotlib scipy
python -B ../meeting12/check_parameter_sweep_fft.py
```

The verification checks solver agreement, the current baseline, FFT convention, clock alignment, slider updates, and plot limits, then refreshes `../meeting12/parameter_sweep.png`.

Local temporary, archive, audit, and older thermoregulation development folders are excluded from Git. They were moved out of the active project during the November 3, 2026 cleanup.
