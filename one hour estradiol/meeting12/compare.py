"""Run this file to show the fixed model and Licinio as overlaid lines."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import detrend

from model import simulate_e2

folder = Path(__file__).resolve().parent

# 1. Sample both traces every 7 minutes for 3 days, starting at 08:00.
step = 7 / 60  # hours
hours = np.arange(0, 72, step)
model_time, e2 = simulate_e2()
model = np.interp(96 + 1 + hours, model_time, e2)  # 96 h warmup, then 07:00 to 08:00

day_hours, day_e2 = np.loadtxt(folder / 'licinio.csv', delimiter=',', skiprows=1, unpack=True)
licinio = np.interp(hours % 24, day_hours, day_e2, period=24)  # Repeat the digitized day

# 2. Ivy's FFT calculation. Drop the zero-frequency and final Nyquist bins.
def power(trace):
    coefficients = np.fft.rfft(detrend(trace))
    return np.abs(coefficients[1:-1])**2 / len(trace)


model_power = power(model)
licinio_power = power(licinio)
frequency = np.fft.rfftfreq(len(hours), d=step)[1:-1]  # cycles per hour
period = 1 / frequency                                 # hours per cycle

# 3. Plot the time series above the full FFT and its four zooms.
fig, axes = plt.subplot_mosaic(
    [['time', 'time'], ['full', 'full'], ['1 h', '2 h'], ['4–6 h', '8–12 h']],
    figsize=(11, 11), constrained_layout=True
)
orange, blue = '#b45c25', '#156ba5'
ax = axes['time']
ax.plot(hours, licinio, color=orange, label='Licinio · digitized day repeated')
ax.plot(hours, model, color=blue, label='Fixed model')
for start in [15, 39, 63]:  # 23:00 to 07:00 each night
    ax.axvspan(start, start + 8, color='lightgray', alpha=0.3)
ax.set(xlim=(0, 72), xlabel='Hours from 08:00', ylabel='E2 (pg/mL)')
ax.legend()

windows = {'full': (0.5, 36), '1 h': (0.9, 1.45), '2 h': (1.65, 2.3),
           '4–6 h': (3.5, 6.5), '8–12 h': (7, 13)}
for name, (low, high) in windows.items():
    ax = axes[name]
    keep = (period >= low) & (period <= high)
    ax.plot(period[keep], licinio_power[keep], '.-', color=orange)
    ax.plot(period[keep], model_power[keep], '.-', color=blue)
    ax.set(xlim=(low, high), ylim=(0, None), xlabel='Period (hours)', ylabel='FFT power',
           title='Full spectrum' if name == 'full' else f'Zoom · {name}')

axes['full'].set_xscale('log')
ticks = [0.5, 1, 2, 4, 8, 12, 24, 36]
axes['full'].set_xticks(ticks, labels=[str(t) for t in ticks])
fig.savefig(folder / 'e2_FFTs.png', dpi=160)
plt.show()
