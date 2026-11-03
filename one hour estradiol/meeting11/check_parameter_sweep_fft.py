"""Check FFT compatibility, clock alignment, and live parameter updates headlessly."""
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import detrend

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "meeting11"))
import parameter_sweep as sweep

# Compare the actual numerical convention with the existing SciPy calculation.
day_hours, day_e2 = sweep.load_licinio()
reference = np.interp(sweep.FFT_HOURS, day_hours, day_e2, period=24)
signals = [reference, np.ones(len(reference)),
           20 + .2 * sweep.FFT_HOURS + 5 * np.sin(2 * np.pi * sweep.FFT_HOURS)]
for signal in signals:
    expected = np.abs(np.fft.rfft(detrend(signal))[1:-1]) ** 2 / len(signal)
    np.testing.assert_allclose(sweep.fft_power(signal), expected, rtol=1e-10, atol=1e-10)
period, power = sweep.fft_spectrum(signals[-1])
assert abs(period[np.argmax(power)] - 1) < .02
print("FFT agrees with meeting12's SciPy detrend and detects the 1-hour signal.", flush=True)

explorer = sweep.ParameterExplorer(plt)
explorer.timer.stop()

def finish():
    deadline = time.monotonic() + 90
    while explorer.pending is not None or explorer.future is not None:
        assert time.monotonic() < deadline, "Live plot update timed out"
        explorer.poll()
        time.sleep(.01)
    assert explorer.result is not None, explorer.status.get_text()
    assert explorer.result["parameters"] == explorer.parameters

try:
    finish()
    # The digitized zero-hour point belongs at 08:00, one hour into this view.
    x, y = explorer.licinio_line.get_data()
    for day in range(3):
        for index in (0, 60, 120):
            position = np.flatnonzero(np.isclose(x, 1 + 24 * day + day_hours[index]))
            if len(position):
                np.testing.assert_allclose(y[position[0]], day_e2[index])
    assert (x[0], x[-1]) == (0., 72.)
    assert explorer.axes[0].get_ylim() == (0., 60.)
    assert explorer.result["trace"]["time_h"][-1] >= sweep.VIEW_START + 73
    reference_lines = {name: axis.lines[0].get_ydata().copy()
                       for name, axis in explorer.fft_axes.items()}
    original = explorer.fft_lines["full"].get_ydata().copy()
    explorer.sliders["constant_inflow"].set_val(20.5)
    finish()
    assert not np.allclose(original, explorer.fft_lines["full"].get_ydata())
    for name, axis in explorer.fft_axes.items():
        np.testing.assert_array_equal(reference_lines[name], axis.lines[0].get_ydata())
        assert axis.lines[0].get_zorder() < explorer.fft_lines[name].get_zorder()
        highest = max(np.max(line.get_ydata()) for line in axis.lines)
        assert axis.get_ylim()[1] > highest
        assert axis.get_ylim()[0] == 0
    for group in sweep.GROUPS:
        explorer.change_group(group)
        explorer.figure.canvas.draw()
    explorer.reset()
    finish()
    np.testing.assert_allclose(explorer.fft_lines["full"].get_ydata(), original)
    explorer.change_group("E2 Constants")
    explorer.figure.savefig(root / "meeting11/parameter_sweep_licinio.png", dpi=130)
    explorer.fft_figure.savefig(root / "meeting11/parameter_sweep_fft.png", dpi=130)
    print("Clock alignment, gray reference, live FFT, reset, all tabs, and y limits passed.", flush=True)
finally:
    explorer.close()
    assert not plt.fignum_exists(explorer.fft_figure.number)
    plt.close(explorer.figure)
