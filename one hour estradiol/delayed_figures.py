"""
Delayed physiological regime (tau=120 min):
  Fig A: E2 time series from a high-orbit cell (n_amp=4, n_freq=3)
  Fig B: sweep n_freq (n_amp=3 fixed)   -- physiological + delay
  Fig C: sweep n_amp  (n_freq=3 fixed)   -- physiological + delay
"""
import numpy as np
import matplotlib.pyplot as plt
from MODEL import run
from orbit_period import e2_orbit_period, period_label

TAU = 120.0

# --- Fig A: high-orbit time series ---
na, nf = 3, 3
T, _, _, E = run(na, nf, TAU)
lab = period_label(e2_orbit_period(T, E, tail=8000.0))
m = T > T[-1] - 3*1440
fig, ax = plt.subplots(figsize=(12, 3.8))
ax.plot((T[m]-T[m][0])/60, E[m], lw=1.0, color="C4")
for x in range(0, 72, 12): ax.axvline(x, color="0.85", lw=0.6)
ax.set_xlabel("time (h)"); ax.set_ylabel("E2 (pg/mL)")
ax.set_title(f"High-orbit E2 time series (3 days):  n_amp={na}, n_freq={nf}, tau={TAU:.0f} min  ->  {lab}")
ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig("fig_delayed_timeseries.png", dpi=130)
print(f"Fig A: n_amp={na},n_freq={nf} -> {lab}")

# --- Figs B/C: sweep n_freq or n_amp, at a given delay tau ---
def sweep_fig(fixed_label, fixed_val, sweep_name, sweep_vals, fname, color, tau):
    delay_txt = f"tau={tau:.0f} min" if tau else "NO DELAY (tau=0)"
    fig, ax = plt.subplots(len(sweep_vals), 1, figsize=(11, 10), sharex=True)
    for k, v in enumerate(sweep_vals):
        na_, nf_ = (fixed_val, v) if sweep_name == "n_freq" else (v, fixed_val)
        T, _, _, E = run(na_, nf_, tau)
        lab = period_label(e2_orbit_period(T, E, tail=8000.0))
        m = T > T[-1] - 48*60
        ax[k].plot((T[m]-T[m][0])/60, E[m], lw=0.9, color=color)
        ax[k].set_ylabel(f"{sweep_name}={v}\nE2"); ax[k].grid(alpha=0.3)
        ax[k].set_title(f"{fixed_label}, {sweep_name}={v}, {delay_txt}:   orbit = {lab}",
                        loc="left", fontsize=9)
        print(f"  {sweep_name}={v}: {lab}")
    ax[-1].set_xlabel("time (h)  (2 days)")
    fig.suptitle(f"sweep {sweep_name} ({fixed_label} fixed) -- {delay_txt}")
    fig.tight_layout(); fig.savefig(fname, dpi=130)

for tau, tag in [(TAU, "delayed"), (0, "nodelay")]:
    print(f"sweep n_freq (n_amp=4), {tag}:")
    sweep_fig("n_amp=4", 4, "n_freq", [0, 1, 2, 3, 4, 5], f"fig_{tag}_freqsweep.png", "C0", tau)
    print(f"sweep n_amp (n_freq=4), {tag}:")
    sweep_fig("n_freq=4", 4, "n_amp", [0, 1, 2, 3, 4, 5], f"fig_{tag}_ampsweep.png", "C1", tau)
print("wrote fig_delayed_timeseries.png + delayed/nodelay freq & amp sweeps")
