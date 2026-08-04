"""
Pure-circadian estradiol model (NO noise, NO E2->GnRH feedback).
================================================================

Test: can an exogenous 24 h circadian clock, modulating the ~1 h GnRH pulse
generator, make estradiol look like panel B of the body-temp/estradiol figure
(ragged ~1 h pulses under a slow midday hump) -- with no stochastic input?

Structure (clean Model 1):
  circadian clock C(t)   -- imposed ~24 h, peaks mid-afternoon
        |  modulates
  GnRH pulse generator   -- faster & taller by day, slower & smaller at night
        |  drives
  E2 leaky filter        -- dE2/dt = A*pulse - kE2*E2   (E2 = passive readout)

GnRH stays the only oscillator; E2 never self-oscillates; the clock is external.
"""

import numpy as np
import matplotlib.pyplot as plt

DAY = 1440.0            # min in 24 h
SIG = 0.04             # pulse width (fraction of interval)

# --- circadian clock: 0..1, peaks at t_peak (mid-afternoon ~14:00) ---------
T_PEAK = 14.0 * 60      # 14:00 in minutes
def circadian(t):
    tod = np.mod(t, DAY)
    return 0.5 * (1.0 + np.cos(2*np.pi * (tod - T_PEAK) / DAY))   # 1 at 14:00, 0 at 02:00

# --- circadian modulation of the GnRH pulse generator ----------------------
P_DAY, P_NIGHT = 48.0, 100.0     # pulse interval (min): fast by day, slow at night
A_MIN, A_MAX   = 0.30, 1.00      # pulse amplitude: small at night, tall by day
def gnrh_period(C):  return P_NIGHT - (P_NIGHT - P_DAY) * C
def gnrh_amp(C):     return A_MIN + (A_MAX - A_MIN) * C

def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIG) ** 2)

# --- E2 leaky filter -------------------------------------------------------
kE2  = 0.022           # 1/min -> tau ~ 45 min (E2 clearance)
A_E2 = 6.5             # pulse->E2 gain (tuned so peaks ~45-50 pg/mL)
BASAL = 0.18           # small basal E2 production -> night floor ~10 pg/mL


def simulate(t_end=5*DAY, dt=0.5, e0=15.0):
    n = int(round(t_end/dt)); T = np.arange(n)*dt
    E2 = np.empty(n); PHI = np.empty(n); G = np.empty(n)
    E2[0], PHI[0] = e0, 0.0
    for i in range(n-1):
        C = circadian(T[i])
        g = gnrh_amp(C) * pulse_kernel(PHI[i]); G[i] = g
        E2[i+1]  = E2[i]  + dt * (A_E2 * g + BASAL - kE2 * E2[i])
        PHI[i+1] = PHI[i] + dt * (1.0 / gnrh_period(C))
    G[-1] = gnrh_amp(circadian(T[-1])) * pulse_kernel(PHI[-1])
    return T, E2, G


if __name__ == "__main__":
    T, E2, G = simulate()
    # show a full day, plotted 08:00 -> 08:00 like the figure
    day_start = 3*DAY + 8*60
    m = (T >= day_start) & (T < day_start + DAY)
    tod = (T[m] - day_start) / 60.0 + 8.0   # hours, 8..32 -> wrap label

    tail = E2[m]
    print(f"E2 over the day: min {tail.min():.1f}  max {tail.max():.1f}  mean {tail.mean():.1f} pg/mL")
    # count pulses (local maxima) as a rough regularity check
    peaks = np.where((E2[m][1:-1] > E2[m][:-2]) & (E2[m][1:-1] > E2[m][2:]))[0]
    print(f"E2 peaks in the day: {len(peaks)}")

    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True,
                           gridspec_kw={"height_ratios": [1, 2.4]})
    ax[0].plot(tod, circadian(T[m]), color="C1", lw=1.6)
    ax[0].fill_between(tod, 0, circadian(T[m]), color="C1", alpha=0.15)
    ax[0].set_ylabel("circadian\nclock C(t)"); ax[0].grid(alpha=0.3)
    ax[0].set_title("Pure-circadian E2 (no noise): clock-modulated GnRH pulses through a leaky filter")

    ax[1].plot(tod, E2[m], color="C3", lw=1.1)
    ax[1].set_ylabel("Oestradiol (pg/mL)"); ax[1].set_xlabel("time of day (h)")
    ax[1].grid(alpha=0.3)
    xt = [8, 14, 20, 2 + 24, 8 + 24]
    ax[1].set_xticks(xt); ax[1].set_xticklabels(["08:00","14:00","20:00","02:00","08:00"])
    fig.tight_layout()
    fig.savefig("fig_pure_circadian.png", dpi=130)
    print("wrote fig_pure_circadian.png")
