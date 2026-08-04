"""
Circadian drive + E2->GnRH feedback (the faithful version).
===========================================================

Same as pure_circadian.py, but GnRH pulse freq/amp are now modulated by BOTH:
  - the exogenous circadian clock C(t)   (as before), AND
  - E2 negative feedback: higher E2 -> longer period (slower) & smaller amplitude
    (the Rasgon-style E2->GnRH arm we had dropped for the pure-circadian demo).

Question: does the E2 feedback REINFORCE or OPPOSE the circadian envelope, and
does it add pulse-to-pulse irregularity? Run with feedback off vs on to compare.
"""

import numpy as np
import matplotlib.pyplot as plt

DAY = 1440.0
SIG = 0.04
T_PEAK = 14.0 * 60
def circadian(t):
    tod = np.mod(t, DAY)
    return 0.5 * (1.0 + np.cos(2*np.pi * (tod - T_PEAK) / DAY))

P_DAY, P_NIGHT = 48.0, 100.0
A_MIN, A_MAX   = 0.30, 1.00
def base_period(C): return P_NIGHT - (P_NIGHT - P_DAY) * C
def base_amp(C):    return A_MIN + (A_MAX - A_MIN) * C

def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIG) ** 2)

kE2  = 0.022
A_E2 = 6.5
BASAL = 0.18

E2_REF = 27.0          # reference E2 for the feedback (≈ no-fb daily mean)

def simulate(t_end=5*DAY, dt=0.5, g_freq=0.0, g_amp=0.0, e0=15.0):
    n = int(round(t_end/dt)); T = np.arange(n)*dt
    E2 = np.empty(n); PHI = np.empty(n)
    E2[0], PHI[0] = e0, 0.0
    for i in range(n-1):
        C = circadian(T[i]); e = E2[i]
        fb = (e - E2_REF) / E2_REF                      # relative E2 deviation
        period = base_period(C) * (1.0 + g_freq * fb)   # higher E2 -> longer period
        period = min(max(period, 20.0), 220.0)
        amp    = max(base_amp(C) * (1.0 - g_amp * fb), 0.05)  # higher E2 -> smaller
        g = amp * pulse_kernel(PHI[i])
        E2[i+1]  = e + dt * (A_E2 * g + BASAL - kE2 * e)
        PHI[i+1] = PHI[i] + dt * (1.0 / period)
    return T, E2


def day_slice(T, E2):
    day_start = 3*DAY + 8*60
    m = (T >= day_start) & (T < day_start + DAY)
    tod = (T[m] - day_start) / 60.0 + 8.0
    return tod, E2[m]


if __name__ == "__main__":
    T0, E0 = simulate(g_freq=0.0, g_amp=0.0)          # feedback OFF (= pure circadian)
    T1, E1 = simulate(g_freq=0.45, g_amp=0.45)        # feedback ON

    for name, (T, E) in [("no fb", (T0, E0)), ("E2->GnRH fb", (T1, E1))]:
        _, e = day_slice(T, E)
        print(f"{name:14s}: E2 {e.min():5.1f}..{e.max():5.1f}  mean {e.mean():5.1f}  "
              f"hump p2p {e.max()-e.min():5.1f}")

    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for k, (title, T, E) in enumerate(
            [("Circadian only (no E2->GnRH feedback)", T0, E0),
             ("Circadian + E2->GnRH negative feedback", T1, E1)]):
        tod, e = day_slice(T, E)
        ax[k].plot(tod, e, color=("C7" if k == 0 else "C3"), lw=1.1)
        ax[k].set_ylabel("Oestradiol\n(pg/mL)"); ax[k].set_title(title, loc="left", fontsize=10)
        ax[k].grid(alpha=0.3)
    ax[-1].set_xlabel("time of day (h)")
    ax[-1].set_xticks([8,14,20,26,32])
    ax[-1].set_xticklabels(["08:00","14:00","20:00","02:00","08:00"])
    fig.suptitle("Does E2->GnRH feedback reinforce or oppose the circadian envelope?")
    fig.tight_layout()
    fig.savefig("fig_circadian_plus_fb.png", dpi=130)
    print("wrote fig_circadian_plus_fb.png")
