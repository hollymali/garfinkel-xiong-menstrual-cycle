"""
Combined cascade (both E2 feedbacks) + exogenous CIRCADIAN modulation.
=====================================================================

Cascade + E2 feedbacks (as in pi_cascade_combined):
    a_n        = A_circ(t) / (1 + (E2/K_amp)^n_amp)          # neg amplitude fb
    period(E2) = Pmin + (Pmax-Pmin)/(1 + (E2/K_freq)^n_freq) # pos frequency fb

NEW: the SCN clock C(t) (0..1, peaks ~14:00) scales the pulse amplitude --
daytime pulses stronger, night pulses weaker:
    A_circ(t) = A * (Anight + (1-Anight)*C(t))
C(t) is an imposed exogenous 24 h signal (the master clock output), NOT a
state variable -- it runs independently of E2.
"""

import numpy as np
import matplotlib.pyplot as plt

DAY = 1440.0
def circadian(t, t_peak=14*60):
    return 0.5*(1.0 + np.cos(2*np.pi*(np.mod(t, DAY) - t_peak)/DAY))

def run(n_amp, n_freq, K_amp, K_freq, Anight=0.35, t_end=8*DAY, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0,
        Pmin=40.0, Pmax=140.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    G = np.empty(steps); F = np.empty(steps); E = np.empty(steps); T = np.arange(steps)*dt
    G[0], F[0], E[0] = y0
    phi = 0.0
    for i in range(steps-1):
        e = max(E[i], 0.0)
        period = Pmin + (Pmax - Pmin)/(1.0 + (e/K_freq)**n_freq)
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi)
        phi = phi_new
        if pulsed:
            A_circ = A*(Anight + (1.0-Anight)*circadian(T[i]))
            a = A_circ/(1.0 + (e/K_amp)**n_amp)
        G[i+1] = G[i] + dt*(-kG*G[i]) + (a if pulsed else 0.0)
        F[i+1] = F[i] + dt*(cF*G[i] - kF*F[i])
        E[i+1] = E[i] + dt*(cE*F[i] - kE*E[i])
    return T, G, F, E


if __name__ == "__main__":
    K = 14.8   # shared feedback midpoint (natural E2), from earlier runs
    T, G, F, E = run(n_amp=10, n_freq=4, K_amp=K, K_freq=K)

    # show one day, 08:00 -> 08:00, past transient
    day_start = 6*DAY + 8*60
    m = (T >= day_start) & (T < day_start + DAY)
    tod = (T[m]-day_start)/60.0 + 8.0
    print(f"E2 over the day: {E[m].min():.1f}..{E[m].max():.1f}  mean {E[m].mean():.1f}")

    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True,
                           gridspec_kw={"height_ratios":[1, 2.6]})
    ax[0].plot(tod, circadian(T[m]), color="C1", lw=1.6)
    ax[0].fill_between(tod, 0, circadian(T[m]), color="C1", alpha=0.15)
    ax[0].set_ylabel("circadian\nC(t)"); ax[0].grid(alpha=0.3)
    ax[0].set_title("Combined cascade (neg amp + pos freq E2 feedback) + circadian amplitude drive")
    ax[1].plot(tod, E[m], color="C4", lw=1.0)
    ax[1].set_ylabel("Oestradiol (pg/mL)"); ax[1].set_xlabel("time of day (h)"); ax[1].grid(alpha=0.3)
    ax[1].set_xticks([8,14,20,26,32]); ax[1].set_xticklabels(["08:00","14:00","20:00","02:00","08:00"])
    fig.tight_layout(); fig.savefig("fig_pi_cascade_circadian.png", dpi=130)
    print("wrote fig_pi_cascade_circadian.png")
