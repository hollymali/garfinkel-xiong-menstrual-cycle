"""
Chaotic E2 time series (n=14, in the chaotic window) -- with and without the
circadian envelope. Shown over 2.5 days: aperiodic, no two days alike, no noise.
"""

import numpy as np
import matplotlib.pyplot as plt

DAY = 1440.0
def circadian(t, t_peak=14*60):
    return 0.5*(1.0 + np.cos(2*np.pi*(np.mod(t, DAY) - t_peak)/DAY))

def run(n, K=14.8, Anight=1.0, t_end=8*DAY, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0,
        Pmin=40.0, Pmax=140.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    E_arr = np.empty(steps); T = np.arange(steps)*dt
    G, F, E = y0; phi = 0.0
    for i in range(steps-1):
        e = max(E, 0.0)
        period = Pmin + (Pmax - Pmin)/(1.0 + (e/K)**n)   # positive freq fb
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi); phi = phi_new
        a = 0.0
        if pulsed:
            A_circ = A*(Anight + (1.0-Anight)*circadian(T[i]))
            a = A_circ/(1.0 + (e/K)**n)                  # negative amp fb (+circadian)
        G = G + dt*(-kG*G) + (a if pulsed else 0.0)
        F = F + dt*(cF*G - kF*F)
        E = E + dt*(cE*F - kE*E)
        E_arr[i+1] = E
    E_arr[0] = y0[2]
    return T, E_arr


if __name__ == "__main__":
    n = 14
    T, E_flat = run(n, Anight=1.0)     # no circadian
    _, E_circ = run(n, Anight=0.35)    # with circadian

    start = 4*DAY
    m = (T >= start) & (T < start + 2.5*DAY)
    th = (T[m]-start)/60.0             # hours from start

    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    ax[0].plot(th, E_flat[m], lw=0.8, color="C0")
    ax[0].set_ylabel("E2 (pg/mL)")
    ax[0].set_title(f"High-period cascade (n={n}), no circadian: ragged-looking but NOT chaos (IC-insensitive)")
    ax[0].grid(alpha=0.3)
    ax[1].plot(th, E_circ[m], lw=0.8, color="C4")
    ax[1].set_ylabel("E2 (pg/mL)")
    ax[1].set_title("Same trace + circadian envelope (the full panel-B effect)")
    ax[1].set_xlabel("time (h)"); ax[1].grid(alpha=0.3)
    for a in ax:
        for d in range(0, 3):
            a.axvspan(d*24+12, d*24+24, color="0.85", alpha=0.5, zorder=0)  # crude night shading
    fig.tight_layout(); fig.savefig("fig_chaos_timeseries.png", dpi=130)
    print("wrote fig_chaos_timeseries.png")

    # quick sensitivity check: two runs with a tiny IC difference should diverge (chaos)
    T1, E1 = run(n, Anight=1.0, y0=(0.0,5.0,15.0))
    T2, E2b = run(n, Anight=1.0, y0=(0.0,5.0,15.0001))
    d = np.abs(E1 - E2b)
    print(f"tiny IC diff 1e-4 -> E2 divergence after 8 days: {d[-4000:].mean():.2f} pg/mL "
          f"(large => sensitive dependence = chaos)")
