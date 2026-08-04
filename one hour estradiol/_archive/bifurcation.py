"""
Bifurcation of the combined cascade as feedback steepness n increases
(n_amp = n_freq = n).  E2-based: plot the E2 PEAK HEIGHTS the system settles onto.
  1 value  -> period-1 ; 2 -> period-2 ; 4,8 -> period-doubling ; band -> high-period.
Bottom panel: E2_period_var (MAD/median of E2 peak spacings).
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from e2metrics import rcv

def run(n, K=14.8, t_end=8000.0, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0,
        Pmin=40.0, Pmax=140.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    E_arr = np.empty(steps); T = np.arange(steps)*dt
    G, F, E = y0; phi = 0.0; E_arr[0] = E
    for i in range(steps-1):
        e = max(E, 0.0)
        period = Pmin + (Pmax - Pmin)/(1.0 + (e/K)**n)
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi); phi = phi_new
        a = A/(1.0 + (e/K)**n) if pulsed else 0.0
        G = G + dt*(-kG*G) + a
        F = F + dt*(cF*G - kF*F)
        E = E + dt*(cE*F - kE*E)
        E_arr[i+1] = E
    return T, E_arr

def e2_peaks(T, E, tail=3500.0):
    m = T > T[-1] - tail
    Tt, Ee = T[m], E[m]
    rng = Ee.max() - Ee.min()
    if rng < 1: return np.array([]), np.array([])
    pk, _ = find_peaks(Ee, prominence=0.15*rng)
    return Ee[pk], np.diff(Tt[pk])


if __name__ == "__main__":
    ns = np.linspace(2.0, 26.0, 120)
    fig, ax = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
    for n in ns:
        T, E = run(n)
        h, iv = e2_peaks(T, E)
        ax[0].plot([n]*len(h),  h,  ".", ms=1.6, color="C4", alpha=0.5)  # E2 peak HEIGHTS
        ax[1].plot([n]*len(iv), iv, ".", ms=1.6, color="C3", alpha=0.5)  # E2 peak-to-peak PERIODS
    ax[0].set_ylabel("settled E2 peak heights\n(attractor)")
    ax[0].set_title("Bifurcation of the combined cascade (E2), n_amp = n_freq = n\n"
                    "top: peak HEIGHTS   bottom: peak-to-peak PERIODS   (1 dot=period-1, 2=period-2, band=high-period)")
    ax[0].grid(alpha=0.3)
    ax[1].set_ylabel("settled E2 peak-to-peak\nperiods (min)")
    ax[1].set_xlabel("feedback steepness  n (= n_amp = n_freq)")
    ax[1].grid(alpha=0.3)
    # denser ticks on both panels
    ax[1].xaxis.set_major_locator(MultipleLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0].yaxis.set_major_locator(MultipleLocator(2)); ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].yaxis.set_major_locator(MultipleLocator(20)); ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))
    for a in ax:
        a.grid(which="minor", alpha=0.12)
        a.tick_params(which="both", direction="out")
    fig.tight_layout(); fig.savefig("fig_bifurcation.png", dpi=130)
    print("wrote fig_bifurcation.png")
