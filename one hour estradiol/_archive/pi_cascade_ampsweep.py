"""
Complement of figure [4]: hold FREQUENCY feedback strong (n_freq = 8) and
sweep the AMPLITUDE feedback steepness n_amp.  (exploratory)
"""
import numpy as np
import matplotlib.pyplot as plt

def run(n_amp, n_freq, K, t_end=12000.0, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0,
        Pmin=40.0, Pmax=140.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    G = np.empty(steps); F = np.empty(steps); E = np.empty(steps); T = np.arange(steps)*dt
    G[0], F[0], E[0] = y0; phi = 0.0; a = 0.0
    pulse_t = []
    for i in range(steps-1):
        e = max(E[i], 0.0)
        period = Pmin + (Pmax - Pmin)/(1.0 + (e/K)**n_freq)   # pos frequency (fixed n_freq)
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi); phi = phi_new
        if pulsed:
            a = A/(1.0 + (e/K)**n_amp)                        # neg amplitude (swept n_amp)
            pulse_t.append(T[i])
        G[i+1] = G[i] + dt*(-kG*G[i]) + (a if pulsed else 0.0)
        F[i+1] = F[i] + dt*(cF*G[i] - kF*F[i])
        E[i+1] = E[i] + dt*(cE*F[i] - kE*E[i])
    return T, E, np.array(pulse_t)


if __name__ == "__main__":
    K = 14.8; n_freq = 8
    from e2metrics import e2_orbit_period, period_label
    n_amps = [0, 3, 8, 16]
    fig, ax = plt.subplots(len(n_amps), 1, figsize=(11, 9), sharex=True)
    for k, na in enumerate(n_amps):
        T, E, pt = run(na, n_freq, K)
        lab = period_label(e2_orbit_period(T, E))
        m = T > T[-1] - 2400
        ax[k].plot(T[m]/60, E[m], lw=0.9, color="C1")
        ax[k].set_ylabel(f"n_amp={na}\nE2"); ax[k].grid(alpha=0.3)
        ax[k].set_title(f"n_freq=8 (fixed), n_amp={na}:   orbit = {lab}", loc="left", fontsize=9)
        print(f"n_amp={na}: orbit {lab}")
    ax[-1].set_xlabel("time (h)")
    fig.suptitle("Hold frequency feedback strong (n_freq=8), sweep amplitude feedback n_amp")
    fig.tight_layout(); fig.savefig("fig_ampsweep.png", dpi=130)
    print("wrote fig_ampsweep.png")
