"""
PI cascade + E2 -> GnRH  POSITIVE  frequency feedback.
=====================================================

Cascade (unchanged):
    dGnRH/dt = -kG*GnRH  (+ A at each pulse)
    dFSH/dt  =  cF*GnRH - kF*FSH
    dE2/dt   =  cE*FSH   - kE*E2

Pulse timing: positive feedback -- higher E2 -> SHORTER period (faster pulses):
    period(E2) = Pmin + (Pmax - Pmin) / (1 + (E2/K)^n)
    (E2=0 -> period=Pmax slow;  E2>>K -> period=Pmin fast;  n = steepness)
    dphi/dt = 1/period(E2);  pulse fires each time phi crosses an integer.

Pulse amplitude is FIXED (A) so we see the frequency feedback ALONE.
Sweep n to see what the positive feedback does.
"""

import numpy as np
import matplotlib.pyplot as plt

def run(n_steep, K, t_end=9000.0, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0,
        Pmin=40.0, Pmax=140.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    G = np.empty(steps); F = np.empty(steps); E = np.empty(steps); T = np.arange(steps)*dt
    G[0], F[0], E[0] = y0
    phi = 0.0
    pulse_t = []
    for i in range(steps-1):
        e = max(E[i], 0.0)
        period = Pmin + (Pmax - Pmin) / (1.0 + (e/K)**n_steep)   # positive: high E2 -> short period
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi)
        phi = phi_new
        if pulsed:
            pulse_t.append(T[i])
        G[i+1] = G[i] + dt*(-kG*G[i]) + (A if pulsed else 0.0)
        F[i+1] = F[i] + dt*(cF*G[i] - kF*F[i])
        E[i+1] = E[i] + dt*(cE*F[i] - kE*E[i])
    return T, G, F, E, np.array(pulse_t)


if __name__ == "__main__":
    # natural E2 with a fixed 90-min period, used as the feedback midpoint K
    Tn, Gn, Fn, En, _ = run(0.0, K=1e9)   # n=0 -> period=Pmin+(Pmax-Pmin)/2 const
    K = En[Tn > Tn[-1]-3000].mean()
    print(f"natural E2 (fixed period) = {K:.1f}  -> feedback midpoint K")

    from e2metrics import e2_orbit_period, period_label
    steeps = [1, 2, 4, 8]
    fig, ax = plt.subplots(len(steeps), 1, figsize=(11, 9), sharex=True)
    for k, n in enumerate(steeps):
        T, G, F, E, pt = run(n, K=K)
        lab = period_label(e2_orbit_period(T, E))
        m = T > T[-1] - 2400
        ax[k].plot(T[m]/60, E[m], lw=1.0, color="C2")
        ax[k].set_ylabel(f"n={n}\nE2"); ax[k].grid(alpha=0.3)
        ax[k].set_title(f"n_freq={n}:   orbit = {lab}", loc="left", fontsize=9)
        print(f"n_freq={n}: orbit {lab}")
    ax[-1].set_xlabel("time (h)")
    fig.suptitle("PI cascade + E2->GnRH POSITIVE frequency feedback (fixed amplitude): effect of steepness n")
    fig.tight_layout()
    fig.savefig("fig_pi_cascade_freqpos.png", dpi=130)
    print("wrote fig_pi_cascade_freqpos.png")
