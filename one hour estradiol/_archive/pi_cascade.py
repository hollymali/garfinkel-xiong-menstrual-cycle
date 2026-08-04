"""
PI's minimal cascade -> deterministic irregular periods & alternating amplitude.
================================================================================

Three first-order stages (this is the key -- 3 lags give >180 deg phase, enough
for negative feedback to destabilize; 2 lags never could):

    dGnRH/dt = pulse(phi) - kG*GnRH        # GnRH pulse + first-order decay
    dFSH/dt  = cF*GnRH    - kF*FSH         # GnRH-dependent FSH inflow + decay
    dE2/dt   = cE*FSH     - kE*E2          # FSH-dependent E2 inflow + decay

Closed by E2 -> GnRH pulse FREQUENCY feedback (period rises with E2, steepness n):
    period(E2) = Pmin + (Pmax-Pmin) * E2^n / (K^n + E2^n)
    dphi/dt    = 1/period(E2);   a pulse fires each time phi crosses an integer.

Crank the steepness n and the frequency-map period-doubles:
  regular pulses -> alternating high/low -> irregular (period-4 / chaotic).
NO noise anywhere; all irregularity is deterministic.
"""

import numpy as np
import matplotlib.pyplot as plt

def run(n_steep, K, t_end=9000.0, dt=0.25, feedback=True,
        kG=0.05, kF=0.05, kE=0.05, cF=0.245, cE=0.245, A=5.0,
        Pmin=20.0, Pmax=180.0, y0=(0.0, 5.0, 40.0)):
    steps = int(t_end/dt)
    G = np.empty(steps); F = np.empty(steps); E = np.empty(steps); T = np.arange(steps)*dt
    G[0], F[0], E[0] = y0
    phi = 0.0
    pulse_t, pulse_E2 = [], []
    for i in range(steps-1):
        e = max(E[i], 0.0)
        if feedback:
            h = e**n_steep / (K**n_steep + e**n_steep)
            period = Pmin + (Pmax - Pmin) * h
        else:
            period = 60.0
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi)
        phi = phi_new
        if pulsed:
            pulse_t.append(T[i]); pulse_E2.append(e)
        G[i+1] = G[i] + dt*(-kG*G[i]) + (A if pulsed else 0.0)
        F[i+1] = F[i] + dt*(cF*G[i] - kF*F[i])
        E[i+1] = E[i] + dt*(cE*F[i] - kE*E[i])
    return T, G, F, E, np.array(pulse_t), np.array(pulse_E2)


if __name__ == "__main__":
    # set the feedback midpoint K to the natural E2 level (feedback off)
    T, G, F, E, *_ = run(1.0, K=100.0, feedback=False)
    K = E[T > T[-1]-3000].mean()
    print(f"natural E2 level (no feedback) = {K:.1f}  -> use as feedback midpoint K")

    steeps = [10, 18, 26, 36]
    fig, ax = plt.subplots(len(steeps), 1, figsize=(11, 9), sharex=True)
    for k, n in enumerate(steeps):
        T, G, F, E, pt, pe = run(n, K=K)
        # pulse-to-pulse intervals over the settled tail
        tail = pt[pt > T[-1]-3000]
        iv = np.diff(tail)
        m = T > T[-1] - 2400   # last 40 h
        ax[k].plot(T[m]/60, E[m], lw=0.9, color="C0")
        ax[k].set_ylabel(f"n={n}\nE2")
        ax[k].grid(alpha=0.3)
        if len(iv) > 2:
            ax[k].set_title(f"steepness n={n}:  interval {iv.min():.0f}-{iv.max():.0f} min, "
                            f"CV={iv.std()/iv.mean():.2f}   E2 {E[m].min():.0f}-{E[m].max():.0f}",
                            loc="left", fontsize=9)
        seq = " ".join(f"{v:.0f}" for v in iv[-8:])
        print(f"n={n:2d}: CV {iv.std()/iv.mean():.3f}  last intervals(min): {seq}")
    ax[-1].set_xlabel("time (h)")
    fig.suptitle("PI cascade (GnRH->FSH->E2, 3 first-order lags + E2->frequency feedback): "
                 "steepness drives period-doubling")
    fig.tight_layout()
    fig.savefig("fig_pi_cascade.png", dpi=130)
    print("wrote fig_pi_cascade.png")
