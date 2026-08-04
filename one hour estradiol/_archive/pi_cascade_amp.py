"""
PI cascade, AMPLITUDE-feedback version -> clean alternating high/low pulses.
===========================================================================

Same 3 first-order stages (GnRH -> FSH -> E2), but the E2 feedback now sets the
next pulse's AMPLITUDE (not its timing):

    a_n = A / (1 + (E2_at_pulse / K)^n)      # higher E2 -> smaller next pulse
    pulses fire at a FIXED interval T0        # regular spacing -> resolved pulses

Steep enough n -> the pulse-amplitude map period-doubles -> tall, short, tall,
short... distinct, fully-resolved pulses (no shoulder bumps). Deterministic.
"""

import numpy as np
import matplotlib.pyplot as plt

def run(n_steep, K, T0=90.0, t_end=9000.0, dt=0.25, feedback=True,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0, y0=(0.0, 5.0, 15.0)):
    steps = int(t_end/dt)
    G = np.empty(steps); F = np.empty(steps); E = np.empty(steps); T = np.arange(steps)*dt
    G[0], F[0], E[0] = y0
    phi = 0.0; a = 0.0
    pulse_t, pulse_amp = [], []
    for i in range(steps-1):
        t = T[i]; e = max(E[i], 0.0)
        phi_new = phi + dt/T0                          # dphi/dt = 1/T0 (constant period)
        pulsed = np.floor(phi_new) > np.floor(phi)     # pulse when phi crosses an integer
        phi = phi_new
        if pulsed:
            a = A / (1.0 + (e/K)**n_steep) if feedback else A/2.0
            pulse_t.append(t); pulse_amp.append(a)
        G[i+1] = G[i] + dt*(-kG*G[i]) + (a if pulsed else 0.0)
        F[i+1] = F[i] + dt*(cF*G[i] - kF*F[i])
        E[i+1] = E[i] + dt*(cE*F[i] - kE*E[i])
    return T, G, F, E, np.array(pulse_t), np.array(pulse_amp)


if __name__ == "__main__":
    T, G, F, E, *_ = run(1.0, K=40.0, feedback=False)
    K = E[T > T[-1]-3000].mean()
    print(f"natural E2 (no feedback) = {K:.1f}  -> feedback midpoint K")

    from e2metrics import e2_orbit_period, period_label
    steeps = [3, 5, 8, 12]
    fig, ax = plt.subplots(len(steeps), 1, figsize=(11, 9), sharex=True)
    for k, n in enumerate(steeps):
        T, G, F, E, pt, pa = run(n, K=K)
        lab = period_label(e2_orbit_period(T, E))
        m = T > T[-1] - 2000
        ax[k].plot(T[m]/60, E[m], lw=1.0, color="C0")
        ax[k].set_ylabel(f"n={n}\nE2"); ax[k].grid(alpha=0.3)
        ax[k].set_title(f"n_amp={n}:   orbit = {lab}", loc="left", fontsize=9)
        print(f"n_amp={n:2d}: orbit {lab}")
    ax[-1].set_xlabel("time (h)")
    fig.suptitle("PI cascade with E2->AMPLITUDE feedback (fixed 75-min pulses): "
                 "clean alternating high/low E2 pulses")
    fig.tight_layout()
    fig.savefig("fig_pi_cascade_amp.png", dpi=130)
    print("wrote fig_pi_cascade_amp.png")
