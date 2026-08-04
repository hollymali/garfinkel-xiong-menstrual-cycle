"""
Is the 'chaotic band' real chaos? Test sensitive dependence on initial
conditions across n. For each n, run two trajectories whose E2 start 1e-6 apart
and track their separation. Chaos -> separation grows to O(attractor). Periodic/
quasiperiodic -> separation stays tiny or ->0.
"""

import numpy as np

def run(n, e0, K=14.8, t_end=10*1440.0, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0, Pmin=40.0, Pmax=140.0):
    steps = int(t_end/dt)
    E_arr = np.empty(steps)
    G, F, E = 0.0, 5.0, e0; phi = 0.0
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
    return E_arr

if __name__ == "__main__":
    print("  n   final |dE2|   verdict")
    for n in [8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 24]:
        E1 = run(n, 15.0)
        E2 = run(n, 15.000001)
        sep = np.abs(E1 - E2)
        tail = sep[-8000:].mean()      # last ~8 h averaged
        peak = sep.max()
        verdict = "CHAOS (sensitive)" if tail > 0.5 else ("periodic/quasi" if peak < 1.0 else "transient-only growth")
        print(f"  {n:2d}   {tail:9.4f}   peak {peak:7.3f}   {verdict}")
