"""
Add a TIME DELAY to the E2 feedback (Mackey-Glass route) and test for GENUINE chaos.

The pulse generator now reads E2(t - tau) instead of instantaneous E2:
    period(t) = Pmin + (Pmax-Pmin)/(1 + (E2(t-tau)/K)^n)     # positive freq
    a_n       = A / (1 + (E2(t-tau)/K)^n)                    # negative amp

Delayed nonlinear feedback is infinite-dimensional -> can be chaotic.
For each tau: run two trajectories 1e-6 apart and measure divergence.
  |dE2| grows to O(amplitude) -> CHAOS ; stays ~0 -> still periodic.
"""

import numpy as np

def run(n, tau, e0, K=14.8, t_end=12*1440.0, dt=0.25,
        kG=0.25, kF=0.10, kE=0.018, cF=0.5, cE=0.3, A=8.0, Pmin=40.0, Pmax=140.0):
    steps = int(t_end/dt); ndelay = int(round(tau/dt))
    E_arr = np.empty(steps); E_arr[0] = e0
    G, F, E = 0.0, 5.0, e0; phi = 0.0
    pulse_t = []
    for i in range(steps-1):
        j = i - ndelay
        e_del = max(E_arr[j] if j >= 0 else e0, 0.0)
        period = Pmin + (Pmax - Pmin)/(1.0 + (e_del/K)**n)
        phi_new = phi + dt/period
        pulsed = np.floor(phi_new) > np.floor(phi); phi = phi_new
        a = A/(1.0 + (e_del/K)**n) if pulsed else 0.0
        if pulsed: pulse_t.append(i*dt)
        G = G + dt*(-kG*G) + a
        F = F + dt*(cF*G - kF*F)
        E = E + dt*(cE*F - kE*E)
        E_arr[i+1] = E
    return E_arr, np.array(pulse_t)

if __name__ == "__main__":
    print("  n  tau(min)   final|dE2|   peak   interval CV   verdict")
    for n in [10, 16, 24]:
        for tau in [60, 120, 240, 480]:
            E1, pt = run(n, tau, 15.0)
            E2, _  = run(n, tau, 15.000001)
            sep = np.abs(E1 - E2)
            tail = sep[-8000:].mean(); peak = sep.max()
            iv = np.diff(pt[pt > pt[-1]-4000]); cv = iv.std()/iv.mean() if len(iv) > 1 else 0
            verdict = "CHAOS" if tail > 0.5 else "periodic/quasi"
            print(f"  {n:2d}  {tau:4d}     {tail:9.4f}   {peak:6.2f}   {cv:9.3f}    {verdict}")
