"""
Does a time delay stop E2 from settling?  Test delays on BOTH feedback arms:
  tau_fsh  : E2 -> FSH   uses  e2_to_fsh( E2(t - tau_fsh) )     [strong neg loop]
  tau_gnrh : E2 -> GnRH  uses  gnrh_period/amp( E2(t - tau_gnrh) )  [weaker loop]

Two-loop model (E2->GnRH freq/amp restored), cE2=0.012 so E2 rests near the
e2_to_fsh midpoint (loop gain present). If E2 stops settling -> tail peak-to-peak
of the slow envelope grows from ~0.
"""

import numpy as np
from efp_two_loop import (pulse_kernel, fsh_to_e2, gnrh_to_fsh, e2_to_fsh,
                          gnrh_period, gnrh_amp, kE2, kFSH, cFSH)

CE2, K_E2 = 0.012, 400.0

def simulate(tau_fsh=0.0, tau_gnrh=0.0, t_end=16000.0, dt=0.5, y0=(300.0, 5.0, 0.0)):
    n = int(round(t_end/dt)); T = np.arange(n)*dt
    E2 = np.empty(n); FSH = np.empty(n); PHI = np.empty(n)
    E2[0], FSH[0], PHI[0] = y0
    h0 = y0[0]
    nf = int(round(tau_fsh/dt)); ng = int(round(tau_gnrh/dt))

    def deriv(e, f, p, e_fsh, e_gnrh):
        G = gnrh_amp(e_gnrh) * pulse_kernel(p)
        dE2  = CE2  * fsh_to_e2(f)                            - kE2  * e
        dFSH = cFSH * (e2_to_fsh(e_fsh, K_E2) + gnrh_to_fsh(G)) - kFSH * f
        dphi = 1.0 / gnrh_period(e_gnrh)
        return dE2, dFSH, dphi

    for i in range(n-1):
        e, f, p = E2[i], FSH[i], PHI[i]
        ef = E2[i-nf] if i-nf >= 0 else h0
        eg = E2[i-ng] if i-ng >= 0 else h0
        k1 = deriv(e,            f,            p,            ef, eg)
        k2 = deriv(e+dt/2*k1[0], f+dt/2*k1[1], p+dt/2*k1[2], ef, eg)
        k3 = deriv(e+dt/2*k2[0], f+dt/2*k2[1], p+dt/2*k2[2], ef, eg)
        k4 = deriv(e+dt*k3[0],   f+dt*k3[1],   p+dt*k3[2],   ef, eg)
        E2[i+1]  = e + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        FSH[i+1] = f + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        PHI[i+1] = p + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    return T, E2, FSH

def p2p(T, E2, tail=4000.0):
    e = E2[T > T[-1]-tail]; return e.max()-e.min()

if __name__ == "__main__":
    print("baseline (no delay):")
    T, E2, _ = simulate(); print(f"  tail p2p = {p2p(T,E2):.3f}  (settles)\n")

    print("E2->FSH delay only (tau_gnrh=0):")
    for tf in [120, 240, 360, 480, 600]:
        T, E2, _ = simulate(tau_fsh=tf)
        print(f"  tau_fsh={tf:4d} min : tail p2p = {p2p(T,E2):7.2f}")

    print("\nE2->GnRH delay only (tau_fsh=0):")
    for tg in [120, 240, 360, 480, 600]:
        T, E2, _ = simulate(tau_gnrh=tg)
        print(f"  tau_gnrh={tg:4d} min : tail p2p = {p2p(T,E2):7.2f}")

    print("\nboth delays together:")
    for td in [120, 240, 360]:
        T, E2, _ = simulate(tau_fsh=td, tau_gnrh=td)
        print(f"  tau_fsh=tau_gnrh={td:4d} min : tail p2p = {p2p(T,E2):7.2f}")
