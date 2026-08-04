"""
FREEZE TEST -- is E2 a leaky integrator, or can it self-oscillate?
==================================================================

Hypothesis (the slow-fast picture): GnRH is the only oscillator; E2 is a leaky
integrator that passively filters the GnRH->FSH drive through its slow decay.
If true, then holding the GnRH drive CONSTANT must make E2 (and FSH) SETTLE --
no intrinsic rhythm of their own.

Test:
  (1) pulsatile : normal model, GnRH pulses with E2-modulated freq/amp.
  (2) frozen    : replace the GnRH->FSH drive with its time-average constant,
                  leaving E2's own decay + e2_to_fsh feedback intact.
Also run (2) from several initial conditions to confirm one fixed point
(monostable), i.e. a settling filter, not a limit cycle.

No delay anywhere -- this is the pure leaky-integrator model.
"""

import numpy as np
import matplotlib.pyplot as plt
from efp_two_loop import (pulse_kernel, fsh_to_e2, gnrh_to_fsh, e2_to_fsh,
                          gnrh_period, gnrh_amp, kE2, kFSH, cFSH)

CE2 = 0.012
K_E2 = 400.0

def run(mode, drive_const=None, t_end=8000.0, dt=0.5, y0=(300.0, 5.0, 0.0)):
    n = int(round(t_end / dt)); T = np.arange(n) * dt
    E2 = np.empty(n); FSH = np.empty(n); PHI = np.empty(n)
    E2[0], FSH[0], PHI[0] = y0

    def deriv(e, f, p):
        drive = drive_const if mode == "frozen" else gnrh_to_fsh(gnrh_amp(e) * pulse_kernel(p))
        dE2  = CE2  * fsh_to_e2(f)                     - kE2  * e
        dFSH = cFSH * (e2_to_fsh(e, K_E2) + drive)     - kFSH * f
        dphi = 1.0 / gnrh_period(e)
        return dE2, dFSH, dphi

    for i in range(n - 1):
        e, f, p = E2[i], FSH[i], PHI[i]
        k1 = deriv(e, f, p)
        k2 = deriv(e+dt/2*k1[0], f+dt/2*k1[1], p+dt/2*k1[2])
        k3 = deriv(e+dt/2*k2[0], f+dt/2*k2[1], p+dt/2*k2[2])
        k4 = deriv(e+dt*k3[0], f+dt*k3[1], p+dt*k3[2])
        E2[i+1]  = e + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        FSH[i+1] = f + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        PHI[i+1] = p + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    return T, E2, FSH


if __name__ == "__main__":
    # (1) pulsatile: run the normal leaky-integrator model
    T, E2, FSH = run("pulsatile")
    tailmask = T > T[-1] - 3000
    # Freeze the GnRH drive at the value implied by the settled means: at steady
    # state dFSH=0 => drive = (kFSH/cFSH)*FSH_mean - e2_to_fsh(E2_mean).
    E2m = E2[tailmask].mean(); FSHm = FSH[tailmask].mean()
    Dbar = (kFSH / cFSH) * FSHm - e2_to_fsh(E2m, K_E2)
    print(f"pulsatile settled means: E2={E2m:.2f}  FSH={FSHm:.3f}  "
          f"implied mean GnRH drive Dbar={Dbar:.3f}")

    # (2) frozen: hold GnRH drive at Dbar, from several ICs
    print("\n--- FROZEN GnRH (constant drive): does E2 settle or oscillate? ---")
    finals = []
    for y0 in [(150., 3., 0.), (300., 5., 0.), (500., 8., 0.), (250., 6., 0.)]:
        Tz, E2z, FSHz = run("frozen", drive_const=Dbar, y0=y0)
        tail = E2z[Tz > Tz[-1]-1500]
        finals.append(E2z[-1])
        print(f"  IC E2={y0[0]:4.0f} -> final E2={E2z[-1]:7.3f}  "
              f"tail p2p={tail.max()-tail.min():.4f}")
    print(f"  spread of finals across ICs = {max(finals)-min(finals):.4f}  "
          f"(->0 means single fixed point, monostable filter)")

    # figure: pulsatile vs frozen
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    m = T > T[-1] - 40*60
    ax[0].plot(T[m]/60, E2[m], lw=0.8, color="C0")
    ax[0].set_ylabel("E2\n(pulsatile)"); ax[0].grid(alpha=0.3)
    ax[0].set_title("Freeze test: E2 is a leaky integrator (settles when GnRH held constant)")
    Tz, E2z, FSHz = run("frozen", drive_const=Dbar, y0=(150., 3., 0.))
    ax[1].plot(Tz/60, E2z, lw=1.2, color="C3")
    ax[1].set_ylabel("E2\n(GnRH frozen)"); ax[1].set_xlabel("time (h)"); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("fig_freeze_test.png", dpi=130)
    print("\nwrote fig_freeze_test.png")
