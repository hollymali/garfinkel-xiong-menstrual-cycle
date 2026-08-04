"""
Two-loop reduced EFP model  --  "one hour estradiol"
=====================================================

Adds the SECOND feedback loop back in: E2 modulates the GnRH pulse generator
(frequency AND amplitude), restored from the Rasgon EFP functions. So now:

  Loop 1 (delayed direct):  E2 --| FSH  (e2_to_fsh, delayed by tau)  -->  E2
  Loop 2 (via GnRH):        E2 --> GnRH freq/amp --> FSH (gnrh_to_fsh) --> E2

Two coupled feedback loops at different timescales is the substrate for
mode-locking / quasiperiodic / skipped-beat (non-uniform) E2 dynamics.

Key fix: the restored GnRH-modulation functions only respond for E2 in
[250, 400]; the reduced model otherwise rests E2 near ~570 (saturated, dead
modulation). We lower cE2 so E2 rests INSIDE the responsive band -- which also
lands E2 near the e2_to_fsh midpoint (400), restoring loop gain at the paper's
own midpoint. GnRH pulse phase `phi` is now a state variable (dphi/dt varies
in time because the period depends on E2).
"""

import numpy as np
import matplotlib.pyplot as plt

SIGMA_PHASE = 0.04
def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIGMA_PHASE) ** 2)

def hill_up(x, k, n, A, C):
    z = (k * x) ** n
    return A * z / (1.0 + z) + C

def fsh_to_e2(FSH):  return hill_up(FSH, 0.15, 7, 150.0, 250.0)   # stimulate
def gnrh_to_fsh(G):  return hill_up(G,   0.2, 15,  5.0,   5.0)    # drive

N_E2 = 15
def e2_to_fsh(E2, K_E2=400.0):
    return 9.0 / (1.0 + (E2 / K_E2) ** N_E2) + 1.0               # suppress

# ---- RESTORED E2 -> GnRH modulation (Rasgon EFP, native band E2 in [250,400]) ----
def gnrh_period(E2):     # minutes between pulses; higher E2 -> slower pulses
    return np.clip(0.267 * E2 - 6.75, 60.0, 100.0)
def gnrh_amp(E2):        # pulse height; higher E2 -> smaller pulses
    return np.clip(-0.013 * E2 + 8.2, 3.0, 5.0)

# ---- Parameters ----------------------------------------------------------
kE2, kFSH = 0.011, 0.009
cFSH = 5e-3
# cE2 LOWERED (was 0.025) so E2 rests inside the [250,400] modulation band.


def deriv(E2, FSH, phi, E2_delayed, K_E2, cE2):
    G = gnrh_amp(E2) * pulse_kernel(phi)
    dE2  = cE2  * fsh_to_e2(FSH)                                  - kE2  * E2
    dFSH = cFSH * (e2_to_fsh(E2_delayed, K_E2) + gnrh_to_fsh(G))  - kFSH * FSH
    dphi = 1.0 / gnrh_period(E2)
    return dE2, dFSH, dphi


def simulate(t_end=16000.0, dt=0.5, tau=0.0, cE2=0.016, K_E2=400.0,
             y0=(300.0, 5.0, 0.0)):
    n = int(round(t_end / dt))
    ndelay = int(round(tau / dt))
    T = np.arange(n) * dt
    E2 = np.empty(n); FSH = np.empty(n); PHI = np.empty(n)
    E2[0], FSH[0], PHI[0] = y0
    hist0 = y0[0]

    for i in range(n - 1):
        e, f, p = E2[i], FSH[i], PHI[i]
        j = i - ndelay
        ed = E2[j] if j >= 0 else hist0
        k1 = deriv(e,              f,              p,              ed, K_E2, cE2)
        k2 = deriv(e+dt/2*k1[0],   f+dt/2*k1[1],   p+dt/2*k1[2],   ed, K_E2, cE2)
        k3 = deriv(e+dt/2*k2[0],   f+dt/2*k2[1],   p+dt/2*k2[2],   ed, K_E2, cE2)
        k4 = deriv(e+dt*k3[0],     f+dt*k3[1],     p+dt*k3[2],     ed, K_E2, cE2)
        E2[i+1]  = e + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        FSH[i+1] = f + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        PHI[i+1] = p + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    return T, E2, FSH, PHI


if __name__ == "__main__":
    cE2 = 0.012          # E2 rests ~372: in the modulation band AND near the loop-gain midpoint
    taus = [0.0, 240.0, 360.0, 480.0, 600.0, 720.0]

    fig, ax = plt.subplots(len(taus), 1, figsize=(11, 12), sharex=True)
    for k, tau in enumerate(taus):
        T, E2, FSH, PHI = simulate(tau=tau, cE2=cE2, t_end=16000.0)
        tail = E2[T > T[-1]-4000]
        per = gnrh_period(E2)
        print(f"tau={tau:5.0f} : E2 tail {tail.min():6.1f}..{tail.max():6.1f} "
              f"(p2p {tail.max()-tail.min():6.1f}) ; "
              f"GnRH period ranges {per[T>T[-1]-4000].min():.1f}..{per[T>T[-1]-4000].max():.1f} min")
        m = T > T[-1] - 60*60          # last 60 h
        ax[k].plot(T[m]/60, E2[m], lw=0.9, color="C0")
        ax[k].set_ylabel(f"tau={tau:.0f}\nE2")
        ax[k].grid(alpha=0.3)
    ax[-1].set_xlabel("time (h)")
    fig.suptitle(f"Two-loop EFP (E2->GnRH freq/amp restored, cE2={cE2}): E2 vs delay")
    fig.tight_layout()
    fig.savefig("fig_two_loop_tau_sweep.png", dpi=130)
    print("wrote fig_two_loop_tau_sweep.png")
