"""
Forced + delayed reduced EFP model  --  "one hour estradiol"
=============================================================

Question: with GnRH held at a clean ~1-hour pulse period (physiological EFP),
can estradiol still show dynamics that are NOT just a boring 1:1 copy of the
pulse train?

Answer we're testing: yes -- if the E2<->FSH loop has internal dynamics of its
own. A periodically forced nonlinear oscillator can answer at a DIFFERENT,
longer, non-repeating timescale (subharmonics, slow envelopes, skipped beats).

Ingredient added here: a time DELAY tau in the E2->FSH negative-feedback arm,
   dFSH/dt uses e2_to_fsh( E2(t - tau) )
which is the Goodwin/Hopf mechanism the Rasgon paper names in prose but never
wrote down. tau=0 recovers the flat fixed point; large enough tau makes the
loop self-oscillate, and the 1h GnRH forcing then rides on / interacts with it.

GnRH is a FIXED-period pulse train (default 60 min) -- we are deliberately
keeping the ~1h period for now (physiological EFP; see sourced note).

Integrator: fixed-step RK4 with a history buffer for the delayed E2 (frozen
across each 0.5-min step -- negligible error at this dt).
"""

import numpy as np
import matplotlib.pyplot as plt

# ---- GnRH pulse kernel (Gaussian bump each period) -----------------------
SIGMA_PHASE = 0.04
def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIGMA_PHASE) ** 2)

# ---- Hill helpers --------------------------------------------------------
def hill_up(x, k, n, A, C):
    z = (k * x) ** n
    return A * z / (1.0 + z) + C

def fsh_to_e2(FSH):  return hill_up(FSH, 0.15, 7, 150.0, 250.0)   # stimulate
def gnrh_to_fsh(G):  return hill_up(G,   0.2, 15,  5.0,   5.0)    # drive

# E2->FSH suppression. K_E2 = the E2 midpoint of the switch. The ORIGINAL model
# used K_E2=400 (k=1/400=0.0025), but the resting E2 sits at ~571 -- out on the
# flat floor -- so the loop has no gain there. Aligning K_E2 with the resting E2
# puts the operating point on the STEEP flank, giving the loop gain a delay can
# destabilize.
N_E2 = 15
def e2_to_fsh(E2, K_E2=400.0):
    return 9.0 / (1.0 + (E2 / K_E2) ** N_E2) + 1.0  # suppress

# ---- Parameters (faithful EFP) -------------------------------------------
kE2, kFSH = 0.011, 0.009
cE2, cFSH = 0.025, 5e-3


def deriv(t, E2, FSH, E2_delayed, gnrh_period, gnrh_amp, K_E2):
    G = gnrh_amp * pulse_kernel(t / gnrh_period)
    dE2  = cE2  * fsh_to_e2(FSH)                                  - kE2  * E2
    dFSH = cFSH * (e2_to_fsh(E2_delayed, K_E2) + gnrh_to_fsh(G))  - kFSH * FSH
    return dE2, dFSH


def simulate(t_end=12000.0, dt=0.5, tau=0.0, gnrh_period=60.0, gnrh_amp=5.0,
             K_E2=400.0, y0=(200.0, 5.0)):
    n = int(round(t_end / dt))
    ndelay = int(round(tau / dt))
    T = np.arange(n) * dt
    E2 = np.empty(n); FSH = np.empty(n)
    E2[0], FSH[0] = y0
    hist0 = y0[0]  # assume E2 constant = y0 before t=0

    for i in range(n - 1):
        t, e, f = T[i], E2[i], FSH[i]
        j = i - ndelay
        ed = E2[j] if j >= 0 else hist0          # delayed E2 (frozen over step)
        k1 = deriv(t,        e,               f,               ed, gnrh_period, gnrh_amp, K_E2)
        k2 = deriv(t+dt/2,   e+dt/2*k1[0],    f+dt/2*k1[1],    ed, gnrh_period, gnrh_amp, K_E2)
        k3 = deriv(t+dt/2,   e+dt/2*k2[0],    f+dt/2*k2[1],    ed, gnrh_period, gnrh_amp, K_E2)
        k4 = deriv(t+dt,     e+dt*k3[0],      f+dt*k3[1],      ed, gnrh_period, gnrh_amp, K_E2)
        E2[i+1]  = e + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        FSH[i+1] = f + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return T, E2, FSH


def _summary(T, E2, tail_min=4000.0):
    m = T > T[-1] - tail_min
    e = E2[m]
    return e.min(), e.max(), e.max() - e.min()


if __name__ == "__main__":
    # Align the suppression midpoint with where E2 actually rests (~571) so the
    # loop has gain, then sweep the delay tau.
    K_E2 = 603.0
    taus = [0.0, 240.0, 360.0, 480.0, 600.0, 720.0]  # minutes
    fig, ax = plt.subplots(len(taus), 1, figsize=(11, 12), sharex=True)

    for k, tau in enumerate(taus):
        T, E2, FSH = simulate(tau=tau, K_E2=K_E2, t_end=16000.0)
        lo, hi, rng = _summary(T, E2)
        print(f"tau={tau:5.0f} min :  E2 tail range {lo:6.1f}..{hi:6.1f}  "
              f"(peak-to-peak {rng:6.1f})")
        # show the last 40 h so slow structure + 1h pulses are both visible
        m = T > T[-1] - 40*60
        ax[k].plot(T[m]/60, E2[m], lw=0.9, color="C3")
        ax[k].set_ylabel(f"tau={tau:.0f}\nE2")
        ax[k].grid(alpha=0.3)
    ax[-1].set_xlabel("time (h)")
    fig.suptitle(f"Reduced EFP forced at 60-min GnRH (K_E2={K_E2:.0f}): "
                 f"delay wakes up E2 dynamics")
    fig.tight_layout()
    fig.savefig("fig_efp_forced_tau_sweep.png", dpi=130)
    print("wrote fig_efp_forced_tau_sweep.png")
