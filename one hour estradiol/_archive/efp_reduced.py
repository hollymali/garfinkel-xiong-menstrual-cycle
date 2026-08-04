"""
Reduced early-follicular-phase (EFP) estradiol model.
=====================================================

Starting point for the "one hour estradiol" project.

We drop LH and P4 from the faithful Rasgon EFP model. This is exact, not an
approximation: in the EFP equations LH and P4 are pure one-way FOLLOWERS
(E2 -> LH -> P4) and neither ever appears in the dE2 or dFSH equations, so
they cannot influence estradiol. The only real feedback loop is E2 <-> FSH,
forced by the GnRH pulse clock.

State (3 variables):
    E2   estradiol
    FSH
    phi  GnRH pulse-clock phase (advances at 1/period(E2); pulses at integers)

GnRH is algebraic, not a state variable:
    GnRH = amplitude(E2) * pulse_kernel(phi)

Behaviour as printed: E2 relaxes to a stable fixed point (damped focus). The
E2-FSH Jacobian has trace = -(kE2+kFSH) < 0, so sustained / aperiodic
oscillation needs an added ingredient (delay, a 3rd slow variable, or external
forcing such as an HPT/temperature rhythm) -- which is what this project will
explore.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# GnRH pulse kernel: Gaussian bump each time phi crosses an integer.
# --------------------------------------------------------------------------
SIGMA_PHASE = 0.04   # pulse half-width as a fraction of the inter-pulse interval

def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIGMA_PHASE) ** 2)

# --------------------------------------------------------------------------
# Generic Hill helpers
#   up   (stimulate): A * x^n/(1+x^n) + C     -> C .. A+C
#   down (suppress) : (A + x^n)/(1+x^n) + C   -> A+C .. 1+C   (A>1)
# --------------------------------------------------------------------------
def hill_up(x, k, n, A, C):
    z = (k * x) ** n
    return A * z / (1.0 + z) + C

# --------------------------------------------------------------------------
# Parameters (EFP, faithful values from the whole-model work)
# --------------------------------------------------------------------------
EFP = dict(kE2=0.011, kFSH=0.009, cE2=0.025, cFSH=5e-3)

# E2-dependent GnRH pulse shape (piecewise-linear, from the paper's Appendix C)
def efp_period(E2):        # minutes between pulses (array-safe)
    return np.clip(0.267 * np.asarray(E2) - 6.75, 60.0, 100.0)

def efp_amplitude(E2):     # pulse height (array-safe)
    return np.clip(-0.013 * np.asarray(E2) + 8.2, 3.0, 5.0)

# Feedback functions
def fsh_to_e2(FSH):  return hill_up(FSH, 0.15, 7, 150.0, 250.0)   # 250 -> 400, stimulate
def gnrh_to_fsh(G):  return hill_up(G,   0.2, 15,  5.0,   5.0)    #   5 -> 10,   drive
def e2_to_fsh(E2):   return 9.0 / (1.0 + (0.0025 * E2) ** 15) + 1.0  # 10 -> 1, suppress

# --------------------------------------------------------------------------
# Right-hand side (E2, FSH, phi)
# --------------------------------------------------------------------------
def rhs(t, y, p=EFP):
    E2, FSH, phi = y
    GnRH = efp_amplitude(E2) * pulse_kernel(phi)

    dE2  = p['cE2']  * fsh_to_e2(FSH)                        - p['kE2']  * E2
    dFSH = p['cFSH'] * (e2_to_fsh(E2) + gnrh_to_fsh(GnRH))   - p['kFSH'] * FSH
    dphi = 1.0 / efp_period(E2)
    return [dE2, dFSH, dphi]


def run(t_end=8000.0, y0=None, max_step=0.25):
    if y0 is None:
        y0 = [200.0, 5.0, 0.0]   # E2, FSH, phi
    sol = solve_ivp(rhs, (0.0, t_end), y0, max_step=max_step,
                    rtol=1e-8, atol=1e-8, dense_output=True)
    return sol


if __name__ == "__main__":
    sol = run()
    t = sol.t
    E2, FSH, phi = sol.y
    GnRH = efp_amplitude(E2) * pulse_kernel(phi)

    print(f"final E2  = {E2[-1]:.3f}")
    print(f"final FSH = {FSH[-1]:.3f}")
    print(f"E2 range over last 2000 min: "
          f"{E2[t > t[-1]-2000].min():.3f} .. {E2[t > t[-1]-2000].max():.3f}")

    fig, ax = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(t/60, E2, lw=1.0);  ax[0].set_ylabel("E2")
    ax[1].plot(t/60, FSH, lw=1.0); ax[1].set_ylabel("FSH")
    ax[2].plot(t/60, GnRH, lw=0.6); ax[2].set_ylabel("GnRH"); ax[2].set_xlabel("time (h)")
    for a in ax: a.grid(alpha=0.3)
    fig.suptitle("Reduced EFP model (E2 <-> FSH, GnRH pulse-forced)")
    fig.tight_layout()
    fig.savefig("fig_efp_reduced.png", dpi=130)
    print("wrote fig_efp_reduced.png")
