"""
Reference reimplementation of the Rasgon et al. (2003) HPG-axis menstrual-cycle model.

Rasgon NL, Pumphrey L, Prolo P, Elman S, Negrao AB, Licinio J, Garfinkel A.
"Emergent Oscillations in Mathematical Model of the Human Menstrual Cycle."
CNS Spectrums. 2003;8(11):805-814.

This is a SCRATCH / VERIFICATION file -- the goal is to get the model running and see
whether it reproduces the paper's qualitative behaviour. The clean teaching version goes
into menstrual_cycle.ipynb afterwards.

All equations/parameters are transcribed in MODEL_EXTRACTED.md (with the four resolved
ambiguity notes). Two design choices that modernize the original Berkeley Madonna model:

  1. Integration: scipy.integrate.solve_ivp (adaptive RK45) instead of fixed-step RK4.
  2. GnRH pulse train: the original used Berkeley Madonna's PULSE() (a DT-scaled spike).
     We replace it with a step-size-independent forcing function built from an internal
     PHASE variable phi that advances at rate 1/Period(state). A pulse is a narrow bump
     of PEAK height = Amplitude(state) centred at each integer phase. Peak-height (not
     area) is the right reading: the GnRH->FSH sigmoid in Appendix C is plotted over
     GnRH in 0-6, i.e. it expects the pulse *height* to equal the Amplitude variable
     (~3-5), not an impulse area.
"""

import numpy as np
from scipy.integrate import solve_ivp


# --------------------------------------------------------------------------------------
# GnRH pulse kernel
# --------------------------------------------------------------------------------------
# A periodic bump with PEAK 1 at every integer value of the phase phi.
# d = signed distance to the nearest integer (in cycles); Gaussian of width SIGMA_PHASE.
SIGMA_PHASE = 0.04   # pulse half-width as a fraction of the inter-pulse interval

def pulse_kernel(phi):
    d = phi - np.round(phi)
    return np.exp(-0.5 * (d / SIGMA_PHASE) ** 2)


# --------------------------------------------------------------------------------------
# Generic Hill helpers (see NOTE 2 in MODEL_EXTRACTED.md)
#   negative feedback (decreasing): (A + x^n)/(1 + x^n) + C   with A>1   -> A+C down to 1+C
#   positive feedback (increasing): A * x^n/(1 + x^n) + C                -> C    up to   A+C
# --------------------------------------------------------------------------------------
def hill_down(x, k, n, A, C):
    z = (k * x) ** n
    return (A + z) / (1.0 + z) + C

def hill_up(x, k, n, A, C):
    z = (k * x) ** n
    return A * z / (1.0 + z) + C


# ======================================================================================
# EARLY FOLLICULAR PHASE (EFP)
# ======================================================================================
EFP = dict(kE2=0.011, kFSH=0.009, kLH=0.002, kP4=9.2e-5,
           cE2=0.01, cFSH=2e-7, cLH=0.1, cP4=0.045)   # cFSH = the suspicious 2e-7


def efp_multiplier(E2):
    if E2 < 128:   return 2.0
    if E2 > 800:   return 0.5
    return -0.002 * E2 + 2.26

def efp_period(E2):
    if E2 < 250:   return 60.0
    if E2 > 400:   return 100.0
    return 0.267 * E2 - 6.75

def efp_amplitude(E2):
    if E2 < 250:   return 5.0
    if E2 > 400:   return 3.0
    return -0.013 * E2 + 8.2

# Feedback functions (EFP) -- forms confirmed against Appendix C plots
def efp_fsh_to_e2(FSH):   return hill_up(FSH, 0.15, 7,  150.0, 250.0)   # 250 -> 400
def efp_gnrh_to_fsh(G):   return hill_up(G,   0.2, 15,   5.0,   5.0)    #   5 -> 10
def efp_e2_to_fsh(E2):    return 9.0 / (1.0 + (0.0025 * E2) ** 15) + 1.0  # 10 -> 1, NOTE 1


def efp_rhs(t, y, p=EFP, ovx=False):
    E2, FSH, LH, P4, phi = y
    # Ovariectomy (Fig 2): remove ovarian steroids -> their feedback on GnRH/FSH is gone.
    E2_fb = 0.0 if ovx else E2
    amp = efp_amplitude(E2_fb)
    GnRH = amp * pulse_kernel(phi)

    dE2  = 0.0 if ovx else p['cE2'] * efp_fsh_to_e2(FSH) - p['kE2'] * E2
    dFSH = p['cFSH'] * (efp_e2_to_fsh(E2_fb) + efp_gnrh_to_fsh(GnRH))   - p['kFSH'] * FSH
    dLH  = p['cLH']  * efp_multiplier(E2_fb) * GnRH                     - p['kLH']  * LH
    dP4  = 0.0 if ovx else p['cP4'] * LH - p['kP4'] * P4
    dphi = 1.0 / efp_period(E2_fb)
    return [dE2, dFSH, dLH, dP4, dphi]


# ======================================================================================
# MID-LUTEAL PHASE (MLP)
# ======================================================================================
MLP = dict(kE2=0.9, kFSH=0.4, kLH=1.9, kP4=0.0835,
           cE2=0.005, cFSH=0.05, cLH=0.0014, cP4=0.08)


def mlp_period(E2, P4):
    s = E2 + P4
    if s < 282.55:  return 120.0
    if s > 293.3:   return 240.0
    return 11.111 * s - 3018.88889

def mlp_amplitude(E2, P4):
    s = E2 + P4
    if s < 275:     return 10.0
    if s > 310:     return 1.0
    return -0.2571 * s + 80.7143

# Feedback functions (MLP) -- checked against Appendix D plots
def mlp_e2_to_fsh(E2):  return hill_down(E2, 0.36, 30, 9.5, 0.5)         # 10 -> 1.5
def mlp_e2_to_lh(E2):   return hill_down(E2, 0.36, 30, 5.0, 5.0)         # 10 -> 6
def mlp_fsh_to_e2(FSH): return hill_up(FSH, 0.2, 8, 0.75, 2.5)          # 2.5 -> 3.25 (NOTE 2: *)
def mlp_p4_to_lh(P4):   return hill_down(P4, 0.00347, 175, 5.0, 5.0)     # 10 -> 6  (NOTE 3: +5)


def mlp_rhs(t, y, p=MLP):
    E2, FSH, LH, P4, phi1, phi2, phi3 = y
    amp = mlp_amplitude(E2, P4)
    # self-priming: primary train + two satellite trains (intervals period, +20, +30)
    GnRH = amp * (pulse_kernel(phi1) + pulse_kernel(phi2) + pulse_kernel(phi3))

    dE2  = p['cE2']  * (mlp_fsh_to_e2(FSH) - p['kE2'] * E2)
    dFSH = p['cFSH'] * (5.0 * GnRH + mlp_e2_to_fsh(E2)                  - p['kFSH'] * FSH)
    dLH  = p['cLH']  * (20.0 * GnRH + mlp_e2_to_lh(E2) + mlp_p4_to_lh(P4) - p['kLH'] * LH)
    dP4  = p['cP4']  * (3.0 * LH                                       - p['kP4'] * P4)

    per = mlp_period(E2, P4)
    return [dE2, dFSH, dLH, dP4, 1.0 / per, 1.0 / (per + 20.0), 1.0 / (per + 30.0)]


# ======================================================================================
# Driver
# ======================================================================================
def run(phase='EFP', t_end=8000.0, max_step=0.25, y0=None):
    if phase == 'EFP':
        if y0 is None:
            y0 = [200.0, 5.0, 10.0, 5.0, 0.0]
        rhs = efp_rhs
    else:
        if y0 is None:
            y0 = [250.0, 5.0, 10.0, 40.0, 0.0, 0.0, 0.0]
        rhs = mlp_rhs

    t_eval = np.arange(0.0, t_end, 0.5)
    sol = solve_ivp(rhs, (0.0, t_end), y0, method='RK45',
                    t_eval=t_eval, max_step=max_step, rtol=1e-6, atol=1e-9)
    return sol


def run_rk4(phase='EFP', t_end=8000.0, dt=0.25, y0=None, p=None):
    """Fixed-step classic RK4 -- matches the original Berkeley Madonna scheme
    (RK4, DT=0.25 min). Sticks to the printed equations; only the integrator
    changes vs run() (which uses adaptive solve_ivp).

    Returns an object with .t and .y (shape [nvars, npts]) so it is a drop-in
    for the solve_ivp solution used by the plotting code.
    """
    if phase == 'EFP':
        if y0 is None:
            y0 = [200.0, 5.0, 10.0, 5.0, 0.0]
        rhs = (lambda t, y: efp_rhs(t, y, p)) if p is not None else efp_rhs
    else:
        if y0 is None:
            y0 = [250.0, 5.0, 10.0, 40.0, 0.0, 0.0, 0.0]
        rhs = (lambda t, y: mlp_rhs(t, y, p)) if p is not None else mlp_rhs

    n = int(round(t_end / dt))
    y = np.array(y0, dtype=float)
    ys = np.empty((n + 1, y.size))
    ts = np.empty(n + 1)
    ys[0] = y
    ts[0] = 0.0
    for i in range(n):
        t = i * dt
        k1 = np.asarray(rhs(t, y), dtype=float)
        k2 = np.asarray(rhs(t + 0.5 * dt, y + 0.5 * dt * k1), dtype=float)
        k3 = np.asarray(rhs(t + 0.5 * dt, y + 0.5 * dt * k2), dtype=float)
        k4 = np.asarray(rhs(t + dt,       y + dt * k3),       dtype=float)
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        ys[i + 1] = y
        ts[i + 1] = t + dt

    class _Sol:
        pass
    sol = _Sol()
    sol.t = ts
    sol.y = ys.T
    sol.success = True
    return sol


def gnrh_trace(phase, sol):
    """Recompute the GnRH forcing along a solution (for plotting)."""
    if phase == 'EFP':
        amp = np.array([efp_amplitude(e) for e in sol.y[0]])
        return amp * pulse_kernel(sol.y[4])
    amp = np.array([mlp_amplitude(e, p) for e, p in zip(sol.y[0], sol.y[3])])
    return amp * (pulse_kernel(sol.y[4]) + pulse_kernel(sol.y[5]) + pulse_kernel(sol.y[6]))


if __name__ == '__main__':
    for phase in ('EFP', 'MLP'):
        sol = run(phase)
        names = ['E2', 'FSH', 'LH', 'P4']
        print(f'\n===== {phase} =====  success={sol.success}  npts={sol.t.size}')
        for i, nm in enumerate(names):
            y = sol.y[i]
            print(f'  {nm:3s}  min={y.min():.4g}  max={y.max():.4g}  '
                  f'last={y[-1]:.4g}  range={y.max()-y.min():.4g}')
        g = gnrh_trace(phase, sol)
        print(f'  GnRH min={g.min():.3g} max={g.max():.3g}')
