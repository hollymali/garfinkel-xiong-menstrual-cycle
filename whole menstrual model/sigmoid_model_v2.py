"""
Redesigned EFP menstrual-cycle model (v2) -- the sigmoid redesign WITH the GnRH pulse generator.

This supersedes sigmoid_model.py (which is the pre-redesign linear-LH / additive-FSH / no-inhibin
version, kept as a documented reference). v2 implements the finalized design:

  * GnRH pulse generator  = impulsive ODE: a phase clock phi ticks at an E2-set frequency; each time
                            phi crosses an integer, GnRH gets an instantaneous kick (+A); between
                            kicks GnRH clears first-order (-k_GnRH*GnRH). -> spike-and-decay train.
  * Michaelis-Menten GnRH coupling  MM(GnRH; K_G) for BOTH gonadotropins (replaces linear LH).
  * Inhibin B  as a NEW state (the real FSH-selective follicular brake; E2 no longer brakes FSH).
  * Shared pituitary priming gain (1 + A*Hp(E2)) on both gonadotropins (E2 SENSITIZES, same sign).

ANCHORING (all EFP constants pinned; see the two project memories):
  * clearances k = ln2/t_half.
  * frequency sigmoid f_min/f_max/K_freq from human LH-pulse intervals (early-foll ~1/120 min floor,
    late-foll ~1/60 min ceiling).
  * K_G = 1.0 = the GnRH pulse PEAK (A=1): MM at half-max AT THE PEAK -> coupling tracks pulse
    amplitude with headroom both ways (NOT the unworkable <MM>=0.5 recipe, which forces K_G~1e-6 and
    saturates the coupling). <MM(GnRH;K_G=1)> ~= 0.025 over the ~120-min EFP cycle.
  * V_L=9.78, V_F=0.574 solved from the steady-state cycle-average budgets  <production> = k*<hormone>.

Native units: E2 pmol/L, FSH/LH IU/L, InhB pg/mL, P4 nmol/L, time min.
"""

import numpy as np


# ---- Hill primitives -----------------------------------------------------------------
def Hp(x, K, n):   # activation sigmoid  (0 -> 1)
    z = (x / K) ** n
    return z / (1.0 + z)

def Hm(x, K, n):   # inhibition sigmoid  (1 -> 0) = 1 - Hp
    return 1.0 / (1.0 + (x / K) ** n)

def MM(x, K):      # Michaelis-Menten (Hill n=1): saturating, non-cooperative
    return x / (K + x)


# ---- Parameters (EFP) ----------------------------------------------------------------
P = dict(
    # clearances (per min): k = ln2 / t_half
    k_E2  = np.log(2)/25,      # t1/2 25 min
    k_FSH = np.log(2)/240,     # t1/2 240 min (sialylated, slow)
    k_LH  = np.log(2)/90,      # t1/2 ~90 min EFFECTIVE (slow glycoform phase). NOT the 20-min fast
                               #   phase: LH clearance is bi-exponential; the SLOW terminal glycoform
                               #   (80-196 min, Bousfield 2022) governs the INTERPULSE tail. t1/2=20
                               #   made LH crash to ~0 between 120-min pulses (40x swing); the slow
                               #   phase sets a realistic nadir while pulses stay resolvable.
    k_P4  = np.log(2)/5,       # t1/2 5 min
    k_inh = np.log(2)/120,     # t1/2 ~120 min (SHAKY: n=2 patients)
    k_GnRH= 0.23,              # t1/2 ~3 min (firmly anchored)

    # GnRH generator
    A      = 1.0,              # A0: reference kick amplitude at EFP E2 (convention, absorbed into K_G)
    f_min  = 0.0083,           # /min  (1 pulse / 120 min: early-foll low-E2 floor)
    f_max  = 0.0167,           # /min  (1 pulse / 60 min:  late-foll high-E2 ceiling)
    K_freq = 500.0,            # pmol/L (E2 half-max of the follicular frequency rise; soft)
    n_freq = 4,
    # E2 -> pulse AMPLITUDE: DECREASING (follicular NEGATIVE feedback, NKB suppression).
    #   kick = A_max * Hm(E2;K_amp,n_amp)   (A_max computed below = A0/Hm(E2_ref) so kick=A0 at EFP E2,
    #   keeping '1'=EFP pulse and the K_G=peak anchor valid). Suppression deepens as E2 rises above EFP.
    #   This is E2's HYPOTHALAMIC negative feedback (distinct from FREQUENCY, which E2 raises, and from
    #   PITUITARY priming, which E2 raises).
    K_amp  = 315.0,            # pmol/L: E2 half-suppression (~Rasgon amp-law midpoint 8.2-0.013*E2; SOFT)
    n_amp  = 2,
    E2_ref = 140.0,            # EFP E2 where amplitude is normalized to A0

    # E2 <- FSH
    V_E   = 7.764, K_FSH_E2 = 6.5, n_FSH_E2 = 4,
    # InhB <- FSH  (structural clone of E2 loop)
    V_inh = 0.924, K_FSHinh = 6.5, n_FSHinh = 4,
    # shared E2 priming (off at follicular E2)
    K_sens = 734.0, n_sens = 8, A_LH = 4.0, A_FSH = 1.0,
    # LH <- MM(GnRH) * priming   (V_L re-solved below from k_LH to hold mean LH=7.1)
    V_L   = None, K_G = 1.0,
    # FSH <- basal + MM(GnRH) * priming * inhibin brake
    V_F   = 0.574, b_FSH = np.log(2)/240 * 4.0,   # nadir ~4 IU/L
    K_inh = 80.0, n_inh = 4,
    # P4 <- basal + LH  (V_P = 0 in EFP: no granulosa LH receptors yet)
    V_P   = 0.0, b_P4 = np.log(2)/5 * 0.21, K_LH_P4 = 7.0, n_LH_P4 = 4,
)

# EFP steady operating levels (initial conditions for the slow states)
EFP = dict(E2=140.0, InhB=80.0, FSH=6.5, LH=7.1, P4=0.21)

# Re-solve V_L from the steady-state cycle-average budget  V_L*<MM>*gain = k_LH*LH  (gain~=1 at EFP).
# <MM(GnRH;K_G=1)> = 0.0252 is fixed by the GENERATOR (A, K_G, k_GnRH, period), independent of k_LH,
# so changing the LH clearance just rescales V_L -- the mean stays pinned at 7.1.
MMBAR_EFP = 0.0252
P['V_L'] = P['k_LH'] * EFP['LH'] / MMBAR_EFP


# Max pulse amplitude (at E2 -> 0), as a MULTIPLE of the EFP pulse. Computed once so that the kick
# equals A0 at EFP E2 -- folds the old normalization (Hm(E2)/Hm(E2_ref)) into a single constant, so
# amplitude is just "A_max * one decreasing sigmoid" instead of a sigmoid divided by a sigmoid.
P['A_max'] = P['A'] / Hm(P['E2_ref'], P['K_amp'], P['n_amp'])   # = 1.198 (unsuppressed pulse = 1.2x EFP)


def amp_gain(E2, p=P):
    """E2 -> GnRH pulse amplitude (follicular NEGATIVE feedback, NKB): a DECREASING sigmoid.
    = A_max * Hm(E2;K_amp,n_amp), with A_max fixed so amplitude = A0 at EFP E2 (so '1' = an
    EFP-strength pulse and the K_G=peak anchor stays valid); shrinks as E2 rises above EFP."""
    return p['A_max'] * Hm(E2, p['K_amp'], p['n_amp'])


# ---- Continuous RHS (GnRH kick handled separately by the phi clock) -------------------
def rhs(y, p=P):
    """Derivatives of the 6 continuous states. Order: [E2, InhB, FSH, LH, P4, GnRH]."""
    E2, InhB, FSH, LH, P4, GnRH = y

    coupling = MM(GnRH, p['K_G'])
    prime_LH  = 1.0 + p['A_LH']  * Hp(E2, p['K_sens'], p['n_sens'])
    prime_FSH = 1.0 + p['A_FSH'] * Hp(E2, p['K_sens'], p['n_sens'])
    brake     = Hm(InhB, p['K_inh'], p['n_inh'])

    dE2   = p['V_E']  * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2'])   - p['k_E2']  * E2
    dInhB = p['V_inh']* Hp(FSH, p['K_FSHinh'], p['n_FSHinh'])   - p['k_inh'] * InhB
    dLH   = p['V_L'] * coupling * prime_LH                      - p['k_LH']  * LH
    dFSH  = p['b_FSH'] + p['V_F'] * coupling * prime_FSH * brake- p['k_FSH'] * FSH
    dP4   = p['b_P4']  + p['V_P'] * Hp(LH, p['K_LH_P4'], p['n_LH_P4']) - p['k_P4'] * P4
    dGnRH = -p['k_GnRH'] * GnRH
    return np.array([dE2, dInhB, dFSH, dLH, dP4, dGnRH])


def simulate(t_max=960.0, dt=0.01, p=P, y0=None):
    """Fixed-step RK4 on the continuous states + impulsive phi-clock kicks on GnRH."""
    if y0 is None:
        # seed GnRH with a pulse at t=0 (phi "just fired") to avoid a one-period startup dead-zone
        y0 = np.array([EFP['E2'], EFP['InhB'], EFP['FSH'], EFP['LH'], EFP['P4'], p['A']])
    n = int(round(t_max / dt))
    t = np.linspace(0.0, t_max, n + 1)
    Y = np.empty((n + 1, 6))
    Y[0] = y0
    y = y0.copy()
    phi = 0.0
    pulse_times = []
    for i in range(1, n + 1):
        # phase clock advances at the E2-set frequency; kick on integer crossing
        dphi = p['f_min'] + (p['f_max'] - p['f_min']) * Hp(y[0], p['K_freq'], p['n_freq'])
        phi_new = phi + dphi * dt
        if np.floor(phi_new) > np.floor(phi):
            y[5] += p['A'] * amp_gain(y[0], p)   # kick; amplitude E2-suppressed (foll neg feedback)
            pulse_times.append(t[i])
        phi = phi_new
        # RK4 for the continuous part
        k1 = rhs(y, p)
        k2 = rhs(y + 0.5*dt*k1, p)
        k3 = rhs(y + 0.5*dt*k2, p)
        k4 = rhs(y + dt*k3, p)
        y = y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        Y[i] = y
    return t, Y, np.array(pulse_times)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t, Y, pulses = simulate(t_max=960.0, dt=0.01)
    E2, InhB, FSH, LH, P4, GnRH = Y.T
    th = t / 60.0   # hours for the x-axis

    # cycle-average coupling check over the last full period
    period = 1.0 / (P['f_min'] + (P['f_max']-P['f_min'])*Hp(140.0, P['K_freq'], P['n_freq']))
    mask = t >= t[-1] - period
    MMbar = MM(GnRH[mask], P['K_G']).mean()
    print(f"EFP pulse period      = {period:.1f} min ({len(pulses)} pulses in {t[-1]/60:.1f} h)")
    print(f"<MM(GnRH;K_G)>        = {MMbar:.4f}   (anchor target ~0.025)")
    print(f"GnRH peak             = {GnRH.max():.3f}")
    print(f"means over run: E2={E2.mean():.1f}  InhB={InhB.mean():.1f}  FSH={FSH.mean():.2f}  "
          f"LH={LH.mean():.2f}  P4={P4.mean():.3f}")

    fig, ax = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    panels = [
        (ax[0,0], GnRH, "GnRH (kick+decay train)", "tab:purple", None, "GnRH (a.u., A=1)"),
        (ax[0,1], LH,   "LH (pulsatile, tracks GnRH)", "tab:red", EFP['LH'], "LH (IU/L)"),
        (ax[1,0], FSH,  "FSH (slow, inhibin-braked)", "tab:blue", EFP['FSH'], "FSH (IU/L)"),
        (ax[1,1], E2,   "E2 (FSH-driven)", "tab:green", EFP['E2'], "E2 (pmol/L)"),
        (ax[2,0], InhB, "Inhibin B", "tab:orange", EFP['InhB'], "Inhibin B (pg/mL)"),
        (ax[2,1], P4,   "P4 (basal only, V_P=0 in EFP)", "tab:brown", EFP['P4'], "P4 (nmol/L)"),
    ]
    for a, series, title, c, level, ylabel in panels:
        a.plot(th, series, color=c, lw=0.9)
        if level is not None:
            a.axhline(level, color="k", ls=":", lw=0.8, alpha=0.6, label=f"EFP level {level:g}")
            a.legend(loc="upper right", fontsize=8)
        a.set_title(title, fontsize=10)
        a.set_ylabel(ylabel, fontsize=9)
        a.grid(alpha=0.25)
    for a in ax[2]:
        a.set_xlabel("time (hours)")
    fig.suptitle("Redesigned EFP model v2 -- GnRH generator + MM coupling + inhibin", fontsize=12)
    fig.tight_layout()
    out = "fig_sigmoid_v2_timeseries.png"
    fig.savefig(out, dpi=130)
    print(f"saved {out}")
