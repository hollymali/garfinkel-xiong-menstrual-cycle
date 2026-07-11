"""
Autonomous menstrual-cycle model (v6) -- GnRH collapsed to a constant.

v5 drove LH & FSH through an exact cycle-average of the fast GnRH pulse train
(the `coupling(E2,P4)` function: E2-set frequency, P4 clock-brake, E2-modulated
amplitude, closed-form saturable-receptor mean).  A direct experiment showed the
autonomous monthly rhythm SURVIVES when that whole apparatus is frozen to a single
constant: the clock lives in the slow ovarian relaxation loop and the pituitary's
own biphasic E2 response, not in the GnRH pulses.  GnRH is a PERMISSIVE CARRIER.

v6 therefore deletes the entire GnRH generator -- the `coupling()` function and its
~13 parameters (A, f_min/f_max/K_freq, K_amp/E2_ref, A_sg, K_G, K_P4f, k_GnRH) -- and
replaces the GnRH drive with one constant, C_GnRH.  This is deliberately CRUDE: it
throws away the pulse frequency, amplitude, and P4-brake modulation.  It is the
minimal honest statement of the model: if the rhythm still emerges (it does), the
rhythm cannot depend on the structure of the GnRH pulse train.

What the crude constant costs (all quantitative, none structural), v5 -> v6:
    period 29.4 -> ~25 d (lost the P4->GnRH frequency brake that slowed the clock)
    LH surge 62 -> ~27   (lost GnRH's amplitude surge arm; surge now fires purely
                          from the pituitary gain gL, i.e. E2 acting on the pituitary)
    luteal LH floor 0.7 -> ~3 (a flat drive can't quiet LH between cycles)

Two further changes vs the first v6 draft (verified to hold the ~25 d cycle):
  - dropped the (Fmax-F) logistic cap on follicle growth (barely engaged: F peaks ~2.4
    of a 10 cap); rescaled k_grow 4e-5 -> 3.2e-4 to keep the follicular-phase pacing.
  - the FSH-independent LH->E2 term is now gated on follicle MATURITY, Hp(F;K_Fmat),
    representing LH-receptor acquisition (FSH/E2-driven, primes in the late follicular
    phase BEFORE the surge) rather than on instantaneous LH level.  Rlh=Hp(F;K_Fmat) is
    algebraic (no new state variable); separates receptor priming from surge activation.

Native units: E2 pmol/L, FSH/LH IU/L, InhB pg/mL, P4 nmol/L, F & Lut dimensionless, time min.
"""

import numpy as np

np.seterr(all="ignore")   # sigmoids at x=0 raise harmless divide/overflow warnings


# ---- Hill primitives -----------------------------------------------------------------
def Hp(x, K, n):   # activation sigmoid (0 -> 1)
    z = (x / K) ** n
    return z / (1.0 + z)

def Hm(x, K, n):   # inhibition sigmoid (1 -> 0)
    return 1.0 / (1.0 + (x / K) ** n)


DAY = 1440.0   # minutes per day


# ---- Parameters ----------------------------------------------------------------------
P = dict(
    # clearances (per min): k = ln2 / t_half
    k_E2  = np.log(2)/25,
    k_FSH = np.log(2)/240,
    k_LH  = np.log(2)/90,
    k_P4  = np.log(2)/5,
    k_inh = np.log(2)/120,

    # GnRH drive: a single constant permissive carrier (was the coupling(E2,P4) average).
    # Value = v5's follicular operating point coupling(150,1) ~ 0.0247.
    C_GnRH = 0.0247,

    # shared biphasic pituitary gain on LH & FSH: E2 negative arm + positive/surge arm
    B_LH = 0.7, B_FSH = 0.7, K_neg = 200.0, n_neg = 3,
    K_sens = 734.0, n_sens = 8, A_LH = 4.0, A_FSH = 1.0,

    # E2 production: FSH-driven (FSH-saturated dominant follicle) + LH-gated two-cell + luteal(Lut)
    V_E    = 10.4, K_FSH_E2 = 3.0, n_FSH_E2 = 4,
    V_E2L  = 3.7,  K_LHe2 = 12.0, n_LHe2 = 3,
    K_Fmat = 1.2,  n_Fmat = 4,                                     # LH-receptor availability: gate LH->E2 on follicle maturity F
    V_EL_lut = 6.0,                                                # corpus-luteum E2 (luteal hump)
    s0 = 0.7, K_LH_E2 = 7.0, n_LH_E2 = 2,                          # LH->E2 substrate gate

    # Inhibin B from the follicle (F); FSH-selective brake on FSH production
    V_inh = 2.0, K_FSHinh = 6.5, n_FSHinh = 4, K_inh = 80.0, n_inh = 4,

    # gonadotropins
    V_L = 2.34, V_F = 0.9, b_FSH = np.log(2)/240 * 3.0,
    K_lutFSH = 1.0, n_lutFSH = 2,                                  # luteal FSH brake (inhibin-A-like, Lut)

    # P4 from the corpus luteum
    b_P4     = np.log(2)/5 * 0.5,
    kP4_prod = 3.4,                                               # luteal P4 peaks ~40 nmol/L

    # slow ovarian states: follicle F, corpus luteum Lut
    k_grow = 3.2e-4, K_Frec = 3.5, n_Frec = 4,                   # FSH-driven recruitment (logistic cap removed)
    k_atr  = 5.0e-5,                                              # baseline follicle loss
    K_P4gate = 2.5, n_P4gate = 4,                                 # low-P4 gate on follicle growth
    k_ov   = 5.0e-3, K_surge = 33.0, n_surge = 6,                 # LH surge ovulates the follicle
    k_lut  = 5.0e-3, k_luteo = np.log(2)/(4.1*DAY),               # CL formation / regression (~14 d luteal)
)

# Default initial condition = menses (low E2, small follicle, no corpus luteum)
Y0 = np.array([150.0, 90.0, 7.0, 6.0, 1.0, 0.5, 0.0])
STATE_NAMES = ["E2", "InhB", "FSH", "LH", "P4", "F", "Lut"]


def rhs(y, p=P):
    """Derivatives of the 7 slow states.  Order: [E2, InhB, FSH, LH, P4, F, Lut].

    GnRH no longer appears: its drive is the constant p['C_GnRH'] (`c` below)."""
    E2, InhB, FSH, LH, P4, F, Lut = y

    c     = p['C_GnRH']                                 # constant permissive GnRH drive
    gL    = max(1 - p['B_LH']  * Hp(E2, p['K_neg'], p['n_neg'])
                  + p['A_LH']  * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    gF    = max(1 - p['B_FSH'] * Hp(E2, p['K_neg'], p['n_neg'])
                  + p['A_FSH'] * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    LH_andro_level = p['s0'] + (1 - p['s0']) * Hp(LH, p['K_LH_E2'], p['n_LH_E2'])   # LH-driven theca androgen substrate
    LH_surge       = Hp(LH, p['K_surge'], p['n_surge'])          # ~0 except during the LH surge
    Rlh            = Hp(F,  p['K_Fmat'],  p['n_Fmat'])           # LH-receptor availability (primes with follicle maturity, before the surge)

    dE2  = ( p['V_E']  * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2']) * LH_andro_level
           + p['V_E2L'] * Rlh * Hp(LH,  p['K_LHe2'],  p['n_LHe2']) ) * F \
           + p['V_EL_lut'] * Lut                                          - p['k_E2']  * E2
    dInhB= p['V_inh'] * Hp(FSH, p['K_FSHinh'], p['n_FSHinh']) * F          - p['k_inh'] * InhB
    dLH  = p['V_L'] * c * gL                                              - p['k_LH']  * LH
    dFSH = ( p['b_FSH'] + p['V_F'] * c * gF
             * Hm(InhB, p['K_inh'], p['n_inh'])
             * Hm(Lut,  p['K_lutFSH'], p['n_lutFSH']) )                   - p['k_FSH'] * FSH
    dP4  = p['b_P4'] + p['kP4_prod'] * Lut                                - p['k_P4']  * P4
    dF   = ( p['k_grow'] * Hp(FSH, p['K_Frec'], p['n_Frec'])
             * Hm(P4, p['K_P4gate'], p['n_P4gate'])
           - p['k_atr'] * F
           - p['k_ov']  * LH_surge * F )
    dLut = p['k_lut'] * LH_surge * F                                     - p['k_luteo'] * Lut
    return np.array([dE2, dInhB, dFSH, dLH, dP4, dF, dLut])


def simulate(days=200.0, dt=1.0, p=P, y0=None):
    """Fixed-step RK4 on the 7 slow states (mean-field; no fast-pulse integration)."""
    if y0 is None:
        y0 = Y0.copy()
    n = int(days * DAY / dt)
    t = np.linspace(0.0, days * DAY, n + 1)
    Y = np.empty((n + 1, 7))
    Y[0] = y0
    y = y0.copy()
    for i in range(1, n + 1):
        k1 = rhs(y, p)
        k2 = rhs(y + 0.5*dt*k1, p)
        k3 = rhs(y + 0.5*dt*k2, p)
        k4 = rhs(y + dt*k3, p)
        y = np.maximum(y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4), 0.0)
        Y[i] = y
    return t, Y


def surge_peaks(t, Y, burn_days=60.0):
    """Indices of the preovulatory E2 peaks after a burn-in (one per cycle)."""
    td = t / DAY
    E2 = Y[:, 0]
    dE = np.diff(E2)
    return [i for i in range(1, len(E2) - 1)
            if dE[i-1] > 0 and dE[i] <= 0 and td[i] > burn_days and E2[i] > 500]


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t, Y = simulate(days=200.0, dt=1.0)
    td = t / DAY
    pk = surge_peaks(t, Y)
    period = td[pk[-1]] - td[pk[-2]]

    print(f"autonomous cycle (constant GnRH): {len(pk)} ovulations; period = {period:.1f} days")
    print("min/max over run:")
    for j, nm in enumerate(STATE_NAMES):
        print(f"  {nm:5s}: {Y[:, j].min():8.2f} .. {Y[:, j].max():8.2f}")

    a = pk[-2]
    seg0 = slice(a, pk[-1] + 1)
    imin = a + int(np.argmin(Y[seg0, 0]))
    seg = slice(imin, imin + int(period * DAY) + 1)
    x = td[seg] - td[imin]
    Ys = Y[seg]
    ovd = x[int(np.argmax(Ys[:, 0]))]

    fig, ax = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    panels = [("E2 (pmol/L)", 0, "tab:green"), ("LH (IU/L)", 3, "tab:red"),
              ("FSH (IU/L)", 2, "tab:blue"),  ("P4 (nmol/L)", 4, "tab:brown"),
              ("Inhibin B (pg/mL)", 1, "tab:orange"), ("Follicle / Corpus luteum", 5, None)]
    for a_, (lab, idx, c) in zip(ax.ravel(), panels):
        if idx == 5:
            a_.plot(x, Ys[:, 5], color="tab:purple", lw=1.6, label="Follicle F")
            a_.plot(x, Ys[:, 6], color="tab:cyan",   lw=1.6, label="Corpus luteum Lut")
            a_.legend(fontsize=8)
        else:
            a_.plot(x, Ys[:, idx], color=c, lw=1.8)
        a_.axvline(ovd, color="gray", ls="--", lw=0.8, alpha=0.6)
        a_.set_title(lab, fontsize=10); a_.set_ylabel(lab, fontsize=9); a_.grid(alpha=0.25)
    for a_ in ax[2]:
        a_.set_xlabel("cycle day (0 = menses; dashed = ovulation)")
    fig.suptitle(f"Autonomous menstrual cycle (v6) -- constant GnRH drive, period {period:.1f} d",
                 fontsize=12)
    fig.tight_layout()
    out = "fig_sigmoid_v6_cycle.png"
    fig.savefig(out, dpi=130, facecolor="white")
    print(f"saved {out}  (ovulation at cycle day {ovd:.1f})")
