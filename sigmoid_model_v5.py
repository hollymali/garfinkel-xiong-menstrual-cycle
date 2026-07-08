"""
Autonomous menstrual-cycle model (v5) -- clinically calibrated 28-day cycle.

v4 produced a stable autonomous ~monthly rhythm from a slow relaxation / positive-feedback loop,
but its hormone LEVELS were off and its phase structure was scrambled (no menstrual nadir; early-
follicular E2 ~360 vs clinical ~150).  v5 keeps the v4 relaxation-loop idea and re-engineers the
ovarian / luteal block so the time series and levels match clinical reference data across the cycle.

Structural changes over v4
--------------------------
1. Follicle growth is GATED on LOW progesterone (Hm(P4)).  The follicle cannot grow while the
   corpus luteum is producing P4, so it stays ~0 through the luteal phase and only regrows after
   luteolysis.  This creates a genuine MENSTRUAL NADIR where E2 and P4 are both low at once.
2. E2 has TWO sources:  the follicle F (follicular rise + explosive preovulatory peak) AND the
   corpus luteum Lut (the luteal secondary hump).  -> biphasic E2 like the real cycle.
3. FSH gets a LUTEAL brake from Lut (inhibin-A-like), so FSH is low in the luteal phase and rises
   at menses (the intercycle FSH rise that recruits the next follicle).
4. The dominant follicle is FSH-SATURATED for BOTH its E2 output (low K_FSH_E2) and its GROWTH
   (low K_Frec): once selected it keeps maturing through the mid-follicular FSH decline instead of
   stalling -- required for E2 to reach the surge threshold.
5. Two-cell LH->E2 gain (V_E2L) is LH-gated, so it lifts the preovulatory E2 SPIKE (LH ~60) but
   not the follicular baseline (LH ~6).

The relaxation cycle
--------------------
  menses     : F small, Lut regressed, P4 & E2 low  ->  FSH intercycle rise
  follicular : FSH recruits F (P4 low, growth gate open) -> E2 climbs to the surge threshold
  ovulation  : E2 positive feedback fires the LH SURGE -> the surge ovulates F (F -> Lut)
  luteal     : Lut -> high P4 (brakes clock, closes the follicle growth gate) + secondary E2 hump;
               Lut -> luteal FSH brake
  reset      : corpus luteum regresses -> P4 falls -> growth gate reopens, FSH rises -> next follicle

Fast GnRH pulses are replaced by their EXACT closed-form cycle-average coupling (see coupling());
this is the honest mean, not a moving-average low-pass, so the month-long run is artifact-free.

Native units: E2 pmol/L, FSH/LH IU/L, InhB pg/mL, P4 nmol/L, F & Lut dimensionless, time min.

Calibration (period 27.6 d, ovulation day ~14) vs clinical reference ranges:
    E2 menses 135 (110-180) | E2 preovulatory 985 (700-1500) | LH base 5-6 (2-10) |
    LH surge 61 (25-100) | FSH intercycle 7.3 (5-9) | FSH luteal 3.2 (2-4) |
    P4 follicular 1-5 (<3) | P4 mid-luteal 43 (20-70) | InhB 13-141 (20-150)
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
    k_GnRH= 0.23,

    # GnRH generator (enters only through the mean-field coupling)
    A      = 1.0,
    f_min  = 0.0083, f_max = 0.0167, K_freq = 500.0, n_freq = 4,   # E2 raises pulse FREQUENCY
    K_amp  = 315.0, n_amp = 2, E2_ref = 140.0,                     # E2 lowers pulse AMPLITUDE (NKB)
    A_sg   = 2.5,                                                  # E2 POSITIVE surge arm (AVPV/kiss)
    K_G    = 1.0,
    K_P4f  = 25.0, n_P4f = 2,                                      # P4 brakes the whole clock

    # shared biphasic pituitary gain on LH & FSH: E2 negative arm + positive/surge arm
    B_LH = 0.7, B_FSH = 0.7, K_neg = 200.0, n_neg = 3,
    K_sens = 734.0, n_sens = 8, A_LH = 4.0, A_FSH = 1.0,

    # E2 production: FSH-driven (FSH-saturated dominant follicle) + LH-gated two-cell + luteal(Lut)
    V_E    = 10.4, K_FSH_E2 = 3.0, n_FSH_E2 = 4,
    V_E2L  = 3.2,  K_LHe2 = 12.0, n_LHe2 = 3,
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
    k_grow = 4.0e-5, K_Frec = 3.5, n_Frec = 4, Fmax = 10.0,       # FSH-saturated recruitment
    k_atr  = 5.0e-5,                                              # baseline follicle loss
    K_P4gate = 2.5, n_P4gate = 4,                                 # low-P4 gate on follicle growth
    k_ov   = 5.0e-3, K_surge = 33.0, n_surge = 6,                 # LH surge ovulates the follicle
    k_lut  = 5.0e-3, k_luteo = np.log(2)/(4.1*DAY),               # CL formation / regression (~14 d luteal)
)

# Default initial condition = menses (low E2, small follicle, no corpus luteum)
Y0 = np.array([150.0, 90.0, 7.0, 6.0, 1.0, 0.5, 0.0])
STATE_NAMES = ["E2", "InhB", "FSH", "LH", "P4", "F", "Lut"]


def coupling(E2, P4, p=P):
    """EXACT cycle-average <MM(GnRH;K_G)> of the impulsive GnRH kick train.

    f = E2-set frequency * P4 clock-brake ;  a = E2-modulated pulse amplitude (NKB negative arm +
    AVPV positive surge arm).  Steady-train mean of MM(GnRH;K_G) is closed-form:
        <MM> = (f/k) * ln((K_G+G_pk)/(K_G+G_pk-a)),  G_pk = a/(1-e^{-k/f}).
    Equals 0.0252 at the EFP operating point (matches the pinned MMBAR of the pulse model).
    """
    f = (p['f_min'] + (p['f_max'] - p['f_min']) * Hp(E2, p['K_freq'], p['n_freq'])) \
        * Hm(P4, p['K_P4f'], p['n_P4f'])
    f = max(f, 1e-9)
    A_max = p['A'] / Hm(p['E2_ref'], p['K_amp'], p['n_amp'])
    a = A_max * (Hm(E2, p['K_amp'], p['n_amp']) + p['A_sg'] * Hp(E2, p['K_sens'], p['n_sens']))
    Gpk = a / (1.0 - np.exp(-p['k_GnRH'] / f))
    return (f / p['k_GnRH']) * np.log((p['K_G'] + Gpk) / (p['K_G'] + Gpk - a))


def rhs(y, p=P):
    """Derivatives of the 7 slow states.  Order: [E2, InhB, FSH, LH, P4, F, Lut]."""
    E2, InhB, FSH, LH, P4, F, Lut = y

    c     = coupling(E2, P4, p)
    gL    = max(1 - p['B_LH']  * Hp(E2, p['K_neg'], p['n_neg'])
                  + p['A_LH']  * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    gF    = max(1 - p['B_FSH'] * Hp(E2, p['K_neg'], p['n_neg'])
                  + p['A_FSH'] * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    gate  = p['s0'] + (1 - p['s0']) * Hp(LH, p['K_LH_E2'], p['n_LH_E2'])
    surge = Hp(LH, p['K_surge'], p['n_surge'])          # ~0 except during the LH surge

    dE2  = ( p['V_E']  * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2']) * gate
           + p['V_E2L'] * Hp(LH,  p['K_LHe2'],  p['n_LHe2']) ) * F \
           + p['V_EL_lut'] * Lut                                          - p['k_E2']  * E2
    dInhB= p['V_inh'] * Hp(FSH, p['K_FSHinh'], p['n_FSHinh']) * F          - p['k_inh'] * InhB
    dLH  = p['V_L'] * c * gL                                              - p['k_LH']  * LH
    dFSH = ( p['b_FSH'] + p['V_F'] * c * gF
             * Hm(InhB, p['K_inh'], p['n_inh'])
             * Hm(Lut,  p['K_lutFSH'], p['n_lutFSH']) )                   - p['k_FSH'] * FSH
    dP4  = p['b_P4'] + p['kP4_prod'] * Lut                                - p['k_P4']  * P4
    dF   = ( p['k_grow'] * Hp(FSH, p['K_Frec'], p['n_Frec'])
             * Hm(P4, p['K_P4gate'], p['n_P4gate']) * (p['Fmax'] - F)
           - p['k_atr'] * F
           - p['k_ov']  * surge * F )
    dLut = p['k_lut'] * surge * F                                        - p['k_luteo'] * Lut
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

    print(f"autonomous cycle: {len(pk)} ovulations; period = {period:.1f} days")
    print("min/max over run:")
    for j, nm in enumerate(STATE_NAMES):
        print(f"  {nm:5s}: {Y[:, j].min():8.2f} .. {Y[:, j].max():8.2f}")

    # one cycle aligned to the menstrual E2 nadir (cycle day 0 = menses)
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
    fig.suptitle(f"Autonomous menstrual cycle (v5) -- clinically calibrated, period {period:.1f} d",
                 fontsize=12)
    fig.tight_layout()
    out = "fig_sigmoid_v5_cycle.png"
    fig.savefig(out, dpi=130, facecolor="white")
    print(f"saved {out}  (ovulation at cycle day {ovd:.1f})")
