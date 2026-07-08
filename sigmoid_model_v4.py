"""
Autonomous menstrual-cycle model (v4) -- the SLOW RELAXATION / POSITIVE-FEEDBACK loop.

WHY v4 exists
-------------
v2 (sigmoid_model_v2.py) is a GnRH-pulse-FORCED oscillator pinned at the early-follicular
phase (EFP): it holds a stable slow fixed point and the only oscillation is the fast pulse
train.  v3 added negative-feedback terms (E2 biphasic pituitary arm, LH->E2 gate, P4->clock
brake).  A full bifurcation hunt then PROVED that negative feedback -- any channel, gain,
steepness or timescale -- only re-tunes the forced pulse train and NEVER produces an
autonomous rhythm (finding-p4gnrh-no-hopf).  The conclusion was that a ~monthly cycle needs a
different ingredient: a slow relaxation / positive-feedback loop.  v4 builds exactly that and
DOES self-oscillate.

The two ingredients v4 adds
---------------------------
1. Two slow ovarian states that make the system EXCITABLE / relaxation-type:
     F   = dominant-follicle E2-secreting capacity.  Grows under FSH over ~2 weeks; the
           preovulatory LH surge dumps it (ovulation).
     Lut = corpus-luteum mass.  Formed FROM the ovulated follicle by the surge; secretes
           progesterone; regresses slowly (luteolysis) over ~2 weeks.
2. The preovulatory POSITIVE-feedback engine that fires the surge:
     * the dominant follicle is FSH-SATURATED (low K_FSH_E2), so its E2 output tracks
       follicle mass F and keeps climbing through the mid-follicular FSH decline instead of
       collapsing with it -- E2 can therefore reach the surge threshold;
     * a two-cell LH->E2 term (V_E2L): the mature follicle converts LH drive to E2, closing
       the fast E2<->LH positive loop that makes the surge explosive;
     * an E2 POSITIVE surge arm on the GnRH pulse amplitude (A_sg, AVPV/kisspeptin), so the
       hypothalamic coupling does not stay suppressed by the follicular NKB negative arm at
       high E2.

The relaxation cycle
--------------------
  follicular : FSH recruits F -> E2 = (FSH- and LH-driven) * F climbs -> E2 crosses the
               positive-feedback threshold -> explosive LH SURGE
  ovulation  : the surge ovulates F (F -> Lut)
  luteal     : Lut -> high P4 -> P4 brakes the GnRH clock AND atreses any new follicle
               (refractory period); E2 falls because the follicle is gone
  reset      : the corpus luteum regresses -> P4 falls -> the brake releases -> FSH rises
               (the intercycle FSH rise) and recruits the next follicle -> repeat.

EFP is no longer a pinned fixed point; it is one PHASE the limit cycle passes through.

Fast pulses -> exact mean field
-------------------------------
The month-long slow dynamics are integrated with the fast GnRH pulse train replaced by its
EXACT cycle-average coupling.  For impulsive kicks of amplitude a at frequency f decaying at
rate k_GnRH, the steady-train mean of MM(GnRH;K_G) is closed-form:

    <MM> = (f / k_GnRH) * ln( (K_G + G_pk) / (K_G + G_pk - a) ),   G_pk = a / (1 - e^{-k/f})

This evaluates to 0.0252 at the EFP operating point -- identical to the pinned MMBAR of the
pulse model -- so it is the honest cycle-average, NOT a moving-average low-pass (which would
inject the filtering-artifact "beat" documented in finding-p4gnrh-no-hopf).

Native units: E2 pmol/L, FSH/LH IU/L, InhB pg/mL, P4 nmol/L, F and Lut dimensionless, time min.
"""

import numpy as np

np.seterr(all="ignore")   # sigmoids with x=0 raise harmless divide/overflow warnings


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

    # GnRH generator (enters only through the mean-field coupling below)
    A      = 1.0,
    f_min  = 0.0083, f_max = 0.0167, K_freq = 500.0, n_freq = 4,   # E2 raises pulse FREQUENCY
    K_amp  = 315.0, n_amp = 2, E2_ref = 140.0,                     # E2 lowers pulse AMPLITUDE (NKB)
    A_sg   = 3.0,                                                  # E2 POSITIVE surge arm (AVPV/kiss)
    K_G    = 1.0,
    K_P4f  = 25.0, n_P4f = 2,                                      # P4 brakes the whole clock

    # shared biphasic pituitary gain on LH & FSH (v3): E2 negative arm + positive/surge arm
    B_LH = 0.7, B_FSH = 0.7, K_neg = 200.0, n_neg = 3,
    K_sens = 734.0, n_sens = 8, A_LH = 4.0, A_FSH = 1.0,

    # E2 production: FSH-driven (dominant follicle FSH-SATURATED: low K so E2 tracks F, not the
    # mid-follicular FSH decline) + two-cell LH-driven arm; all scaled by follicle mass F.
    V_E   = 4.77, K_FSH_E2 = 3.0, n_FSH_E2 = 4,
    V_E2L = 0.9,  K_LHe2  = 12.0, n_LHe2  = 3,
    s0 = 0.7, K_LH_E2 = 7.0, n_LH_E2 = 2,                          # LH->E2 substrate gate (v3)

    # Inhibin B (granulosa product, scaled by F): FSH-selective brake
    V_inh = 0.924, K_FSHinh = 6.5, n_FSHinh = 4, K_inh = 80.0, n_inh = 4,

    # gonadotropins  (V_L, V_F, b_FSH carried over from the v3 EFP anchoring)
    V_L = 2.505, V_F = 0.672, b_FSH = np.log(2)/240 * 4.0,

    # ---- P4 from the corpus luteum ----
    b_P4    = np.log(2)/5 * 0.21,   # tiny follicular basal
    kP4_prod= 2.1,                  # P4 secretion per unit Lut (luteal P4 peaks ~30-35 nmol/L)

    # ---- slow ovarian states: follicle F, corpus luteum Lut ----
    k_grow = 5.5e-5, K_Frec = 6.0, n_Frec = 4, Fmax = 10.0,   # FSH-driven follicle growth (~2 wk)
    k_atr  = 8.0e-4, K_P4F = 22.0, n_P4F = 4,                 # luteal P4 atreses new follicle (refractory)
    k_ov   = 5.0e-3, K_surge = 22.0, n_surge = 6,             # LH surge ovulates the follicle (fast)
    k_lut  = 5.0e-3,                                          # surge converts follicle F -> luteal Lut
    k_luteo= np.log(2)/(3.6*DAY),                             # corpus-luteum regression t1/2 ~3.6 d
)

# Reference EFP operating levels (also the default initial condition, mid-follicular follicle).
EFP = dict(E2=140.0, InhB=80.0, FSH=6.5, LH=7.1, P4=0.21, F=1.0, Lut=0.0)

STATE_NAMES = ["E2", "InhB", "FSH", "LH", "P4", "F", "Lut"]


def coupling(E2, P4, p=P):
    """EXACT cycle-average <MM(GnRH;K_G)> of the impulsive kick train -- no pulses integrated.

    f = E2-set frequency * P4 clock-brake ;  a = E2-modulated pulse amplitude (NKB negative arm
    + AVPV positive surge arm).  Returns the mean coupling the slow gonadotropin ODEs integrate.
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

    c    = coupling(E2, P4, p)
    gL   = max(1 - p['B_LH']  * Hp(E2, p['K_neg'], p['n_neg'])
                 + p['A_LH']  * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    gF   = max(1 - p['B_FSH'] * Hp(E2, p['K_neg'], p['n_neg'])
                 + p['A_FSH'] * Hp(E2, p['K_sens'], p['n_sens']), 0.0)
    gate = p['s0'] + (1 - p['s0']) * Hp(LH, p['K_LH_E2'], p['n_LH_E2'])
    surge = Hp(LH, p['K_surge'], p['n_surge'])          # ~0 except during the LH surge

    dE2  = ( p['V_E']  * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2']) * gate
           + p['V_E2L'] * Hp(LH,  p['K_LHe2'],  p['n_LHe2']) ) * F        - p['k_E2']  * E2
    dInhB= p['V_inh'] * Hp(FSH, p['K_FSHinh'], p['n_FSHinh']) * F          - p['k_inh'] * InhB
    dLH  = p['V_L'] * c * gL                                              - p['k_LH']  * LH
    dFSH = p['b_FSH'] + p['V_F'] * c * gF * Hm(InhB, p['K_inh'], p['n_inh'])- p['k_FSH'] * FSH
    dP4  = p['b_P4'] + p['kP4_prod'] * Lut                                - p['k_P4']  * P4
    dF   = ( p['k_grow'] * Hp(FSH, p['K_Frec'], p['n_Frec']) * (p['Fmax'] - F)
           - p['k_atr']  * Hp(P4, p['K_P4F'], p['n_P4F']) * F
           - p['k_ov']   * surge * F )
    dLut = p['k_lut'] * surge * F                                        - p['k_luteo'] * Lut
    return np.array([dE2, dInhB, dFSH, dLH, dP4, dF, dLut])


def simulate(days=150.0, dt=1.0, p=P, y0=None):
    """Fixed-step RK4 on the 7 slow states (mean-field; no fast-pulse integration)."""
    if y0 is None:
        y0 = np.array([EFP['E2'], EFP['InhB'], EFP['FSH'], EFP['LH'], EFP['P4'], EFP['F'], EFP['Lut']])
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
        y = y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        y = np.maximum(y, 0.0)
        Y[i] = y
    return t, Y


def cycle_period(t, Y, burn_days=30.0, thresh=300.0):
    """Cycle length = spacing of E2 preovulatory peaks after a burn-in (returns list of days)."""
    td = t / DAY
    E2 = Y[:, 0]
    dE = np.diff(E2)
    peaks = [i for i in range(1, len(E2) - 1)
             if dE[i-1] > 0 and dE[i] <= 0 and td[i] > burn_days and E2[i] > thresh]
    return [td[peaks[j+1]] - td[peaks[j]] for j in range(len(peaks) - 1)], [td[i] for i in peaks]


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t, Y = simulate(days=150.0, dt=1.0)
    td = t / DAY

    print("min/max over 150-day run:")
    for j, nm in enumerate(STATE_NAMES):
        print(f"  {nm:5s}: min={Y[:, j].min():8.2f}  max={Y[:, j].max():8.2f}")

    periods, peak_days = cycle_period(t, Y)
    if periods:
        print(f"\nautonomous cycle: {len(peak_days)} ovulatory E2 peaks; "
              f"period(s) = {[round(x, 1) for x in periods]} days")
    else:
        print("\nno autonomous cycle detected")

    # ---- figure: two representative cycles ----
    per = periods[-1] if periods else 30.0
    w = (td >= 45) & (td <= 45 + 2*per)
    fig, ax = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    series = [("E2 (pmol/L)", 0, "tab:green"), ("LH (IU/L)", 3, "tab:red"),
              ("FSH (IU/L)", 2, "tab:blue"),  ("P4 (nmol/L)", 4, "tab:brown"),
              ("Follicle F", 5, "tab:purple"),("Corpus luteum Lut", 6, "tab:orange")]
    for a, (lab, idx, c) in zip(ax.ravel(), series):
        a.plot(td[w] - 45, Y[w, idx], color=c, lw=1.6)
        a.set_title(lab, fontsize=10); a.set_ylabel(lab, fontsize=9); a.grid(alpha=0.25)
    for a in ax[2]:
        a.set_xlabel("time (days)")
    fig.suptitle("Autonomous menstrual cycle (v4) -- follicle-maturation relaxation loop, "
                 f"period ~{per:.0f} d", fontsize=12)
    fig.tight_layout()
    out = "fig_sigmoid_v4_autonomous_cycle.png"
    fig.savefig(out, dpi=130, facecolor="white")
    print(f"saved {out}")
