"""
Redesigned "more natural" sigmoid version of the Rasgon (2003) menstrual-cycle model.

This is the PARALLEL model to menstrual_model.py (the faithful reference). Here every
feedback is a proper saturating Hill function instead of the paper's linear/piecewise
forms, and the constants are anchored to physiology rather than tuned for oscillation.

DESIGN PRINCIPLE
    production = saturating Hill feedback(s)  -  first-order clearance
    (No basal term: b=0 for all hormones. Dropped for simplicity -- E2/LH/P4 basal is
    negligible/off-axis, and though FSH IS constitutively secreted (a real GnRH-independent
    floor), we have no clean number; reintroduce b_FSH at the GnRH/budget-split step if needed.)
    Combination rule chosen from biology:
      * ADDITIVE        when two regulators act independently        -> FSH (GnRH + E2)
      * MULTIPLICATIVE  when one regulator sets the GAIN on another  -> LH (E2 PRIMES the pituitary),
                        and the optional two-cell E2 term (LH gates FSH-driven aromatase)
    E2 -> LH is POSITIVE pituitary self-priming: rising/sustained E2 SENSITIZES the gonadotrope, so
    LH-per-GnRH INCREASES with E2 (the surge is priming taken to high E2). This replaces an earlier
    "brake" (a negative pituitary gain), which was biologically misplaced. E2's follicular NEGATIVE
    feedback is NOT a pituitary gain cut: it is HYPOTHALAMIC -- E2 lowers GnRH PULSE AMPLITUDE (and
    slows the pulse generator) -> that belongs in the GnRH generator (pending), entering via vl_gnrh.

CALIBRATION (see solve_amplitudes() -- it reproduces every amplitude below)
    * clearances k = ln2 / t_half        (half-lives from literature; NOT the paper's values)
    * thresholds K                       (anchored to physiological hormone levels)
    * amplitudes V, b                    (steady-state pinning: production = k*X at quasi-
                                          steady phase points -- follicular for EFP, luteal MLP)
    * exponents n                        (FREE -- receptor cooperativity ~2-4; surge steep ~8)

PENDING GnRH (Holly builds the GnRH generator separately): everything multiplied by GnRH
(V_GnRH_FSH, the V_E2_FSH split, V_L) is only identifiable once a GnRH unit convention is fixed --
GnRH and its gain are an entangled product. The steady-state PRODUCTION BUDGETS are solved
here (FSH 0.0144, LH 0.246 IU/L/min at EFP) ready to be split when GnRH is defined.

KEY GOTCHAS
    * The paper's clearances are NOT physiological (k_P4 9.2e-5 => 5-DAY P4 half-life vs real
      ~5 min; real LH clears ~17x faster). Real k's change behaviour a lot -> re-tune around them.
    * The model is genuinely PHASE-SPECIFIC: luteal E2/P4 come from the corpus luteum, not from
      their drivers, so one parameter set can't span follicular->luteal. Calibrate EFP to
      follicular levels, MLP to luteal levels (note V_P jumps 0.044 -> 7.49, ~170x, = luteal mass).
    * The surge needs SUSTAINED E2 (~50 h), not instantaneous -- a pure Hill on E2 may fire too
      eagerly; a duration mechanism may be needed later.

Data sources: Elecsys cycle study (PMC8042396, E2/LH/P4 ranges); FSH ranges (Medscape 2089048);
half-lives (JCEM 2022, PMC9731043); E2 surge threshold 200 pg/mL = 734 pmol/L
(Frontiers Neurosci 2022, 10.3389/fnins.2022.953252).
"""

import numpy as np

# --------------------------------------------------------------------------------------
# Hill primitives (K = half-saturation point; H+ rises 0->1, H- falls 1->0, both =1/2 at x=K)
# --------------------------------------------------------------------------------------
def Hp(x, K, n):   # activation
    z = (x / K) ** n
    return z / (1.0 + z)

def Hm(x, K, n):   # inhibition  (= 1 - Hp)
    return 1.0 / (1.0 + (x / K) ** n)


# --------------------------------------------------------------------------------------
# Physiological inputs (native units: E2 pmol/L, FSH/LH IU/L, P4 nmol/L, time min)
# --------------------------------------------------------------------------------------
# Median hormone levels per phase (use quasi-steady phases to pin amplitudes; NOT ovulation,
# which is violently non-steady).
#   NOTE on 'foll': this model is the EARLY follicular phase (EFP), so E2 is anchored to the EARLY-
#   follicular level (~140; day-3 band ~80-200 pmol/L), NOT a whole-follicular median (~198, which is
#   really mid-follicular). FSH/LH/P4 EFP levels (5, 7.1, 0.21) are already early-follicular-consistent.
LEVELS = dict(
    foll=dict(E2=140, FSH=6.5, LH=7.1, P4=0.21),   # EFP early follicular: E2 was 198(mid-foll), FSH was 5
    ovul=dict(E2=757, FSH=10, LH=22.6, P4=1.81),   # for surge cross-checks only
    lut =dict(E2=412, FSH=4,  LH=6.2,  P4=28.8),
)

# Plasma half-lives (min) -> clearances k = ln2/t_half. NOT the paper's (its k_P4 => 5-DAY P4 t1/2).
# Chosen for the cycle-dynamics timescale: the RATIO of half-lives matters more than absolutes, and the
# essential separation is FSH slow (sialylated) >> LH fast (pulses must stay resolvable) > E2 > P4.
#   ---------------------------------------------------------------------------------------------------
#   hormone | t_half | k=ln2/t  | source for value used               | alternative / note
#   --------|--------|----------|-------------------------------------|----------------------------------
#   E2      |  25    | 0.0277   | IV estradiol terminal t1/2          | measured 27.45+-5.65 min (n=8) [3];
#           |        |          |   27.45 min, n=8 women [3]          |   our 25 ~ this directly measured value
#   FSH     | 240    | 0.0029   | canonical clinical ~4h (textbook)  | glycoform 343-757 min [1]; long
#           |        |          |                                     |   t1/2 = sialylation [2]
#   LH      |  20    | 0.0347   | classic fast-phase ~20 min (bulk    | glycoform terminal 80-196 min [1][2];
#           |        |          |   LH; fast clearance => pulses      |   NOT used: would smear LH pulses
#           |        |          |   stay resolvable)                  |
#   P4      |   5    | 0.1386   | rapid clearance: MCR 2100-2800      | IV t1/2 spans 3-90 min across studies;
#           |        |          |   L/day [4] -> short t1/2; 5 = fast | 5 = fast-phase end (high MCR confirms)
#   ---------------------------------------------------------------------------------------------------
#   SOURCES
#   [1] Bousfield GR, et al. "Determination of Half-lives of Circulating FSH and LH Glycoforms in Women
#       During GnRH Receptor Blockade." J Clin Endocrinol Metab. 2022;107(10):e4058-e4070.
#       https://academic.oup.com/jcem/article/107/10/e4058/6653068  (PMC9731043)
#   [2] Wide L, Eriksson K, et al. "Serum Half-Life of Pituitary Gonadotropins Is Decreased by
#       Sulfonation and Increased by Sialylation in Women." J Clin Endocrinol Metab. 2009;94(3):958-964.
#       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2681272/  (PubMed 19116233)
#   [3] White CM, et al. "The Pharmacokinetics of Intravenous Estradiol - A Preliminary Study."
#       Pharmacotherapy. 1998;18. IV bolus, 8 postmenopausal women; terminal t1/2 27.45 +- 5.65 min.
#       DOI 10.1002/j.1875-9114.1998.tb03157.x  (pagination not verified - abstract paywalled)
#   [4] Little B, Tait JF, Tait SA, Erlenmeyer F. "The metabolic clearance rate of progesterone in
#       males and ovariectomized females." J Clin Invest. 1966;45(6):901-912. DOI 10.1172/JCI105405.
#       (MCR 2100-2800 L/day => rapid clearance; IV t1/2 reported 3-90 min across the literature)
#   ---------------------------------------------------------------------------------------------------
HALF_LIVES_MIN = dict(E2=25, FSH=240, LH=20, P4=5)

def clearances():
    """k = ln2 / t_half  (per minute)."""
    return {h: np.log(2) / t for h, t in HALF_LIVES_MIN.items()}

# Thresholds K -- anchored to physiological levels (these were already ~right in the paper)
K_FSH_E2 = 6.5     # FSH -> E2     half-max FSH (IU/L) = EFP FSH median (day-3 ~3.5-12.5, healthy ~6-7).
                   #   Re-anchored from the paper's inherited 6 to the real EARLY-follicular FSH median
                   #   so E2 production is half-max at the EFP operating FSH (symmetric responsiveness).
                   #   Still the softest threshold: E2 actually rises while FSH FALLS through follicular
                   #   (follicle gains aromatase) -- the lumped FSH->E2 Hill can't capture that decoupling.
K_E2_FSH = 150.0   # E2  -> FSH    EFP brake midpoint (pmol/L), centered in the EARLY-follicular E2
                   #               band (~80-200, median ~140): the brake's steep point sits in the EFP
                   #               operating range so E2 negative feedback on FSH is LIVE during EFP.
                   #   WHY 150 not 200: 200 was anchored to a follicular MEDIAN of 198, but 198 is
                   #   mid-follicular; true EARLY-follicular E2 centers ~140 (day-3 ~80-200). K=150
                   #   re-centers the brake on the actual EFP excursion (half-max mid-EFP).
                   #   WHY NOT 400: K~400 is the WHOLE-CYCLE midpoint (E2 100->757) -> brake parked on
                   #   its flat top (~0.94) through EFP -> no EFP FSH decline. One sigmoid can't be
                   #   steep across both spans -> the phase-specific trade.
                   #   BIOLOGY NOTE: real EFP FSH decline is driven largely by INHIBIN B, which this
                   #   model lacks; with E2 the only FSH suppressor, E2 must carry the decline, forcing
                   #   K down into the EFP band. If inhibin is added later, K could rise back toward 400.
# E2 -> LH = PITUITARY SELF-PRIMING (POSITIVE). Rising/sustained E2 sensitizes the gonadotrope, so
# the gain of GnRH on LH RISES with E2 (Veldhuis JCI: priming absent without E2, unmasked over days;
# Evans & Veldhuis: E2 doubles LH burst duration).
# gain = G_LH0 + A_SENSITIVITY * H+(E2; K_E2_LH_sensitivity, n_E2_LH_sensitivity).
#   WHY NOT A BRAKE: a down-in-E2 pituitary gain is the WRONG SIGN/NODE. E2's follicular NEGATIVE
#   feedback is HYPOTHALAMIC -- it DECREASES GnRH PULSE AMPLITUDE (and lengthens the period); that
#   enters through vl_gnrh/GnRH (the pulse generator, pending), NOT as a pituitary gain reduction.
G_LH0        = 1.0     # baseline (unprimed) pituitary gain = the +C OFFSET on the priming sigmoid:
                       # gain = G_LH0 + A_sensitivity*H+(E2) runs G_LH0 -> G_LH0+A_sensitivity (1 -> 5), so it
                       # FLOORS AT 1, not 0. Unlike production sigmoids (V*H+, ->0 with no driver),
                       # the pituitary keeps a baseline GnRH responsiveness even unprimed; E2 only
                       # AMPLIFIES it. C=1 makes the gain a pure multiplier on V_L*GnRH.
K_E2_LH_sensitivity = 734.0   # priming/surge E2 threshold = 200 pg/mL (pmol/L): E2 sensitizes the pituitary
K_LH_P4      = 7.1     # LH  -> P4     half-max LH (IU/L) = follicular LH median 7.14 (PMC8042396),
                       #   half-max at the EFP LH operating point (consistent with K_FSH_E2 anchoring).
                       #   CAVEAT: fundamentally a LUTEAL parameter. In EFP the LH->P4 coupling is weak
                       #   (granulosa lack LH receptors until post-ovulatory luteinization) and EFP V_P
                       #   (0.044) makes at most ~0.32 nmol/L P4 even at saturating LH (vs ovulation
                       #   1.81) -> EFP does NOT constrain this K. Re-calibrate with MLP, where LH->P4
                       #   is the strong dominant coupling.

# Exponents n -- FREE parameters (tune via sensitivity analysis)
N_DEFAULT = 4
N_SENSITIVITY   = 8          # steep: priming/surge is switch-like once E2 crosses threshold

# Priming amplification -- FREE knob (start 3-5): gain rises to (G_LH0 + A_SENSITIVITY) at high E2
A_SENSITIVITY = 4.0
# TODO(duration): real priming needs E2 SUSTAINED ~days; this instantaneous Hill fires on LEVEL only
# and will prime too eagerly. Add a sustained-E2 / integrating gate when modeling the EFP->surge exit.


# --------------------------------------------------------------------------------------
# Steady-state amplitude solver:  production = k*X  =>  V = k*X / Hill(driver)   (basal b=0)
# Re-run this whenever you change a K or n and it re-pins V_E / V_P for you.
# --------------------------------------------------------------------------------------
def solve_amplitudes(nF=N_DEFAULT, nP=N_DEFAULT):
    k = clearances()
    f, l = LEVELS['foll'], LEVELS['lut']

    # E2 driven by FSH (follicular pin); no basal
    V_E = k['E2'] * f['E2'] / Hp(f['FSH'], K_FSH_E2, nF)

    # P4 driven by LH; no basal.  EFP pinned at follicular, MLP at luteal (corpus-luteum mass).
    V_P_efp = k['P4'] * f['P4'] / Hp(f['LH'], K_LH_P4, nP)
    V_P_mlp = k['P4'] * l['P4'] / Hp(l['LH'], K_LH_P4, nP)

    # GnRH-coupled production BUDGETS (solved totals; split once GnRH convention is set)
    budget_FSH = k['FSH'] * f['FSH']                                   # = V_GnRH_FSH*H+(GnRH) + V_E2_FSH*H-(E2)
    # EFP placeholder amplitude: treat follicular FSH production as fully E2-suppressible (GnRH
    # floor ~0 in follicular) and PIN it so production = k*FSH at the follicular steady state for
    # whatever K_E2_FSH is set. This keeps FSH balanced while the brake genuinely sweeps with E2.
    V_E2_FSH_eff = budget_FSH / Hm(f['E2'], K_E2_FSH, N_DEFAULT)
    gain_foll  = G_LH0 + A_SENSITIVITY * Hp(f['E2'], K_E2_LH_sensitivity, N_SENSITIVITY)  # priming ~0 at follicular E2
    budget_LH  = k['LH'] * f['LH']                                     # = V_L*<GnRH>*gain
    VL_times_GnRH = budget_LH / gain_foll

    return dict(V_E=V_E, V_P_efp=V_P_efp, V_P_mlp=V_P_mlp,
                budget_FSH=budget_FSH, V_E2_FSH_eff=V_E2_FSH_eff, VL_times_GnRH=VL_times_GnRH,
                **{f'k_{h}': v for h, v in k.items()})


# Solved parameter set (call solve_amplitudes() to re-derive; values below are the frozen result).
_A = solve_amplitudes()

PARAMS_EFP = dict(
    # clearances (per min)
    k_E2=_A['k_E2'], k_FSH=_A['k_FSH'], k_LH=_A['k_LH'], k_P4=_A['k_P4'],
    # E2:  V_E*H+(FSH; K, n)
    V_E=_A['V_E'], K_FSH_E2=K_FSH_E2, n_FSH_E2=N_DEFAULT,
    # FSH: V_GnRH_FSH*H+(GnRH) + V_E2_FSH*H-(E2; K, n)   -- V_GnRH_FSH/V_E2_FSH PENDING GnRH; budget solved.
    #   NOTE: basal b_FSH dropped for simplicity (no clean number), but FSH IS constitutively secreted
    #   (GnRH-independent floor; ovariectomy+GnRH-antiserum: LH vanishes, FSH persists). Reintroduce a
    #   small b_FSH from a GnRH-blockade residual when splitting budget_FSH at the GnRH step if needed.
    V_GnRH_FSH=None, V_E2_FSH=None, budget_FSH=_A['budget_FSH'], V_E2_FSH_eff=_A['V_E2_FSH_eff'],
    K_E2_FSH=K_E2_FSH, n_E2_FSH=N_DEFAULT,
    # LH:  V_L * GnRH * ( G_LH0 + A_sensitivity*H+(E2; K_E2_LH_sensitivity, n_E2_LH_sensitivity) )
    #      V_L and GnRH multiply EXPLICITLY in the RHS; only their PRODUCT V_L*<GnRH> is identifiable
    #      until the GnRH convention is set, so VL_times_GnRH stands in. Pituitary self-priming is
    #      POSITIVE; E2's NEGATIVE feedback (GnRH pulse amplitude) lives in the GnRH generator.
    V_L=None, VL_times_GnRH=_A['VL_times_GnRH'],
    G_LH0=G_LH0, K_E2_LH_sensitivity=K_E2_LH_sensitivity, n_E2_LH_sensitivity=N_SENSITIVITY, A_sensitivity=A_SENSITIVITY,
    # P4:  V_P*H+(LH; K, n)
    V_P=_A['V_P_efp'], K_LH_P4=K_LH_P4, n_LH_P4=N_DEFAULT,
)

PARAMS_MLP = dict(PARAMS_EFP, V_P=_A['V_P_mlp'])   # luteal P4 amplitude (corpus-luteum mass)


# --------------------------------------------------------------------------------------
# Redesigned EFP right-hand side.
# GnRH is passed in as a value (Holly's generator, built later). The LH carrier is V_L * GnRH:
# once V_L is set and a real GnRH(t) is passed, the two multiply explicitly below. Until then the
# solved PRODUCT VL_times_GnRH stands in so the LH equation balances at the follicular steady state.
# --------------------------------------------------------------------------------------
def efp_sigmoid_rhs(t, y, GnRH, p=PARAMS_EFP, vl_gnrh=None, v_gnrh_fsh=None, v_e2_fsh=None):
    E2, FSH, LH, P4 = y

    # E2 <- FSH (positive, saturating)
    dE2 = p['V_E'] * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2']) - p['k_E2'] * E2

    # FSH <- GnRH (+) and E2 (-), ADDITIVE.  Needs the GnRH split; until then approximate the
    # GnRH-driven part by its budget remainder so FSH still balances at steady state.
    fsh_neg = Hm(E2, p['K_E2_FSH'], p['n_E2_FSH'])
    if v_gnrh_fsh is not None and v_e2_fsh is not None:
        dFSH = v_gnrh_fsh * Hp(GnRH, 1.0, 1) + v_e2_fsh * fsh_neg - p['k_FSH'] * FSH
    else:
        dFSH = p['V_E2_FSH_eff'] * fsh_neg - p['k_FSH'] * FSH          # placeholder: pinned E2-suppressible budget

    # LH <- carrier (V_L * GnRH)  x  PITUITARY SELF-PRIMING gain (POSITIVE in E2: more LH per GnRH
    # as E2 rises). E2's follicular NEGATIVE feedback is NOT here -- it acts at the hypothalamus by
    # DECREASING GnRH PULSE AMPLITUDE (and slowing the generator), entering through GnRH itself.
    #   carrier = V_L * GnRH  (explicit product once both are available; else the solved stand-in)
    if vl_gnrh is not None:
        carrier = vl_gnrh                              # caller passed the product directly
    elif p['V_L'] is not None and GnRH is not None:
        carrier = p['V_L'] * GnRH                      # EXPLICIT V_L * GnRH (post-GnRH convention)
    else:
        carrier = p['VL_times_GnRH']                   # entangled-product stand-in until GnRH set
    gain = p['G_LH0'] + p['A_sensitivity'] * Hp(E2, p['K_E2_LH_sensitivity'], p['n_E2_LH_sensitivity'])
    dLH = carrier * gain - p['k_LH'] * LH

    # P4 <- LH (positive, saturating)
    dP4 = p['V_P'] * Hp(LH, p['K_LH_P4'], p['n_LH_P4']) - p['k_P4'] * P4

    return [dE2, dFSH, dLH, dP4]


if __name__ == "__main__":
    a = solve_amplitudes()
    print("clearances /min:", {h: round(a[f'k_{h}'], 4) for h in HALF_LIVES_MIN})
    print(f"V_E       = {a['V_E']:.2f} pmol/L/min")
    print(f"V_P (EFP) = {a['V_P_efp']:.4f}   V_P (MLP) = {a['V_P_mlp']:.3f} nmol/L/min")
    print(f"FSH budget (foll) = {a['budget_FSH']:.4f}   LH V_L*<GnRH> = {a['VL_times_GnRH']:.4f} IU/L/min")
    print("surge threshold 200 pg/mL =", round(200 / 272.4 * 1000), "pmol/L")
