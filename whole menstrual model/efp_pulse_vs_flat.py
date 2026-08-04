"""EFP submodel time series: GnRH pulses vs a flat GnRH of the same mean.

EFP (the follicular operating point, sigmoid_model_v2.py) is a stable EQUILIBRIUM, not a
cycle.  So the pulsed-vs-flat comparison here shows: the impulsive GnRH train makes the
gonadotropins (LH especially) RIPPLE around the EFP level, while a flat GnRH carrier of the
same cycle-average coupling drives them straight to the same operating point with no ripple.
Same equilibrium either way -> at the follicular set point, GnRH pulse structure is carrier,
not clock (the EFP analog of the whole-cycle result).

Left  = pulsatile GnRH (the v2 generator).   Right = flat constant GnRH (matched mean coupling).
Run from the repo root:  python efp_pulse_vs_flat.py  ->  fig_efp_pulse_vs_flat.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigmoid_model_v2 as m

P, Hp, Hm, MM, EFP = m.P, m.Hp, m.Hm, m.MM, m.EFP
T_MAX, DT = 1440.0, 0.01   # 24 h


def rhs_flat(y, coupling, p=P):
    """Continuous hormone RHS with GnRH coupling held at a fixed constant (no pulses)."""
    E2, InhB, FSH, LH, P4 = y
    prime_LH  = 1.0 + p['A_LH']  * Hp(E2, p['K_sens'], p['n_sens'])
    prime_FSH = 1.0 + p['A_FSH'] * Hp(E2, p['K_sens'], p['n_sens'])
    brake     = Hm(InhB, p['K_inh'], p['n_inh'])
    dE2   = p['V_E']  * Hp(FSH, p['K_FSH_E2'], p['n_FSH_E2'])          - p['k_E2']  * E2
    dInhB = p['V_inh']* Hp(FSH, p['K_FSHinh'], p['n_FSHinh'])          - p['k_inh'] * InhB
    dLH   = p['V_L'] * coupling * prime_LH                            - p['k_LH']  * LH
    dFSH  = p['b_FSH'] + p['V_F'] * coupling * prime_FSH * brake       - p['k_FSH'] * FSH
    dP4   = p['b_P4']  + p['V_P'] * Hp(LH, p['K_LH_P4'], p['n_LH_P4']) - p['k_P4']  * P4
    return np.array([dE2, dInhB, dFSH, dLH, dP4])


def simulate_flat(coupling, t_max=T_MAX, dt=DT):
    n = int(round(t_max / dt))
    t = np.linspace(0.0, t_max, n + 1)
    Y = np.empty((n + 1, 5))
    y = np.array([EFP['E2'], EFP['InhB'], EFP['FSH'], EFP['LH'], EFP['P4']])
    Y[0] = y
    for i in range(1, n + 1):
        k1 = rhs_flat(y, coupling); k2 = rhs_flat(y + 0.5*dt*k1, coupling)
        k3 = rhs_flat(y + 0.5*dt*k2, coupling); k4 = rhs_flat(y + dt*k3, coupling)
        y = y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        Y[i] = y
    return t, Y


if __name__ == "__main__":
    # --- pulsed EFP (the v2 generator) ---
    t, Yp, pulses = m.simulate(t_max=T_MAX, dt=DT)
    E2p, InhBp, FSHp, LHp, P4p, GnRHp = Yp.T
    coup_p = MM(GnRHp, P['K_G'])

    # matched-mean flat coupling = cycle-average of the pulsed coupling over the last period
    period = 1.0 / (P['f_min'] + (P['f_max']-P['f_min'])*Hp(EFP['E2'], P['K_freq'], P['n_freq']))
    MMbar = coup_p[t >= t[-1] - period].mean()
    tf, Yf = simulate_flat(MMbar)
    E2f, InhBf, FSHf, LHf, P4f = Yf.T
    coup_f = np.full_like(tf, MMbar)

    th = t / 60.0
    # rows: (label, pulsed series, flat series, EFP level)
    ROWS = [("GnRH -> pituitary\ndrive MM(GnRH;K_G)", coup_p, coup_f, MMbar),
            ("LH (IU/L)",         LHp,   LHf,   EFP['LH']),
            ("FSH (IU/L)",        FSHp,  FSHf,  EFP['FSH']),
            ("E2 (pmol/L)",       E2p,   E2f,   EFP['E2']),
            ("Inhibin B (pg/mL)", InhBp, InhBf, EFP['InhB']),
            ("P4 (nmol/L)",       P4p,   P4f,   EFP['P4'])]
    COLS = ["EFP with GnRH pulses\n(v2 impulsive generator)",
            "EFP with flat GnRH\n(constant, matched mean)"]

    fig, ax = plt.subplots(len(ROWS), 2, figsize=(12, 12), sharex=True)
    for j, col_lab in enumerate(COLS):
        ax[0, j].set_title(col_lab, fontsize=11, fontweight="bold")
    for i, (lab, sp, sf, lvl) in enumerate(ROWS):
        for j, s in enumerate((sp, sf)):
            a = ax[i, j]
            a.plot(th, s, color=f"C{i}", lw=0.9)
            if i > 0:
                a.axhline(lvl, color="k", ls=":", lw=0.9, alpha=0.6, label=f"EFP {lvl:g}")
            else:
                a.axhline(lvl, color="k", ls=":", lw=0.9, alpha=0.6, label=f"mean {lvl:.3f}")
            a.legend(loc="upper right", fontsize=7); a.grid(alpha=0.25)
        ax[i, 0].set_ylabel(lab, fontsize=9)
        lo = min(ax[i, 0].get_ylim()[0], ax[i, 1].get_ylim()[0])
        hi = max(ax[i, 0].get_ylim()[1], ax[i, 1].get_ylim()[1])
        ax[i, 0].set_ylim(lo, hi); ax[i, 1].set_ylim(lo, hi)
    for a in ax[-1]:
        a.set_xlabel("time (hours)")
    fig.suptitle("EFP follicular operating point: GnRH pulses add ripple; flat GnRH gives the "
                 "same equilibrium", fontsize=12, y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out = "fig_efp_pulse_vs_flat.png"
    fig.savefig(out, dpi=130, facecolor="white")
    print(f"saved {out}")
    print(f"  {len(pulses)} pulses / {T_MAX/60:.0f} h, period {period:.0f} min, matched mean coupling {MMbar:.4f}")
    print(f"  pulsed means: LH {LHp.mean():.2f} FSH {FSHp.mean():.2f} E2 {E2p.mean():.1f}")
    print(f"  flat  ends  : LH {LHf[-1]:.2f} FSH {FSHf[-1]:.2f} E2 {E2f[-1]:.1f}")
