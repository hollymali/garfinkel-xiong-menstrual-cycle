"""Reconstruct the fast GnRH-driven LH PULSES underneath the v5 slow cycle.

v5 (sigmoid_model_v5.py) replaces the fast GnRH pulse train with its EXACT cycle-average coupling,
so its LH curve is smooth (the running mean of the pulses).  This script recovers the spikes: it
takes v5's slow E2 & P4 at three cycle phases, FREEZES them, and integrates the underlying impulsive
GnRH pulse generator + instantaneous LH for 24 h at each phase.  The average of these pulses is
exactly the smooth v5 LH -- here we show the pulses themselves and how E2 (pulse frequency + surge
amplitude) and P4 (luteal clock brake) reshape them across the cycle.

Run from the repo root:  python reconstruct_pulses.py  ->  fig_v5_pulses.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sigmoid_model_v5 as m

P = m.P
Hp, Hm = m.Hp, m.Hm
MM = lambda x, K: x / (K + x)

# ---- v5 slow cycle, aligned to menses (E2 nadir = cycle day 0) ----
t, Y = m.simulate(days=200.0, dt=0.5)
td = t / m.DAY
pk = m.surge_peaks(t, Y)
a0 = pk[-2]
imin = a0 + int(np.argmin(Y[slice(a0, pk[-1] + 1), 0]))
period = td[pk[-1]] - td[pk[-2]]
ov_day = float(np.argmax(Y[slice(imin, imin + int(period * m.DAY / 0.5)), 0])) * 0.5 / m.DAY

def state_at(day_off):
    return Y[imin + int(day_off * m.DAY / 0.5)]

phases = [("Early follicular (day 3)", 3.0),
          (f"Preovulatory SURGE (day {ov_day:.1f})", ov_day),
          (f"Mid-luteal (day {ov_day + 6:.0f})", ov_day + 6.0)]


def run_pulses(E2, P4, LH0, hours=24.0, dt=0.02):
    """Freeze the slow state; integrate impulsive GnRH kicks + instantaneous LH on the fast scale."""
    A_max = P['A'] / Hm(P['E2_ref'], P['K_amp'], P['n_amp'])
    a = A_max * (Hm(E2, P['K_amp'], P['n_amp']) + P['A_sg'] * Hp(E2, P['K_sens'], P['n_sens']))  # kick
    f = (P['f_min'] + (P['f_max'] - P['f_min']) * Hp(E2, P['K_freq'], P['n_freq'])) \
        * Hm(P4, P['K_P4f'], P['n_P4f'])                                                          # freq
    gL = max(1 - P['B_LH'] * Hp(E2, P['K_neg'], P['n_neg'])
               + P['A_LH'] * Hp(E2, P['K_sens'], P['n_sens']), 0.0)
    n = int(hours * 60 / dt)
    T = np.linspace(0, hours * 60, n)
    L = np.empty(n)
    phi, g, lh = 0.0, a, LH0
    for i in range(n):
        phin = phi + f * dt
        if np.floor(phin) > np.floor(phi):
            g += a
        phi = phin
        g += -P['k_GnRH'] * g * dt
        lh += (P['V_L'] * MM(g, P['K_G']) * gL - P['k_LH'] * lh) * dt
        L[i] = lh
    return T / 60.0, L, a, 1.0 / f


if __name__ == "__main__":
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for k, (lab, doff) in enumerate(phases):
        E2, InhB, FSH, LH, P4, F, Lut = state_at(doff)
        Th, L, a, interval = run_pulses(E2, P4, LH)
        ax[k].plot(Th, L, color="tab:red", lw=1.1)
        ax[k].set_ylabel("LH (IU/L)")
        ax[k].set_title(f"{lab}:  E2={E2:.0f} pmol/L, P4={P4:.1f} nmol/L  ->  "
                        f"1 pulse / {interval:.0f} min, kick={a:.2f}", fontsize=10)
        ax[k].grid(alpha=0.25)
        print(f"{lab}: E2={E2:.0f} P4={P4:.1f} interval={interval:.0f}min kick={a:.2f} meanLH={L.mean():.1f}")
    ax[2].set_xlabel("time (hours)  [24 h window with the slow state frozen]")
    fig.suptitle("Fast GnRH-driven LH pulses underneath the v5 slow cycle (reconstructed at 3 phases)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig("fig_v5_pulses.png", dpi=130, facecolor="white")
    print("saved fig_v5_pulses.png")
