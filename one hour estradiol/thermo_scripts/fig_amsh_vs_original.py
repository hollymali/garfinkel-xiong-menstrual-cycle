"""Time series: the one-node opposing-SCN-arm model vs the original two-arm fit.

Three panels (one per animal), digitised mouse mean +/- SD across days behind,
three fitted curves on top:

  ORIGINAL (dual, RMS)  - the two-opposed-sigmoid model behind fig_ts_all3.png
  AMSH (RMS)            - one node + opposing SCN arms, fitted on waveform RMS
  AMSH (amplitudes)     - same model, fitted on the three amplitudes

The point of the figure is the amplitude-vs-shape tension: AMSH-amplitudes lands
the peak-to-trough exactly and loses waveform shape; AMSH-RMS does the reverse.
"""
import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
TITLES = ["Non-pregnant (Cir8)", "Pregnant-CTL (Di36)", "Pregnant-Esr1i (Di7)"]
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
SD = np.vstack([fold[l][2] for l in LABELS])
AMP = OBS.max(1) - OBS.min(1)

src = io.open("freephase.py", encoding="utf-8").read().split("ALL = np.array")[0]
ns = {}
exec(src, ns)
plant, sig, clocks = ns["plant"], ns["sig"], ns["clocks"]


def amsh(x):
    R, Q, x0, wE = x[:4]; A_T, T0, rho = x[4:7]; ph = x[7:10]
    E2 = np.array([0., 1., rho])[:, None]
    return plant(sig(R * (1. - Q * E2) * clocks(ph) + x0 + wE * E2), A_T, 1e-6, T0)


def dual(x):
    b, g, w, _ = x[:4]; A_T, T0, rho = x[4:7]; ph = x[7:10]
    E2 = np.array([0., 1., rho])[:, None]; c = clocks(ph)
    return plant(sig(b + g * c + w * E2) - sig(b - g * c + w * E2), A_T, 1e-6, T0)


def rms(M):
    r = M - OBS
    return float(np.sqrt(np.mean((np.sqrt(np.mean(r ** 2, 1)) / AMP) ** 2)))


def amperr(M):
    a = M.max(1) - M.min(1)
    return float(np.sqrt(np.mean(((a - AMP) / AMP) ** 2)))


BA = [(0, 40), (-2, 2), (-8, 8), (-8, 8), (0.1, 15), (30, 42), (0., 1.)] + [(0., 24.)] * 3
BD = [(-8, 8), (0, 40), (-25, 25), (-1, 1), (0.1, 15), (30, 42), (0., 1.)] + [(0., 24.)] * 3


def fit(fn, B, obj, nstart=120, seed=31):
    rng = np.random.default_rng(seed)
    lo = np.array([q[0] for q in B]); hi = np.array([q[1] for q in B])
    best = None
    for _ in range(nstart):
        x0 = lo + rng.random(len(B)) * (hi - lo)
        r = minimize(lambda z: obj(fn(z)), x0, bounds=B, method="L-BFGS-B",
                     options=dict(maxiter=400, ftol=1e-14))
        if best is None or r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best[1]


x_dual = fit(dual, BD, rms)
x_arms = fit(amsh, BA, rms)
x_amp = fit(amsh, BA, lambda M: rms(M) + 3.0 * amperr(M))   # amplitude-only leaves
#          T0 and phase unconstrained (rms ~2.4) AND is degenerate in Q: both
#          Q=+0.66 (arms cancel) and Q=-2.0 (arms reinforce into saturation)
#          reproduce the amplitudes exactly. The combined objective pins the curve.
M_dual, M_arms, M_amp = dual(x_dual), amsh(x_arms), amsh(x_amp)

# categorical hues, fixed order, CVD-validated (see validate_palette.js)
C_DUAL, C_ARMS, C_AMP = "#0173B2", "#DE8F05", "#029E73"
INK, MUTED = "#22201d", "#6b6660"

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
for i, ax in enumerate(axes):
    ax.fill_between(hrs, OBS[i] - SD[i], OBS[i] + SD[i], color="0.80", alpha=0.55,
                    lw=0, label="mouse ± SD across days" if i == 0 else None,
                    zorder=1)
    ax.plot(hrs, OBS[i], "o", ms=3.2, color=INK, zorder=5,
            label="mouse (digitised)" if i == 0 else None)
    ax.plot(hrs, M_dual[i], lw=2, color=C_DUAL, zorder=3,
            label="ORIGINAL: two output arms (RMS fit)" if i == 0 else None)
    ax.plot(hrs, M_arms[i], lw=2, color=C_ARMS, zorder=3,
            label="opposing SCN arms (RMS fit)" if i == 0 else None)
    ax.plot(hrs, M_amp[i], lw=2.4, color=C_AMP, zorder=4,
            label="opposing SCN arms (RMS + amplitude)" if i == 0 else None)

    a = [AMP[i], np.ptp(M_dual[i]), np.ptp(M_arms[i]), np.ptp(M_amp[i])]
    ax.text(0.03, 0.03,
            f"amplitude (C)\nmouse      {a[0]:.2f}\noriginal   {a[1]:.2f}\n"
            f"SCN arms   {a[2]:.2f}\namp fit    {a[3]:.2f}",
            transform=ax.transAxes, fontsize=8, va="bottom", family="monospace",
            color=MUTED, bbox=dict(boxstyle="round,pad=0.4", fc="white",
                                   ec="0.85", alpha=0.9))
    ax.set_title(TITLES[i], fontsize=10, color=INK)
    ax.set_xlabel("hour of day", fontsize=9, color=MUTED)
    ax.set_xticks(range(0, 25, 6))
    ax.set_xlim(0, 24)
    ax.grid(alpha=0.18, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("0.8")
    ax.tick_params(colors=MUTED, labelsize=8)

axes[0].set_ylabel("core body temperature (C)", fontsize=9, color=INK)
fig.legend(loc="lower center", ncol=5, frameon=False, fontsize=8.5,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Opposing-SCN-arm model vs the original two-arm fit  "
             "(P4 dropped; phases free per animal)", fontsize=11, color=INK)
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
fig.savefig("fig_amsh_vs_original.png", dpi=145, facecolor="white",
            bbox_inches="tight")

print(f"{'model':32s}{'amps':>22}{'ratios':>22}{'rms':>9}")
print(f"{'mouse (observed)':32s}{str(np.round(AMP,2)):>22}"
      f"{str(np.round(AMP/AMP[0],2)):>22}{'-':>9}")
for nm, M in (("ORIGINAL two output arms", M_dual),
              ("opposing SCN arms (RMS)", M_arms),
              ("opposing SCN arms (RMS+amp)", M_amp)):
    a = M.max(1) - M.min(1)
    print(f"{nm:32s}{str(np.round(a,2)):>22}{str(np.round(a/a[0],2)):>22}{rms(M):9.4f}")
print(f"\nQ (RMS fit) = {x_arms[1]:.3f}   Q (amplitude fit) = {x_amp[1]:.3f}")
print("saved fig_amsh_vs_original.png")
