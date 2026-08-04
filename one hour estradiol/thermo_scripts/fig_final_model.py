"""FINAL temperature model figure.

Top row    : folded 24 h waveform, model vs mouse mean +/- SD across days
Bottom row : same model tiled across days 13-18, against the raw digitised trace

Model = noon-locked SCN, all-excitatory coefficients into both arms, no delay,
no P4, constant E2 per animal.  Coefficients from best_noon_excite.npy.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
TITLES = ["Non-pregnant (Cir8)", "Pregnant-CTL (Di36)", "Pregnant-Esr1i (Di7)"]
hrs = fold[LABELS[0]][0]
OBS = np.vstack([fold[l][1] for l in LABELS])
SD = np.vstack([fold[l][2] for l in LABELS])

x = np.load("best_noon_excite.npy")
bL, gL, wL, bK, gK, wK, rho, A_T, T0 = x
c = np.cos(2 * np.pi * (hrs - 12.0) / 24.0)


def sig(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


E2s = [0.0, 1.0, rho]
M = np.vstack([T0 + A_T * (sig(bK + gK * c + wK * E2s[i]) - sig(bL + gL * c + wL * E2s[i]))
               for i in range(3)])

SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert("RGB"))
Rc, Gc, Bc = [img[..., k].astype(float) for k in range(3)]
red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False
XSPAN = [(116, 662), (739, 1285), (1362, 1907)]
YANCH = [[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)],
         [(312.5, 38.0), (367.5, 37.0), (423.0, 36.0)],
         [(320.5, 38.0), (372.0, 37.0), (423.0, 36.0)]]

INK, MUTED = "#22201d", "#6b6660"
CMOD = "#0173B2"
days5 = np.concatenate([13 + k + hrs / 24.0 for k in range(5)])

fig, axes = plt.subplots(2, 3, figsize=(16, 8.2))
for i in range(3):
    ax = axes[0, i]
    ax.fill_between(hrs, OBS[i] - SD[i], OBS[i] + SD[i], color="0.80", alpha=0.55, lw=0,
                    label="mouse ± SD across days" if i == 0 else None)
    ax.plot(hrs, OBS[i], "o", ms=3.2, color=INK, label="mouse (folded)" if i == 0 else None)
    ax.plot(hrs, M[i], lw=2.4, color=CMOD, label="model" if i == 0 else None)
    ax.axvline(12, color="0.6", ls="--", lw=0.8)
    ax.set_title(TITLES[i], fontsize=10, color=INK)
    ax.set_xlabel("hour of day", fontsize=9, color=MUTED)
    ax.set_xticks(range(0, 25, 6)); ax.set_xlim(0, 24); ax.grid(alpha=0.18, lw=0.6)
    ax.text(0.03, 0.03, f"amp mouse {np.ptp(OBS[i]):.2f}\namp model {np.ptp(M[i]):.2f}",
            transform=ax.transAxes, fontsize=8, va="bottom", family="monospace",
            color=MUTED, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.85"))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)

    x0, x1 = XSPAN[i]
    rows, temps = zip(*YANCH[i])
    m, cc = np.polyfit(rows, temps, 1)
    xs, ys = [], []
    for xc in range(x0, x1 + 1):
        col = np.nonzero(red[:, xc])[0]
        if col.size:
            xs.append(xc); ys.append(col.mean())
    tm = 13.0 + 5.0 * (np.array(xs, float) - x0) / (x1 - x0)
    Tm = m * np.array(ys) + cc
    ax = axes[1, i]
    ax.plot(tm, Tm, lw=1.0, color="0.55", label="mouse (raw)" if i == 0 else None)
    ax.plot(days5, np.tile(M[i], 5), lw=1.6, color=CMOD, label="model" if i == 0 else None)
    ax.set_xlabel("day", fontsize=9, color=MUTED)
    ax.set_xticks(range(13, 19)); ax.set_xlim(13, 18); ax.grid(alpha=0.18, lw=0.6)
    ax.set_title(f"raw: mouse amp {np.ptp(Tm):.2f} C vs model {np.ptp(M[i]):.2f} C",
                 fontsize=9, color=MUTED)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)

axes[0, 0].set_ylabel("CBT (C) — folded 24 h", fontsize=9, color=INK)
axes[1, 0].set_ylabel("CBT (C) — raw, days 13–18", fontsize=9, color=INK)
axes[0, 0].legend(fontsize=8, loc="upper left")
axes[1, 0].legend(fontsize=8, loc="lower left")
fig.suptitle("Final model: noon-locked SCN, excitation-only into two saturating arms, "
             "constant E2 per animal", fontsize=12, color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("fig_FINAL_model.png", dpi=145, facecolor="white")
print("saved fig_FINAL_model.png")
print(f"\nbL={bL:.2f} gL={gL:.2f} wL={wL:.2f} | bK={bK:.2f} gK={gK:.2f} wK={wK:.2f}")
print(f"rho={rho:.3f}  A_T={A_T:.2f}  T0={T0:.2f}")
for i, l in enumerate(LABELS):
    print(f"  {l:10s} amp model {np.ptp(M[i]):.2f}  mouse folded {np.ptp(OBS[i]):.2f}")
