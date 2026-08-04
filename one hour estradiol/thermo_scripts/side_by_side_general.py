"""Model beside the ORIGINAL raw multi-day figure (not just the folded average).

Top row  : untouched crop of each CBT-MA panel from mousy temp.png (days 13-18).
Bottom   : the independent-arm model's 24 h waveform tiled across the same 5 days,
           with the digitised mouse trace (raw multi-day, not folded) behind it.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
TITLES = ["Non-pregnant (Cir8)", "Pregnant-CTL (Di36)", "Pregnant-Esr1i (Di7)"]
hrs = fold[LABELS[0]][0]

x = np.load("best_general_arms.npy")
bL, gL, wL, bK, gK, wK, rho, A_T, T0 = x[:9]
ph = x[9:12]


def sig(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


def model_24h(E2, phase):
    c = np.cos(2 * np.pi * (hrs - phase) / 24.0)
    return T0 + A_T * (sig(bK + gK * c + wK * E2) - sig(bL + gL * c + wL * E2))


E2s = [0.0, 1.0, rho]
M = np.vstack([model_24h(E2s[i], ph[i]) for i in range(3)])

SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert("RGB"))
Rc, Gc, Bc = [img[..., k].astype(float) for k in range(3)]
red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False

XSPAN = [(116, 662), (739, 1285), (1362, 1907)]
YANCH = [[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)],
         [(312.5, 38.0), (367.5, 37.0), (423.0, 36.0)],
         [(320.5, 38.0), (372.0, 37.0), (423.0, 36.0)]]
ROW_TOP, ROW_BOT = 313, 423

days5 = np.linspace(13, 18, 5 * len(hrs), endpoint=False)

fig, axes = plt.subplots(2, 3, figsize=(16, 6.6))
for i in range(3):
    x0, x1 = XSPAN[i]
    rows, temps = zip(*YANCH[i])
    m, c = np.polyfit(rows, temps, 1)
    ytop, ybot = m * ROW_TOP + c, m * ROW_BOT + c

    crop = img[ROW_TOP:ROW_BOT + 1, x0:x1 + 1]
    axes[0, i].imshow(crop, extent=[13, 18, ybot, ytop], aspect="auto")
    axes[0, i].set_title(f"{TITLES[i]} -- original figure (raw, not folded)", fontsize=10)
    axes[0, i].set_ylim(ybot, ytop)

    xs, ys = [], []
    for xc in range(x0, x1 + 1):
        col = np.nonzero(red[:, xc])[0]
        if col.size:
            xs.append(xc); ys.append(col.mean())
    tmouse = 13.0 + 5.0 * (np.array(xs, float) - x0) / (x1 - x0)
    Tmouse = m * np.array(ys) + c
    axes[1, i].plot(tmouse, Tmouse, lw=1.0, color="0.55", label="mouse (raw, digitised)")
    axes[1, i].plot(days5, np.tile(M[i], 5), lw=1.8, color="C2", label="independent-arm model")
    axes[1, i].set_ylim(ybot, ytop); axes[1, i].set_xlim(13, 18)
    axes[1, i].grid(alpha=0.3)
    axes[1, i].set_title(f"model (amp {M[i].max()-M[i].min():.2f} C) vs mouse "
                         f"(amp {Tmouse.max()-Tmouse.min():.2f} C raw)", fontsize=9)
    axes[1, i].set_xlabel("day")
    for ax in (axes[0, i], axes[1, i]):
        ax.set_xticks(range(13, 19))
axes[0, 0].set_ylabel("CBT-MA (C)")
axes[1, 0].set_ylabel("CBT-MA (C)")
axes[1, 0].legend(fontsize=8, loc="lower left")
plt.tight_layout()
plt.savefig("fig_side_by_side_general.png", dpi=145)
print("wrote fig_side_by_side_general.png")
for i in range(3):
    print(f"  {LABELS[i]:9s} model daily amp {M[i].max()-M[i].min():.2f} C")
