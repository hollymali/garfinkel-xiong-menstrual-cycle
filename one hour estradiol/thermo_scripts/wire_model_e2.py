"""Drive the non-pregnant arm with REAL E2(t) from MODEL.py (the pulsatile
mid-late-follicular cascade), instead of the constant E2=0 used so far.

Everything else is held fixed from the already-fitted independent-arm model
(best_general_arms.npy): pregnant and Esr1i stay at their constant E2=1, rho;
only the non-pregnant arm's E2 input changes, from 0 to a real time series.
This is a SIMULATION with the existing fit, not a refit -- an honest
out-of-sample check of whether ragged E2 explains what constant E2=0 couldn't.

E2 unit bridge (same physiological anchor used earlier in coupled_e2_temp.py):
thermo E2=0 is non-pregnant, E2=1 is pregnant. preg_ratio = pregnant E2 / cycling
mean E2 maps the cascade's pg/mL swing into that scale:
    E2_in(t) = (E2_cascade(t) - mean) / (mean * (preg_ratio - 1))
preg_ratio is an ASSUMPTION (no hard citation on hand); reported at several values.
"""
import sys
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol")
import MODEL

DAY = 1440.0

# ---------------- fitted mechanism (frozen, not refit) ----------------
x = np.load("best_general_arms.npy")
bL, gL, wL, bK, gK, wK, rho, A_T, T0 = x[:9]
ph = x[9:12]           # phases for [non-preg, pregnant, Esr1i]


def sig(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


def cbt(E2, phase, t_min):
    c = np.cos(2 * np.pi * (t_min / 60.0 - phase) / 24.0)
    return T0 + A_T * (sig(bK + gK * c + wK * E2) - sig(bL + gL * c + wL * E2))


# ---------------- real E2(t) from the cascade, band (55,65), tau=60, n=4 ----------------
BURN_D, NDAYS = 4.0, 5.0
t_end = int((BURN_D + NDAYS) * DAY)
t_c, G, L, E_pgml = MODEL.run(4, 4, tau=60, t_end=t_end, dt=0.25,
                              period_range=(55, 65), e0=15)
m = t_c >= BURN_D * DAY
t_rel = t_c[m] - t_c[m][0]          # minutes, 0..NDAYS*DAY
E_pgml = E_pgml[m]
mean_e = E_pgml.mean()
print(f"cascade E2 over {NDAYS:.0f} d: {E_pgml.min():.1f}-{E_pgml.max():.1f} pg/mL, mean {mean_e:.1f}")

# ---------------- raw digitised non-pregnant trace (for comparison) ----------------
SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert("RGB"))
Rc, Gc, Bc = [img[..., k].astype(float) for k in range(3)]
red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False
x0, x1 = 116, 662
rows, temps = zip(*[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)])
mcal, ccal = np.polyfit(rows, temps, 1)
xs, ys = [], []
for xc in range(x0, x1 + 1):
    col = np.nonzero(red[:, xc])[0]
    if col.size:
        xs.append(xc); ys.append(col.mean())
t_mouse_day = 13.0 + 5.0 * (np.array(xs, float) - x0) / (x1 - x0)
T_mouse = mcal * np.array(ys) + ccal
mouse_amp_raw = T_mouse.max() - T_mouse.min()

fig, axes = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
ax0, ax1 = axes
day_grid = 13.0 + t_rel / DAY

# baseline for reference: constant E2=0 (what we had before)
cbt0 = cbt(0.0, ph[0], t_rel)
ax0.plot(t_mouse_day, T_mouse, lw=1.0, color="0.55", label="mouse (raw, digitised)")
ax0.plot(day_grid, cbt0, lw=1.6, color="C0", alpha=0.7, label="constant E2=0 (previous)")
ax0.set_ylabel("CBT (C)")
ax0.set_title(f"non-pregnant: constant E2=0  |  amp {cbt0.max()-cbt0.min():.2f} C "
              f"(mouse raw {mouse_amp_raw:.2f} C)")
ax0.legend(fontsize=8, loc="lower left")
ax0.grid(alpha=0.25)

print(f"\n{'preg_ratio':>10} {'E2_in swing':>14} {'CBT amp (C)':>12} {'gain over E2=0':>15}")
colors = ["C2", "C3", "C4"]
for i, pr in enumerate((2.0, 4.0, 8.0)):
    e2_in = (E_pgml - mean_e) / (mean_e * (pr - 1.0))
    cbt_r = cbt(e2_in, ph[0], t_rel)
    amp = cbt_r.max() - cbt_r.min()
    print(f"{pr:10.1f} {e2_in.max()-e2_in.min():14.4f} {amp:12.2f} {amp/(cbt0.max()-cbt0.min()):15.2f}x")
    if pr == 4.0:
        ax1.plot(day_grid, cbt_r, lw=1.3, color=colors[i], alpha=0.9,
                label=f"MODEL.py E2(t), preg_ratio={pr:.0f}")

ax1.plot(t_mouse_day, T_mouse, lw=1.0, color="0.55", label="mouse (raw, digitised)")
ax1.set_ylabel("CBT (C)")
ax1.set_xlabel("day")
ax1.set_title("non-pregnant: real ragged E2(t) from MODEL.py drives the SAME sigmoid arms")
ax1.legend(fontsize=8, loc="lower left")
ax1.grid(alpha=0.25)
for ax in axes:
    ax.set_xticks(range(13, 19))
    ax.set_xlim(13, 18)

fig.tight_layout()
fig.savefig("fig_model_e2_noncycle.png", dpi=145, facecolor="white")
print("\nsaved fig_model_e2_noncycle.png")
