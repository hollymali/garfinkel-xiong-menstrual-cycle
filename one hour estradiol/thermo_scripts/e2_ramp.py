"""Let E2 RAMP across gestational days 13->18 instead of sitting at a constant 1.

The mechanism (coefficients of the no-delay, noon-locked, all-excitatory model) is
held FIXED from best_noon_excite.npy -- nothing is refit except the ramp slope.

    E2_pregnant(day) = 1 + s*(day - 15.5)        midpoint 15.5 keeps the 5-day mean at 1,
    E2_Esr1i(day)    = rho * E2_pregnant(day)    so the original fold-fit stays consistent
    E2_nonpreg(day)  = 0                          (no gestational rise)

TEST: fit s on the PREGNANT animal's daily means only, then PREDICT the Esr1i
decline out-of-sample. rho is not refit -- it comes from the original fit.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

fold = np.load("fold.npy", allow_pickle=True).item()
LABELS = ["non-preg", "pregnant", "Esr1i"]
TITLES = ["Non-pregnant (Cir8)", "Pregnant-CTL (Di36)", "Pregnant-Esr1i (Di7)"]
hrs = fold[LABELS[0]][0]

x = np.load("best_noon_excite.npy")
bL, gL, wL, bK, gK, wK, rho, A_T, T0 = x
c = np.cos(2 * np.pi * (hrs - 12.0) / 24.0)


def sig(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


def day_wave(E2):
    return T0 + A_T * (sig(bK + gK * c + wK * E2) - sig(bL + gL * c + wL * E2))


DAYS = [13, 14, 15, 16, 17]


def e2_for(animal, day, s):
    ramp = 1.0 + s * (day + 0.5 - 15.5)          # day+0.5 = mid-day
    if animal == 0:
        return 0.0
    if animal == 1:
        return ramp
    return rho * ramp


def series(animal, s):
    """5 days of waveform, each day at that day's E2."""
    return np.concatenate([day_wave(e2_for(animal, d, s)) for d in DAYS])


def daily_means(animal, s):
    return np.array([day_wave(e2_for(animal, d, s)).mean() for d in DAYS])


# ---- observed daily means (from the raw digitised traces) ----
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

raw_t, raw_T, obs_means = [], [], []
for i in range(3):
    x0, x1 = XSPAN[i]
    rows, temps = zip(*YANCH[i])
    m, cc = np.polyfit(rows, temps, 1)
    xs, ys = [], []
    for xc in range(x0, x1 + 1):
        col = np.nonzero(red[:, xc])[0]
        if col.size:
            xs.append(xc); ys.append(col.mean())
    t = 13.0 + 5.0 * (np.array(xs, float) - x0) / (x1 - x0)
    T = m * np.array(ys) + cc
    raw_t.append(t); raw_T.append(T)
    obs_means.append(np.array([T[(t >= d) & (t < d + 1)].mean() for d in DAYS]))
obs_means = np.vstack(obs_means)

# ---- fit s on the PREGNANT animal's daily means only ----
def cost(s):
    pred = daily_means(1, s)
    return float(np.sum(((pred - pred.mean()) - (obs_means[1] - obs_means[1].mean())) ** 2))


res = minimize_scalar(cost, bounds=(-1.0, 1.0), method="bounded")
s = float(res.x)

print(f"fitted ramp slope s = {s:+.4f} per day   "
      f"(E2 goes {1+s*(13.5-15.5):.2f} -> {1+s*(17.5-15.5):.2f} across days 13->18)")
print()
print(f"{'animal':10s} {'observed decline':>18} {'model decline':>15}   {'fitted or predicted':>20}")
for i, l in enumerate(LABELS):
    od = obs_means[i][-1] - obs_means[i][0]
    md = daily_means(i, s)[-1] - daily_means(i, s)[0]
    tag = "FITTED" if i == 1 else ("PREDICTED (out of sample)" if i == 2 else "no E2 ramp -> 0")
    print(f"{l:10s} {od:+18.2f} {md:+15.2f}   {tag:>20}")

# ---- figure ----
days_grid = np.concatenate([13 + k + hrs / 24.0 for k in range(5)])
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharey=False)
for i, ax in enumerate(axes):
    ax.plot(raw_t[i], raw_T[i], lw=1.0, color="0.55", label="mouse (raw, digitised)")
    ax.plot(days_grid, series(i, s), lw=1.7, color="#029E73",
            label="model, E2 ramping across gestation")
    ax.plot(days_grid, np.concatenate([day_wave(e2_for(i, d, 0.0)) for d in DAYS]),
            lw=1.2, color="#0173B2", alpha=0.65, ls="--", label="model, constant E2")
    ax.set_title(TITLES[i], fontsize=10)
    ax.set_xlabel("day"); ax.set_xticks(range(13, 19)); ax.set_xlim(13, 18)
    ax.grid(alpha=0.25)
axes[0].set_ylabel("CBT (C)")
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle(f"E2 ramping across gestation (slope fit on PREGNANT only, s={s:+.3f}/day)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("fig_e2_ramp_side_by_side.png", dpi=145, facecolor="white")
print("\nsaved fig_e2_ramp_side_by_side.png")
