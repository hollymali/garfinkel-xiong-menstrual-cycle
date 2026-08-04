"""
Full extraction of the three CBT-MA traces, calibrated, folded to 24 h.

Calibration anchors are the panel's own horizontal GRIDLINES (known integer
temperatures) plus the frame, and the vertical gridlines (day boundaries).
Each panel has its own y-limits -- they are NOT the same across panels.
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert('RGB')).astype(float)
Rc, Gc, Bc = img[..., 0], img[..., 1], img[..., 2]
gray = img.mean(axis=2)

red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False                       # drop the red underline in the title

LABELS = ['non-preg', 'pregnant', 'Esr1i']
XSPAN = [(116, 662), (739, 1285), (1362, 1907)]

# y calibration: (row, temperature) anchors read off the detected gridlines
#   panel 1 has three interior gridlines (38/37/36)
#   panel 2 has one (37) + frame top/bottom at 38/36
#   panel 3 has two (38/37) + frame bottom at 36
YANCH = [[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)],
         [(312.5, 38.0), (367.5, 37.0), (423.0, 36.0)],
         [(320.5, 38.0), (372.0, 37.0), (423.0, 36.0)]]

def vertical_gridlines(x0, x1):
    """Day boundaries: light-grey vertical lines inside the plot box."""
    band = gray[318:418, x0 - 30:x1 + 30]
    greyfrac = ((band > 170) & (band < 245)).mean(axis=0)
    cand = np.nonzero(greyfrac > 0.55)[0] + x0 - 30
    groups, cur = [], [cand[0]]
    for v in cand[1:]:
        if v - cur[-1] <= 2:
            cur.append(v)
        else:
            groups.append(np.mean(cur)); cur = [v]
    groups.append(np.mean(cur))
    return np.array(groups)

series = {}
for i, (lab, (x0, x1)) in enumerate(zip(LABELS, XSPAN)):
    rows, temps = zip(*YANCH[i])
    m, c = np.polyfit(rows, temps, 1)          # row -> temperature
    vg = vertical_gridlines(x0, x1)
    print(f"\n{lab}: {m:.5f} C/px, gridlines at days -> {np.round(vg,1)}")

    # trace: mean row of red pixels in each column
    xs, ys = [], []
    for x in range(x0, x1 + 1):
        col = np.nonzero(red[:, x])[0]
        if col.size:
            xs.append(x); ys.append(col.mean())
    xs, ys = np.array(xs, float), np.array(ys)
    T = m * ys + c

    # x -> days. Every panel spans days 13-18 across its plot box, so map the
    # frame directly. (Fitting the detected vertical lines is unsafe: panels 2
    # and 3 carry a dashed event marker that is not a day gridline.)
    t_days = 13.0 + 5.0 * (xs - x0) / (x1 - x0)
    interior = vg[(vg > x0 + 10) & (vg < x1 - 10)]
    pred = 13.0 + 5.0 * (interior - x0) / (x1 - x0)
    print(f"   {(x1-x0)/5:.1f} px/day; detected verticals fall at days "
          f"{np.round(pred, 2)}  (day marks should be integers)")
    print(f"   T range {T.min():.2f} - {T.max():.2f} C   mean {T.mean():.2f}")
    series[lab] = (t_days, T)

np.save('series.npy', np.array(series, dtype=object), allow_pickle=True)

# ------------------------------------------------------- fold to 24 h
NB = 48                                        # 30-min bins
fold = {}
fig, axes = plt.subplots(2, 3, figsize=(15, 6))
for j, lab in enumerate(LABELS):
    td, T = series[lab]
    axes[0, j].plot(td, T, lw=1, color='C3')
    axes[0, j].set_title(f"{lab}  (extracted)"); axes[0, j].grid(alpha=0.3)
    axes[0, j].set_xlabel("day")

    ph = (td % 1.0) * 24.0
    idx = np.clip((ph / 24.0 * NB).astype(int), 0, NB - 1)
    mean = np.array([T[idx == k].mean() if (idx == k).any() else np.nan
                     for k in range(NB)])
    sd = np.array([T[idx == k].std() if (idx == k).sum() > 1 else np.nan
                   for k in range(NB)])
    hrs = (np.arange(NB) + 0.5) * 24 / NB
    axes[1, j].plot(hrs, mean, color='C0')
    axes[1, j].fill_between(hrs, mean - sd, mean + sd, alpha=0.25, color='C0')
    axes[1, j].set_title(f"{lab}  24 h fold"); axes[1, j].grid(alpha=0.3)
    axes[1, j].set_xlabel("hour of day"); axes[1, j].set_ylim(35.4, 38.5)
    fold[lab] = (hrs, mean, sd)

    amp = np.nanmax(mean) - np.nanmin(mean)
    mn = np.nanmean(mean)
    tr = np.nanmin(mean)
    flat = np.nanmean(np.abs(mean - mn)) / (amp / 2)   # 0.64 sine, 1.0 square
    print(f"{lab:9s} folded: amp {amp:.2f}  mean {mn:.2f}  trough {tr:.2f}  "
          f"shape {(mn-tr)/amp:.2f}  FLATNESS {flat:.2f}  "
          f"(mean within-bin SD {np.nanmean(sd):.2f})")

axes[0, 0].set_ylabel("CBT-MA (C)"); axes[1, 0].set_ylabel("CBT-MA (C)")
plt.tight_layout(); plt.savefig("fig_digitized.png", dpi=130)
np.save('fold.npy', np.array(fold, dtype=object), allow_pickle=True)
print("\nwrote fig_digitized.png")
