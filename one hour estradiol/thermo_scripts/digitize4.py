"""
Robustness check on the fold: remove sub-circadian drift before folding.

The pregnant animal drifts down ~0.4 C across days 13-18, and there is an
unexplained dashed event marker at ~day 14.5 in panels 2 and 3. Folding across
a trend suppresses the apparent circadian amplitude, so redo it after
subtracting a centred 24 h moving average (which removes everything slower
than the circadian band and nothing at 24 h).
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = r"C:\Users\holly\source\garfinkel-xiong-menstrual-cycle\one hour estradiol\mousy temp.png"
img = np.asarray(Image.open(SRC).convert('RGB')).astype(float)
Rc, Gc, Bc = img[..., 0], img[..., 1], img[..., 2]
red = (Rc > 110) & (Rc - Gc > 45) & (Rc - Bc > 45)
red[:250, :] = False

LABELS = ['non-preg', 'pregnant', 'Esr1i']
XSPAN = [(116, 662), (739, 1285), (1362, 1907)]
YANCH = [[(326.5, 38.0), (365.5, 37.0), (403.5, 36.0)],
         [(312.5, 38.0), (367.5, 37.0), (423.0, 36.0)],
         [(320.5, 38.0), (372.0, 37.0), (423.0, 36.0)]]

NB = 48
GRID = np.arange(13.0, 18.0, 1.0 / (24 * 4))       # 15-min common grid

def extract(i):
    x0, x1 = XSPAN[i]
    rows, temps = zip(*YANCH[i])
    m, c = np.polyfit(rows, temps, 1)
    xs, ys = [], []
    for x in range(x0, x1 + 1):
        col = np.nonzero(red[:, x])[0]
        if col.size:
            xs.append(x); ys.append(col.mean())
    xs = np.array(xs, float)
    T = m * np.array(ys) + c
    days = 13.0 + 5.0 * (xs - x0) / (x1 - x0)
    return np.interp(GRID, days, T)

def fold(T, hrs_of):
    idx = np.clip((hrs_of / 24.0 * NB).astype(int), 0, NB - 1)
    mean = np.array([T[idx == k].mean() for k in range(NB)])
    sd = np.array([T[idx == k].std() for k in range(NB)])
    return mean, sd

hrs_of = (GRID % 1.0) * 24.0
win = 24 * 4                                        # 24 h in samples

out_raw, out_det = {}, {}
fig, axes = plt.subplots(2, 3, figsize=(15, 6), sharex='row')
print(f"{'cond':10s}{'raw amp':>9s}{'detr amp':>10s}{'ratio raw':>11s}"
      f"{'ratio detr':>12s}{'peak hr':>9s}{'flatness':>10s}")
for i, lab in enumerate(LABELS):
    T = extract(i)
    base = np.convolve(T, np.ones(win) / win, mode='same')
    edge = win // 2
    Td = T.copy()
    Td[edge:-edge] = (T - base + T.mean())[edge:-edge]
    Td = Td[edge:-edge]; hd = hrs_of[edge:-edge]

    m_raw, s_raw = fold(T, hrs_of)
    m_det, s_det = fold(Td, hd)
    out_raw[lab] = (np.arange(NB) * 24 / NB + 0.25, m_raw, s_raw)
    out_det[lab] = (np.arange(NB) * 24 / NB + 0.25, m_det, s_det)

    a_raw = m_raw.max() - m_raw.min()
    a_det = m_det.max() - m_det.min()
    flat = np.mean(np.abs(m_det - m_det.mean())) / (a_det / 2)
    if i == 0:
        A0r, A0d = a_raw, a_det
    print(f"{lab:10s}{a_raw:9.2f}{a_det:10.2f}{a_raw/A0r:11.2f}{a_det/A0d:12.2f}"
          f"{(np.arange(NB)*24/NB)[m_det.argmax()]:9.1f}{flat:10.2f}")

    axes[0, i].plot(GRID, T, lw=0.9, color='C3'); axes[0, i].plot(GRID, base, lw=1.5, color='k')
    axes[0, i].set_title(f"{lab}: trace + 24h baseline"); axes[0, i].grid(alpha=0.3)
    axes[1, i].plot(out_raw[lab][0], m_raw, label='raw fold')
    axes[1, i].plot(out_det[lab][0], m_det, label='detrended fold')
    axes[1, i].fill_between(out_det[lab][0], m_det-s_det, m_det+s_det, alpha=0.2)
    axes[1, i].set_title(f"{lab}: 24 h fold"); axes[1, i].grid(alpha=0.3)
    axes[1, i].set_xlabel("hour of day"); axes[1, i].legend(fontsize=7)

print("\n(sinusoid flatness = 0.64, square wave = 1.00, spiky < 0.5)")
plt.tight_layout(); plt.savefig("fig_detrended.png", dpi=130)
np.save('fold_detr.npy', np.array(out_det, dtype=object), allow_pickle=True)
print("wrote fig_detrended.png, fold_detr.npy")
