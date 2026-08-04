"""Model vs the ORIGINAL (unfolded) mouse trace, not the 5-day average it was fitted to."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

s = np.load('series.npy', allow_pickle=True).item()
fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
p = np.load('best_lambda_notch2.npy')
bK, gK, wK, pK, bL, gL, wL, pL, phi, lam, rho, A_T = p
E2 = np.array([0., 1., 0.])
midd = np.array([13, 14, 15, 16, 17], float) + 0.5
p4b = np.maximum(0., 1. - (midd - 13.) / 6.5).mean()
P4v = np.array([0., p4b, p4b])
OBS = np.vstack([fold[l][1] for l in L])


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def model(i, t_hours):
    R = rho if i == 2 else 1.0
    c = np.cos(2 * np.pi * (t_hours - phi) / 24.)
    return A_T * (sig(bK + gK * c + wK * E2[i] + pK * P4v[i])
                  - lam * R * sig(bL + gL * c + wL * E2[i] + pL * P4v[i]))


# common offset, exactly as fitted (fold-mean matched)
Wf = np.stack([model(i, hrs) for i in range(3)])
off = OBS.mean() - Wf.mean()

fig, ax = plt.subplots(3, 2, figsize=(15, 9),
                       gridspec_kw=dict(width_ratios=[2.4, 1]))
for i, l in enumerate(L):
    t, y = s[l]
    # ---- left: the ORIGINAL continuous 5-day trace, model repeated each day
    a = ax[i, 0]
    a.plot(t, y, 'k-', lw=1.0, label='mouse (raw, as recorded)')
    tm = np.linspace(t.min(), t.max(), 4000)
    a.plot(tm, model(i, (tm % 1.0) * 24.) + off, 'r-', lw=1.6, alpha=.85,
           label='model (same shape every day)')
    a.set_ylabel(f'{l}\nCBT (C)')
    a.set_xlim(t.min(), t.max())
    for d in range(13, 19):
        a.axvline(d, color='0.85', lw=.8, zorder=0)
    if i == 0:
        a.legend(fontsize=8, ncol=2)
        a.set_title('ORIGINAL time series (gestational day 13-18)')
    if i == 2:
        a.set_xlabel('gestational day')
    # ---- right: every single day folded, vs fold-mean vs model
    b = ax[i, 1]
    for d in range(13, 18):
        m = (t >= d) & (t < d + 1)
        if m.sum() < 20:
            continue
        b.plot((t[m] - d) * 24., y[m], '-', color='0.65', lw=.9)
    b.plot(hrs, OBS[i], 'k-', lw=2.2, label='5-day mean (fitted)')
    b.plot(hrs, model(i, hrs) + off, 'r--', lw=2, label='model')
    b.set_xticks([0, 6, 12, 18, 24]); b.set_xlim(0, 24)
    if i == 0:
        b.legend(fontsize=8)
        b.set_title('each day (grey) vs mean vs model')
    if i == 2:
        b.set_xlabel('hour of day')
fig.suptitle('Model (lambda=0.62, rho=0.31) against the ORIGINAL mouse recording', y=.995)
fig.tight_layout()
fig.savefig('fig_vs_raw.png', dpi=125, bbox_inches='tight')

print('day-to-day spread vs model error (C):', flush=True)
print(f'{"animal":10s}{"model RMS":>11}{"mouse day-to-day SD":>22}{"verdict":>28}', flush=True)
for i, l in enumerate(L):
    t, y = s[l]
    days = []
    for d in range(13, 18):
        m = (t >= d) & (t < d + 1)
        if m.sum() < 20:
            continue
        days.append(np.interp(hrs, (t[m] - d) * 24., y[m]))
    D = np.vstack(days)
    sd = float(np.mean(D.std(0)))
    rms = float(np.sqrt(np.mean((model(i, hrs) + off - OBS[i]) ** 2)))
    v = 'within animal variability' if rms <= sd else f'{rms/sd:.1f}x the animal spread'
    print(f'{l:10s}{rms:11.3f}{sd:22.3f}{v:>28}', flush=True)
print('\nwrote fig_vs_raw.png', flush=True)
