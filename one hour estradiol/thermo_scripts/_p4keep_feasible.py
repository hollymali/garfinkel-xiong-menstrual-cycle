"""Is (non-preg bimodal, Esr1i UNIMODAL) reachable at all with P4 in the KEEP arm only?

Structural test, not a fit.  Under noon-lock the whole daily waveform is a function
of c alone, c = cos(2pi(t-12)/24), so the peak count follows analytically:

    f(c)  = sig(bK + pK*P4 + wK*E2 + gK*c) - sig(bL + wL*E2 + gL*c)
    f'(c) = gK*sig'(uK) - gL*sig'(uL)

    #daily peaks = 2*(interior maxima of f) + [f'(1)>0] + [f'(-1)<0]

(c sweeps -1 -> +1 -> -1 each day, so every interior max is visited twice; noon is
c=+1 and midnight c=-1.)  UNIMODAL therefore means f is monotone on [-1,1]:
increasing -> the single peak sits at NOON, decreasing -> at MIDNIGHT.  Observed
Esr1i peak is 13.2 h, so it has to be the increasing branch.

Scan the parameters that decide non-preg vs Esr1i (E2=0 in both, so the E2 weights
drop out entirely), all positive per Holly's spec, then add the pregnant animal.
"""
import numpy as np

rng = np.random.default_rng(20260803)
N = 4_000_000
CH = 50_000
C = np.linspace(-1, 1, 241)[None, :]
OBS = dict(np=2.10, pr=0.72, es=1.16)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def dsig(z):
    s = sig(z)
    return s * (1 - s)


def classify(b_keep, gK, b_loss, gL):
    """peaks-per-day and amplitude of f(c)=sig(b_keep+gK c)-sig(b_loss+gL c)."""
    npk = np.empty(len(b_keep), np.int8)
    amp = np.empty(len(b_keep))
    noon = np.empty(len(b_keep), bool)
    for s in range(0, len(b_keep), CH):
        e = slice(s, s + CH)
        uK = b_keep[e][:, None] + gK[e][:, None] * C
        uL = b_loss[e][:, None] + gL[e][:, None] * C
        f = sig(uK) - sig(uL)
        d = gK[e][:, None] * dsig(uK) - gL[e][:, None] * dsig(uL)
        pos = d > 0
        drops = np.count_nonzero(pos[:, :-1] & ~pos[:, 1:], axis=1)   # + -> - crossings
        up_noon = pos[:, -1]
        dn_mid = ~pos[:, 0]
        npk[e] = 2 * drops + up_noon + dn_mid
        amp[e] = f.max(1) - f.min(1)
        noon[e] = up_noon & ~dn_mid
    return npk, amp, noon


bK = rng.uniform(-30, 10, N)
gK = rng.uniform(0, 60, N)
bL = rng.uniform(-30, 10, N)
gL = rng.uniform(0, 60, N)
pK = rng.uniform(0, 25, N)

n_np, a_np, _ = classify(bK, gK, bL, gL)                # non-preg: E2=0, P4=0
n_es, a_es, noon_es = classify(bK + pK, gK, bL, gL)     # Esr1i:    E2=0, P4=1

ok_np = n_np == 2
ok_es = (n_es == 1) & noon_es
both = ok_np & ok_es
print(f'draws                                 {N:,}')
print(f'non-preg bimodal                      {ok_np.sum():,}')
print(f'Esr1i unimodal, peak at noon          {ok_es.sum():,}')
print(f'  (unimodal but peak at MIDNIGHT)     {(n_es == 1).sum() - ok_es.sum():,}')
print(f'BOTH                                  {both.sum():,}')

if not both.sum():
    raise SystemExit('\n=> P4-in-keep-arm-only cannot make Esr1i unimodal. IMPOSSIBLE as specified.')

ratio = np.where(both, a_es / np.maximum(a_np, 1e-12), np.nan)
good = both & (np.abs(ratio - OBS['es'] / OBS['np']) < 0.10)
print(f'BOTH + Esr1i/non-preg amp 0.55+-0.10  {good.sum():,}')

sel = np.flatnonzero(good if good.sum() else both)
print(f'\n  {"bK":>8}{"gK":>7}{"bL":>8}{"gL":>7}{"pK":>7}{"amp_np":>8}{"amp_es":>8}{"ratio":>7}')
for i in sel[:8]:
    print(f'  {bK[i]:8.2f}{gK[i]:7.2f}{bL[i]:8.2f}{gL[i]:7.2f}{pK[i]:7.2f}'
          f'{a_np[i]:8.3f}{a_es[i]:8.3f}{a_es[i]/a_np[i]:7.2f}')

# ---- now demand the pregnant animal be bimodal too (adds wK,wL>0) ----
M = min(len(sel), 20_000)
src = rng.choice(sel, M, replace=False)
K = 100
hit = 0
ex = None
for _ in range(K):
    wK = rng.uniform(0, 60, M)
    wL = rng.uniform(0, 60, M)
    n_pr, a_pr, _ = classify(bK[src] + pK[src] + wK, gK[src], bL[src] + wL, gL[src])
    m = (n_pr == 2) & (np.abs(a_pr / np.maximum(a_np[src], 1e-12)
                             - OBS['pr'] / OBS['np']) < 0.10)
    hit += int(m.sum())
    if ex is None and m.any():
        j = int(np.flatnonzero(m)[0])
        ex = (src[j], wK[j], wL[j], a_pr[j])
print(f'\n  + pregnant bimodal & amp ratio 0.34+-0.10:  {hit:,} of {M*K:,} (wK,wL) draws')
if ex:
    i, wKx, wLx, apr = ex
    print(f'    example: bK={bK[i]:.2f} gK={gK[i]:.2f} bL={bL[i]:.2f} gL={gL[i]:.2f} '
          f'pK={pK[i]:.2f} wK={wKx:.2f} wL={wLx:.2f}')
    print(f'    amps (arm units)  non-preg {a_np[i]:.3f}  pregnant {apr:.3f}  Esr1i {a_es[i]:.3f}')
