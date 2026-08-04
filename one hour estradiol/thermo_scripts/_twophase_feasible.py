"""Two-phase clock: does letting the KEEP and LOSS arms see the SCN at DIFFERENT phases
break the amplitude-vs-modality wall?

WHY THIS FORM (and not a naive morning/evening sum):
two sinusoids of the SAME 24 h frequency added into the SAME arm collapse to a single
sinusoid (g1*cos(wt-p1) + g2*cos(wt-p2) = R*cos(wt-psi)), so an M+E sum inside one arm buys
exactly nothing.  The dual-oscillator SCN only helps if the two arms receive DIFFERENT phase
mixtures -- i.e. separate M/E cell groups projecting differentially to the heat-keep and
heat-loss limbs.  That is one extra parameter (the phase gap) and it is what this tests.

    keep = sig(bK + gK*cos(2pi(t-phiK)/24) + wK*E2 + pK*P4)
    loss = sig(bL + gL*cos(2pi(t-phiL)/24) + wL*E2 + pL*P4)
    CBT  = T0 + A_T*(keep - loss)

E2 = 0/1/0 and P4 = 0/line/line for non-preg / pregnant / Esr1i (Holly's encoding: Esr1i is
the ERalpha-blocked PREGNANT animal, so it carries P4 but no E2 action).  All coefficients
kept >= 0 per Holly's all-positive call; phases are free.

WHY IT COULD WORK where one shared phase cannot: with phiK != phiL the output stops being a
function of a single variable c(t), which is what forced (i) mirror symmetry about noon,
(ii) exactly-equal peak heights, and (iii) the amplitude-poor "fold" route to bimodality.
Here the loss arm firing mid-plateau CARVES A NOTCH in the keep plateau -- two humps at full
amplitude, notch depth set independently by the loss arm.

Reports feasibility at progressively stricter targets so we learn WHERE it breaks, not just
whether it breaks.
"""
import numpy as np

fold = np.load('fold.npy', allow_pickle=True).item()
L = ['non-preg', 'pregnant', 'Esr1i']
hrs = fold[L[0]][0]
OBS = np.vstack([fold[l][1] for l in L])
AMP = OBS.max(1) - OBS.min(1)
R_OBS = AMP / AMP[0]

# fine grid for honest peak geometry (0.2 h); generic phases make exact ties measure-zero
tg = np.arange(0, 24, 0.2)
NT = tg.size
E2 = np.array([0., 1., 0.])
PREG = np.array([0., 1., 1.])
mid = np.array([13, 14, 15, 16, 17], float) + 0.5
P4d = np.maximum(0., 1. - (mid - 13.) / 6.5)
p4bar = P4d.mean()
P4v = np.array([0., p4bar, p4bar])          # day-averaged P4 per animal

PROM = 0.10
OBS_GAP = 7.0        # non-preg 8.2 -> 15.2
OBS_MIDNP = 11.7     # non-preg pair midpoint
OBS_MIDPR = 18.8     # pregnant pair midpoint  (noon-lock forbids this outright)


def sig(z):
    return 1 / (1 + np.exp(-np.clip(z, -60, 60)))


def waveform(b, g, phi, w, p, e2, p4):
    """(M,) params + scalar condition -> (M,NT)"""
    c = np.cos(2 * np.pi * (tg[None, :] - phi[:, None]) / 24.)
    return sig(b[:, None] + g[:, None] * c + w[:, None] * e2 + p[:, None] * p4)


def peak_stats(W):
    """circular peak geometry with a REAL prominence test.
    returns n_signif_peaks, gap_h, midpoint_h, amp   (gap/mid = 0 unless exactly 2 peaks)"""
    amp = W.max(1) - W.min(1)
    lo = np.roll(W, 1, 1); hi = np.roll(W, -1, 1)
    ismax = (W > lo) & (W > hi)
    ismin = (W < lo) & (W < hi)
    nmax = ismax.sum(1)
    gap = np.zeros(len(W)); midp = np.zeros(len(W)); nsig = nmax.copy()

    two = np.flatnonzero((nmax == 2) & (ismin.sum(1) == 2) & (amp > 1e-9))
    if two.size:
        rm, cm = np.nonzero(ismax[two])
        pk = cm.reshape(-1, 2)
        rn, cn = np.nonzero(ismin[two])
        vl = cn.reshape(-1, 2)
        Wt = W[two]
        pv = np.take_along_axis(Wt, pk, 1)
        vv = np.take_along_axis(Wt, vl, 1)
        # the SHALLOWER dip is the one that decides 2 peaks vs 1 broad peak
        shallow = vv.max(1)
        depth = pv.min(1) - shallow
        real2 = depth > PROM * amp[two]
        t1 = tg[pk[:, 0]]; t2 = tg[pk[:, 1]]
        g_lin = t2 - t1
        g_circ = np.minimum(g_lin, 24 - g_lin)
        # midpoint of the pair on the side that is actually the short arc
        m_lin = (t1 + t2) / 2.
        m_alt = (m_lin + 12.) % 24.
        midp[two] = np.where(g_lin <= 12., m_lin, m_alt)
        gap[two] = g_circ
        nsig[two] = np.where(real2, 2, 1)
    return nsig, gap, midp, amp


rng = np.random.default_rng(90210)
N, CH = 3_000_000, 20_000
cnt = dict(np2=0, np2gap=0, es1=0, both=0, both_pr2=0, all_topo=0, all_topo_amp=0)
best = None
best_err = np.inf
keep_seeds = []
for s in range(0, N, CH):
    bK = rng.uniform(-15, 8, CH);  gK = rng.uniform(0, 30, CH)
    bL = rng.uniform(-15, 8, CH);  gL = rng.uniform(0, 30, CH)
    phiK = rng.uniform(0, 24, CH); phiL = rng.uniform(0, 24, CH)
    wK = rng.uniform(0, 30, CH);   wL = rng.uniform(0, 30, CH)
    pK = rng.uniform(0, 15, CH);   pL = rng.uniform(0, 15, CH)

    st = []
    for i in range(3):
        H = (waveform(bK, gK, phiK, wK, pK, E2[i], P4v[i])
             - waveform(bL, gL, phiL, wL, pL, E2[i], P4v[i]))
        st.append(peak_stats(H) + (H,))

    n_np, gap_np, mid_np, a_np, H_np = st[0]
    n_pr, gap_pr, mid_pr, a_pr, H_pr = st[1]
    n_es, gap_es, mid_es, a_es, H_es = st[2]

    m_np2 = n_np == 2
    m_np2g = m_np2 & (np.abs(gap_np - OBS_GAP) < 2.0)
    m_es1 = n_es == 1
    m_both = m_np2g & m_es1
    m_bpr = m_both & (n_pr == 2)
    m_topo = m_bpr & (np.abs(mid_pr - OBS_MIDPR) < 2.0)
    cnt['np2'] += int(m_np2.sum()); cnt['np2gap'] += int(m_np2g.sum())
    cnt['es1'] += int(m_es1.sum()); cnt['both'] += int(m_both.sum())
    cnt['both_pr2'] += int(m_bpr.sum()); cnt['all_topo'] += int(m_topo.sum())

    ok = np.flatnonzero(m_topo)
    if ok.size:
        an = np.maximum(a_np[ok], 1e-12)
        err = (np.abs(a_pr[ok] / an - R_OBS[1]) + np.abs(a_es[ok] / an - R_OBS[2]))
        good = err < 0.15
        cnt['all_topo_amp'] += int(good.sum())
        j = int(np.argmin(err))
        if err[j] < best_err:
            k = ok[j]
            best_err = float(err[j])
            best = dict(bK=bK[k], gK=gK[k], phiK=phiK[k], wK=wK[k], pK=pK[k],
                        bL=bL[k], gL=gL[k], phiL=phiL[k], wL=wL[k], pL=pL[k],
                        amps=(a_np[k], a_pr[k], a_es[k]),
                        gap=gap_np[k], midnp=mid_np[k], midpr=mid_pr[k])
        if good.any():
            kk = ok[good]
            keep_seeds.append(np.stack([bL[kk], gL[kk], phiL[kk], wL[kk], pL[kk],
                                        bK[kk], gK[kk], phiK[kk], wK[kk], pK[kk]], 1))

print(f'draws {N:,}   (two free phases, all coeffs >= 0)')
print(f'  non-preg 2 significant peaks                       {cnt["np2"]:,}')
print(f'  ... with gap {OBS_GAP}+-2 h                              {cnt["np2gap"]:,}')
print(f'  Esr1i 1 peak                                       {cnt["es1"]:,}')
print(f'  non-preg bimodal(gap) AND Esr1i unimodal           {cnt["both"]:,}')
print(f'  ... and pregnant also bimodal                      {cnt["both_pr2"]:,}')
print(f'  ... and pregnant pair centred on {OBS_MIDPR} h +-2       {cnt["all_topo"]:,}'
      f'   <-- NOON-LOCK FORBIDS THIS ENTIRELY')
print(f'  ... AND amplitude ratios within 0.15               {cnt["all_topo_amp"]:,}')

if best is not None:
    an, ap, ae = best['amps']
    print(f'\n  best amplitude-ratio match among full-topology draws:')
    print(f'    ratios 1.00/{ap/an:.2f}/{ae/an:.2f}   observed 1.00/{R_OBS[1]:.2f}/{R_OBS[2]:.2f}'
          f'   |err| {best_err:.3f}')
    print(f'    non-preg gap {best["gap"]:.1f} h (obs {OBS_GAP}), midpoint {best["midnp"]:.1f} h'
          f' (obs {OBS_MIDNP}); pregnant midpoint {best["midpr"]:.1f} h (obs {OBS_MIDPR})')
    print(f'    phiK={best["phiK"]:.2f} h  phiL={best["phiL"]:.2f} h'
          f'  -> PHASE GAP {abs(best["phiK"]-best["phiL"]):.2f} h')
    print(f'    KEEP b={best["bK"]:+.2f} g={best["gK"]:.2f} w={best["wK"]:.2f} p={best["pK"]:.2f}')
    print(f'    LOSS b={best["bL"]:+.2f} g={best["gL"]:.2f} w={best["wL"]:.2f} p={best["pL"]:.2f}')
    print(f'    non-preg arm amplitude {an:.4f} -> needs A_T = {AMP[0]/an:.1f} C'
          f'  (physical bound ~<=10 C)')
else:
    print('\n  NO draw satisfied the full topology target.')

if keep_seeds:
    S = np.concatenate(keep_seeds)
    np.save('_twophase_seeds.npy', S)
    print(f'\n  saved {len(S):,} seeds -> _twophase_seeds.npy  (ready for Nelder-Mead polish)')
