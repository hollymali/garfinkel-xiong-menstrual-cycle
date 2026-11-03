"""
Mid-late-follicular GnRH -> LH -> E2 cascade + BOTH E2 feedbacks (NO delay here).

    dGnRH/dt = -kG*GnRH (+ a_n at each pulse)
    dLH/dt   =  cL*GnRH - kL*LH
    dE2/dt   =  cE*LH   - kE*E2
    a_n        = A / (1 + (E2/K_amp)^n_amp)                  # NEGATIVE amplitude fb
    period(E2) = Pmin + (Pmax-Pmin)/(1 + (E2/K_freq)^n_freq) # POSITIVE frequency fb

Decay rates are grounded in real hormone half-lives (t1/2 = ln2/k, minutes):
    GnRH t1/2 = 4 min,  LH t1/2 = 20 min,  E2 t1/2 = 30 min.
In mid-late follicular phase LH (pulsatile, t1/2 ~20 min) -- not FSH (tonic,
t1/2 ~3.5 h) -- is the gonadotropin that transmits hourly pulses to E2.
cL, cE are units-only gains (E2/K is invariant under cE rescaling once K is
re-anchored), tuned so E2 sits at ~22 pg/mL (panel-B mid-follicular level).

This is the NO-DELAY reference: across the whole physiological Hill range
(n <= 5) every orbit is period-1/2 -- raggedness requires the E2->GnRH delay
(see make_heatmap_delayed.py).
"""

import numpy as np
import matplotlib.pyplot as plt
from MODEL import run, K   # no-delay case = run(..., tau=0); integrator lives in MODEL.py


if __name__ == "__main__":
    # natural E2 with the feedbacks flat (n=0 -> E2-independent drive); K is anchored to it inside MODEL
    Tn, Gn, Ln, En = run(0.0, 0.0)
    print(f"natural E2 (feedbacks flat) = {En[Tn > Tn[-1]-3000].mean():.1f} pg/mL  (MODEL K = {K})")

    from orbit_period import e2_orbit_period, period_label
    n_amp = 4
    n_freqs = [0, 2, 4, 8]
    fig, ax = plt.subplots(len(n_freqs), 1, figsize=(11, 9), sharex=True)
    for k, nf in enumerate(n_freqs):
        T, G, L, E = run(n_amp, nf)   # tau defaults to 0 -> no delay
        lab = period_label(e2_orbit_period(T, E))
        m = T > T[-1] - 3000
        ax[k].plot(T[m]/60, E[m], lw=0.9, color="C4")
        ax[k].set_ylabel(f"n_freq={nf}\nE2"); ax[k].grid(alpha=0.3)
        ax[k].set_title(f"n_amp={n_amp} (fixed), n_freq={nf}:   orbit = {lab}",
                        loc="left", fontsize=9)
        print(f"n_freq={nf}: orbit {lab}")
    ax[-1].set_xlabel("time (h)")
    fig.suptitle(f"Cascade (GnRH->LH->E2): negative amplitude (n_amp={n_amp}) + positive frequency feedback")
    fig.tight_layout()
    fig.savefig("fig_pi_cascade_combined.png", dpi=130)
    print("wrote fig_pi_cascade_combined.png")
