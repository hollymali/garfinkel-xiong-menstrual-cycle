"""
E2 doesn't settle when its INPUT doesn't settle.
================================================

Illustrates the readout principle behind the leaky-integrator framework:
E2 is a filter, so E2 is non-uniform exactly when the GnRH pulse train is
non-uniform. E2 never self-oscillates here -- all structure is inherited.

Minimal leaky integrator (one line of biology):
    dE2/dt = A * pulse_input(t) - kE2 * E2
with kE2 giving tau ~ 50 min, so ~1h pulses partially relax between beats and
stay visible (as in the estradiol panel of the body-temp/estradiol figure).

Three inputs, same filter:
  (A) regular pulses            -> E2 settles to a uniform ripple  (the metronome)
  (B) noisy pulses              -> E2 perpetually ragged           (source a)
  (C) noisy + slow AM envelope  -> ragged pulses under a slow hump (a + b: like panel B)
"""

import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(3)

DAY = 1440.0          # minutes in 24 h
kE2 = 0.02            # 1/min  -> tau = 50 min (E2 clearance)
SIG = 6.0            # pulse width (min)


def pulse_train(t, period=60.0, amp=1.0, jitter=0.0, amp_cv=0.0,
                slow_am=0.0, am_phase=0.0):
    """Sum of Gaussian pulses. Optional interval jitter, amplitude noise, and a
    slow (circadian) amplitude envelope."""
    # build pulse times across the window with lognormal interval jitter
    times, amps = [], []
    tt = 30.0
    while tt < t[-1] + period:
        times.append(tt)
        # per-pulse amplitude: lognormal noise * slow circadian envelope
        a = amp
        if amp_cv > 0:
            a *= rng.lognormal(mean=0.0, sigma=amp_cv)
        if slow_am > 0:
            a *= 1.0 + slow_am * np.sin(2*np.pi*tt/DAY - am_phase)
        amps.append(a)
        step = period * (rng.lognormal(0.0, jitter) if jitter > 0 else 1.0)
        tt += step
    times = np.array(times); amps = np.array(amps)
    # forcing on the time grid
    f = np.zeros_like(t)
    for ti, ai in zip(times, amps):
        f += ai * np.exp(-0.5 * ((t - ti) / SIG) ** 2)
    return f


def leaky_e2(t, forcing, A=1.2, e0=15.0):
    """Integrate dE2/dt = A*forcing - kE2*E2 (fixed-step Euler on a fine grid)."""
    dt = t[1] - t[0]
    E2 = np.empty_like(t); E2[0] = e0
    for i in range(len(t) - 1):
        E2[i+1] = E2[i] + dt * (A * forcing[i] - kE2 * E2[i])
    return E2


if __name__ == "__main__":
    t = np.arange(0.0, 2*DAY, 0.5)   # 48 h, show the 2nd day (past transient)
    show = t >= DAY

    cases = {
        "A  regular pulses":                dict(jitter=0.0, amp_cv=0.0, slow_am=0.0),
        "B  noisy pulses":                  dict(jitter=0.22, amp_cv=0.30, slow_am=0.0),
        "C  noisy + slow (circadian) AM":   dict(jitter=0.22, amp_cv=0.30, slow_am=0.5,
                                                 am_phase=2*np.pi*8/24),  # peak ~midday
    }

    fig, ax = plt.subplots(len(cases), 1, figsize=(11, 8), sharex=True)
    for k, (label, kw) in enumerate(cases.items()):
        f = pulse_train(t, period=60.0, amp=1.0, **kw)
        E2 = leaky_e2(t, f)
        tail = E2[show]
        print(f"{label:32s}: E2 range {tail.min():5.1f}..{tail.max():5.1f}  "
              f"mean {tail.mean():5.1f}")
        ax[k].plot((t[show]-DAY)/60, E2[show], lw=1.0, color=f"C{k}")
        ax[k].set_ylabel("E2"); ax[k].set_title(label, loc="left", fontsize=10)
        ax[k].grid(alpha=0.3)
    ax[-1].set_xlabel("time of day (h)")
    fig.suptitle("Same leaky filter, different GnRH input: E2 settles only when its input does")
    fig.tight_layout()
    fig.savefig("fig_feed_nonuniform.png", dpi=130)
    print("wrote fig_feed_nonuniform.png")
