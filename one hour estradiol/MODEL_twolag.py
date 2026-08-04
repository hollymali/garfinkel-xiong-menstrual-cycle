"""
Variant of MODEL.run with three independently-placed lags + an optional KNDy slow
recovery variable (dynorphin), so we can ask WHERE the loop lag has to live and
whether an extra *state* can substitute for a pure delay.

    dGnRH/dt = -kG*GnRH                       (+ a_n at each pulse)
    dLH/dt   =  cL*GnRH - kL*LH
    dE2/dt   =  cE*LH(t - tau_prod) - kE*E2          <- production (steroidogenesis) lag
    a_n      = A / (1 + (E2(t-tau_fb)/K)^n_amp)      <- feedback read lag
    period   = Pmin + (Pmax-Pmin)/(1 + (E2(t-tau_fb)/K)^n_freq)

  KNDy option (dyn=True): the generator is no longer a pure metronome.  A slow
  dynorphin-like inhibitor D is dumped at every pulse and decays; it slows the
  phase velocity, i.e. a refractory / recovery state.

    dD/dt      = -kD*D          (D += 1 at each pulse)
    dphi/dt    = (1/period) * 1/(1 + (D/K_D)^m_D)
"""
import numpy as np

kG = np.log(2) / 4.0
kL = np.log(2) / 20.0
kE = np.log(2) / 30.0
cL = 0.5
cE = 0.1145
K = 22.0
A = 8.0


DAY = 1440.0


def circadian(t, t_peak=14 * 60.0):
    return 0.5 * (1 + np.cos(2 * np.pi * (t - t_peak) / DAY))


def run(n_amp, n_freq, tau_fb=0.0, tau_prod=0.0, t_end=12960, dt=0.25, e0=15,
        period_range=(70, 80), dyn=False, kD=np.log(2) / 40.0, K_D=1.0, m_D=2.0,
        aD=1.0, circ=False, circ_depth=0.65):
    steps = int(t_end / dt)
    n_fb = int(round(tau_fb / dt))
    n_pr = int(round(tau_prod / dt))
    time = np.arange(steps + 1) * dt

    GnRH = np.empty(steps + 1)
    LH = np.empty(steps + 1)
    E2 = np.empty(steps + 1)
    GnRH[0], LH[0], E2[0] = 0.0, 5.0, e0
    phi = 0.0
    D = 0.0
    period_min, period_max = period_range

    for i in range(steps):
        j = i - n_fb
        e2d = E2[j] if j >= 0 else e0
        k = i - n_pr
        lhd = LH[k] if k >= 0 else LH[0]

        period = period_min + (period_max - period_min) / (1 + (e2d / K) ** n_freq)
        rate = 1.0 / period
        if dyn:
            rate /= (1.0 + (D / K_D) ** m_D)
        new_phi = phi + dt * rate
        fire = np.floor(new_phi) > np.floor(phi)
        phi = new_phi

        if fire:
            A_eff = A * ((1 - circ_depth) + circ_depth * circadian(time[i])) if circ else A
            amp = A_eff / (1 + (e2d / K) ** n_amp)
        else:
            amp = 0.0

        GnRH[i + 1] = GnRH[i] + dt * (-kG * GnRH[i]) + amp
        LH[i + 1] = LH[i] + dt * (cL * GnRH[i] - kL * LH[i])
        E2[i + 1] = E2[i] + dt * (cE * lhd - kE * E2[i])
        if dyn:
            D = D + dt * (-kD * D) + (aD if fire else 0.0)

    return time, GnRH, LH, E2
