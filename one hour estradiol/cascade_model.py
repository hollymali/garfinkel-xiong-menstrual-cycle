"""
Single source of the GnRH -> LH -> E2 cascade integrator (the "one hour estradiol"
model).  Every figure/test script imports `run` or `run_delay` from here, so the
model lives in ONE place: edit the loop below and the change applies everywhere.

    dGnRH/dt = -kG*GnRH   (+ a_n at each pulse)
    dLH/dt   =  cL*GnRH - kL*LH
    dE2/dt   =  cE*LH   - kE*E2
    a_n      = A*[circ] / (1 + (E2d/K_amp)^n_amp)                # NEGATIVE amplitude fb
    period   = Pmin + (Pmax-Pmin)/(1 + (E2d/K_freq)^n_freq)      # POSITIVE frequency fb

E2d = E2(t - tau) drives BOTH feedbacks.  tau=0 recovers the instantaneous
(no-delay) case; tau>0 is the physiological E2->GnRH delay that makes the orbit
ragged / high-period.

Integrator: explicit Euler, fixed step dt (min), sequential (Gauss-Seidel) update
within a step (E uses the freshly-updated L, L the freshly-updated G).  Decay rates
are grounded in hormone half-lives.  Companion period-# metrics live in e2metrics.py.
"""
import numpy as np

# anchored natural E2 mean (pg/mL): the shared feedback threshold K_amp = K_freq = K
K = 22.0
DAY = 1440.0

# grounded decay rates (1/min) from hormone half-lives (t1/2 = ln2 / k)
kG = np.log(2)/4.0     # GnRH t1/2 = 4 min
kL = np.log(2)/20.0    # LH   t1/2 = 20 min
kE = np.log(2)/30.0    # E2   t1/2 = 30 min


def circadian(t, t_peak=14*60):
    """Daily envelope (0..1), peaking at t_peak minutes into the day."""
    return 0.5*(1.0 + np.cos(2*np.pi*(np.mod(t, DAY) - t_peak)/DAY))


def run(n_amp, n_freq, K_amp=K, K_freq=K, tau=0.0, t_end=16000.0, dt=0.25,
        e0=15.0, cL=0.5, cE=0.1145, A=8.0, Pmin=40.0, Pmax=110.0, circ=False):
    """Integrate the cascade.  E2(t - tau) drives both E2 feedbacks.

    Returns T, G, L, E, pulse_t, pulse_amp  (arrays).  `pulse_t`/`pulse_amp` are
    the times and amplitudes of the fired GnRH pulses.  If you only need the E2
    time series, use run_delay(), which returns just (T, E).
    """
    steps = int(t_end/dt)
    ndelay = int(round(tau/dt))
    T = np.arange(steps)*dt
    G_arr = np.empty(steps); L_arr = np.empty(steps); E_arr = np.empty(steps)
    G, L, E = 0.0, 5.0, e0
    phi = 0.0
    G_arr[0], L_arr[0], E_arr[0] = G, L, E
    pulse_t, pulse_amp = [], []
    for i in range(steps-1):
        j = i - ndelay
        ed = max(E_arr[j] if j >= 0 else e0, 0.0)              # E2(t - tau)
        period = Pmin + (Pmax - Pmin)/(1.0 + (ed/K_freq)**n_freq)
        phin = phi + dt/period
        pulsed = np.floor(phin) > np.floor(phi)
        phi = phin
        if pulsed:
            amp = A*(0.35 + 0.65*circadian(T[i])) if circ else A
            a = amp/(1.0 + (ed/K_amp)**n_amp)
            pulse_t.append(T[i]); pulse_amp.append(a)
        else:
            a = 0.0
        G = G + dt*(-kG*G) + a
        L = L + dt*(cL*G - kL*L)
        E = E + dt*(cE*L - kE*E)
        G_arr[i+1], L_arr[i+1], E_arr[i+1] = G, L, E
    return T, G_arr, L_arr, E_arr, np.array(pulse_t), np.array(pulse_amp)


def run_delay(n_amp, n_freq, tau, K_amp=K, K_freq=K, t_end=16000.0, dt=0.25,
              e0=15.0, cL=0.5, cE=0.1145, A=8.0, Pmin=40.0, Pmax=110.0, circ=False):
    """Convenience wrapper returning just (T, E) from run()."""
    T, G, L, E, pt, pa = run(n_amp, n_freq, K_amp=K_amp, K_freq=K_freq, tau=tau,
                             t_end=t_end, dt=dt, e0=e0, cL=cL, cE=cE, A=A,
                             Pmin=Pmin, Pmax=Pmax, circ=circ)
    return T, E
