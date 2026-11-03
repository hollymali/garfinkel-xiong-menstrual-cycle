"""KNDy -> GnRH -> LH -> E2 pulses with optional circadian amplitude."""

import numpy as np


# Hormone decay and cascade constants
kG = np.log(2) / 4
kL = np.log(2) / 20
E2_HALF_LIFE_MINUTES = 30.0  # Rounded from 27.45 +/- 5.65 min (PMID 9855336)
kE = np.log(2) / E2_HALF_LIFE_MINUTES
cL = 0.5
cE = 0.1145

# Feedback and pulse constants
K = 20.1
A = 6.8
N_HALF = 0.28
DAY_MINUTES = 1440.0
CIRCADIAN_PEAK_MINUTES = 14.0 * 60.0
CIRCADIAN_MODULATION = 0.7


def circadian(time_minutes):
    """Return a simple 24-hour cosine (-1..1), peaking at 14:00."""
    time_of_day = np.mod(time_minutes, DAY_MINUTES)
    return np.cos(
        2.0
        * np.pi
        * (time_of_day - CIRCADIAN_PEAK_MINUTES)
        / DAY_MINUTES
    )


def run(n_amp=4, n_freq=2, tau=0, t_end=9 * 1440, dt=0.1, e0=15,
        LH0=5, period_range=(63, 114), circ=True, return_n=False):
    """Simulate the model and return time, GnRH, LH, and E2 (plus N if asked)."""
    if tau != 0:
        raise ValueError("N dynamics replace the explicit time delay, so tau must be 0")

    steps = int(t_end / dt)
    time = np.arange(steps + 1) * dt
    N = np.empty(steps + 1)
    GnRH = np.empty(steps + 1)
    LH = np.empty(steps + 1)
    E2 = np.empty(steps + 1)

    N[0] = 1 / (1 + (e0 / K) ** 5)
    GnRH[0], LH[0], E2[0] = 0, LH0, e0
    phase = 0
    pmin, pmax = period_range

    for i in range(steps):
        amp_drive = N[i] ** n_amp / (N_HALF ** n_amp + N[i] ** n_amp)
        freq_drive = N[i] ** n_freq / (N_HALF ** n_freq + N[i] ** n_freq)
        period = pmin + (pmax - pmin) * freq_drive

        new_phase = phase + dt / period
        # N is continuous; GnRH secretion is still represented as discrete pulses.
        fire = np.floor(new_phase) > np.floor(phase)
        amplitude_multiplier = 1.0
        if circ:
            amplitude_multiplier += (
                CIRCADIAN_MODULATION * (1.0 + circadian(time[i])) / 2.0
            )
        pulse = A * amp_drive * amplitude_multiplier if fire else 0
        phase = new_phase

        # dN/dt = E2-suppressed production - N decay
        N[i + 1] = N[i] + dt * (0.25 / (1 + (E2[i] / K) ** 5) - 0.25 * N[i])
        GnRH[i + 1] = GnRH[i] + dt * (-kG * GnRH[i]) + pulse
        LH[i + 1] = LH[i] + dt * (cL * GnRH[i] - kL * LH[i])
        E2[i + 1] = E2[i] + dt * (cE * LH[i] - kE * E2[i])

    if return_n:
        return time, N, GnRH, LH, E2
    return time, GnRH, LH, E2
