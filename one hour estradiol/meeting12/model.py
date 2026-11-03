"""Fixed HPG model. Time is in hours, starting at 07:00."""
import numpy as np


def simulate_e2():
    # Basal production, responsive production, and the two ovarian Hill powers
    basal, production = 6, 210.7
    activation_power, inhibition_power = 5, 5

    # Central feedback and pulse timing
    midpoint, central_power = 30.2, 2.36
    fast_minutes, slow_minutes = 67.5, 76.73
    sleep_multiplier = 1.55
    scn_min, scn_max, scn_peak = 35, 277, 16.8  # Peak at 23:48

    # Half-lives are in minutes. Convert to clearance rates per hour.
    gnrh_decay, lh_decay, e2_decay = 60 * np.log(2) / np.array([4.09, 29, 22.53])

    dt = 1.5 / 3600  # 1.5-second Euler step
    time = np.arange(round((96 + 1 + 72) / dt) + 1) * dt
    e2 = np.empty(len(time))
    e2[0] = 19
    gnrh, lh, phase = 0, 3, 0
    n = midpoint**central_power / (midpoint**central_power + e2[0]**central_power)
    response = 1 / (1 + lh**inhibition_power)

    for i, t in enumerate(time[:-1]):
        # A GnRH pulse occurs whenever phase crosses a whole number.
        drive = n**central_power / (1 + n**central_power)
        interval = slow_minutes - (slow_minutes - fast_minutes) * drive
        clock = (7 + t) % 24
        if clock >= 23 or clock < 7:
            interval *= sleep_multiplier
        next_phase = phase + 60 * dt / interval
        pulse = int(next_phase) - int(phase)
        scn = (scn_max + scn_min) / 2 + (scn_max - scn_min) / 2 * np.cos(
            2 * np.pi * (t % 24 - scn_peak) / 24
        )

        # Compute all five rates from the current state before updating anything.
        dn = midpoint**central_power / (midpoint**central_power + e2[i]**central_power) - n
        dg = -gnrh_decay * gnrh
        dl = gnrh - lh_decay * lh
        dr = 1 / (1 + lh**inhibition_power) - response
        activation = lh**activation_power / (1 + lh**activation_power)
        de = basal + production * activation * response - e2_decay * e2[i]

        n += dt * dn
        gnrh += dt * dg + drive * scn * pulse
        lh += dt * dl
        response += dt * dr
        e2[i + 1] = e2[i] + dt * de
        phase = next_phase

    return time, e2
