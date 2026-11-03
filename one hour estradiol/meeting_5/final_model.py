"""Simple in-arm model with visibly deep non-pregnant and pregnant troughs.

    Conservation = H5(1.1*Clock + 0.8*E2 + 0.1*P4)
    Loss         = H5(0.6*Clock + 0.6*E2)
    CBT          = 35.6 + 4*(Conservation - Loss)
    Clock centers = [12, 7.5, 12] h

The E2 output from MODEL.py is normalized to mean 0.706 and used for the
non-pregnant row. Pregnant E2 remains constant at 1 and Esr1i E2 remains 0.

Restoring the shared loss-arm clock coefficient to 0.6 lets loss catch up near
the clock maximum and carve a clear trough between the two temperature peaks.
The shared temperature gain keeps the outer peaks tall. The condition-specific
clock centers remain unchanged and the model retains 2/2/1 modality.
"""

import numpy as np

from MODEL import run


CONDITIONS = ("non-pregnant", "pregnant", "Esr1i")
E2 = np.array([0.706, 1.0, 0.0])
P4 = np.array([0.0, 1.0, 1.0])
CLOCK_CENTERS = np.array([12.0, 7.5, 12.0])


def hill5(value):
    """Increasing Hill response with exponent 5 and unit half-saturation."""
    return value**5 / (1 + value**5)


def calculate_model(start, stop, step):
    """Generate time points and return E2, Clock, both arms, and absolute CBT."""
    hours = np.arange(start, stop, step)

    minutes, _, _, e2_trace = run(
        4,
        4,
        t_end=(3 * 24 + stop - start) * 60,
    )
    use = minutes >= minutes[-1] - (stop - start) * 60
    e2_trace = np.interp(
        (hours - start) * 60,
        minutes[use] - minutes[use][0],
        e2_trace[use],
    )
    e2_input = np.tile(E2[:, None], (1, hours.size))
    e2_input[0] = E2[0] * e2_trace / e2_trace.mean()

    clock = (
        1
        + np.cos(
            2 * np.pi * (hours[None, :] - CLOCK_CENTERS[:, None]) / 24
        )
    ) / 2

    conservation = hill5(
        1.1 * clock + 0.8 * e2_input + 0.1 * P4[:, None]
    )
    loss = hill5(0.6 * clock + 0.6 * e2_input)
    cbt = 35.6 + 4 * (conservation - loss)

    return hours, e2_input, clock, conservation, loss, cbt
