"""Final endocrine-temperature model and a small fixed-step solver.

The model equations and fitted parameter values are documented in README.md.
``run()`` returns the final 120-hour series used by ``plot_model.py``; this
module intentionally contains no plotting or file-export code.
"""

from __future__ import annotations

import numpy as np

# Simulation grid (hours): one shared 6-second step for the whole model.
STEP_H = 1.0 / 600.0
WARMUP_H = 96.0
OUTPUT_H = 120.0

# Endocrine model (hours, pg/mL)
GNRH_HALF_LIFE_H = 1.0 / 15.0
LH_HALF_LIFE_H = 1.0 / 3.0
E2_HALF_LIFE_H = 1.0 / 2.0
K_GNRH = np.log(2.0) / GNRH_HALF_LIFE_H
K_LH = np.log(2.0) / LH_HALF_LIFE_H
K_E2 = np.log(2.0) / E2_HALF_LIFE_H
LH_GAIN = 30.0
E2_GAIN = 6.87
N_RATE = 15.0
E2_FEEDBACK_HALF = 20.1
PULSE_AMPLITUDE = 6.8
N_PULSE_HALF = 0.28
PULSE_PERIOD_LOW_H = 1.05
PULSE_PERIOD_HIGH_H = 1.90
ENDOCRINE_CLOCK_PEAK_H = 14.0
CLOCK_GAIN = 0.7

# Final temperature model (hours, degrees Celsius)
THERMAL_CLOCK_PEAK_H = 14.0
CONSERVATION_CLOCK_WEIGHT = 2.27
CONSERVATION_N_WEIGHT = 0.60
PROGESTERONE_WEIGHT = 0.0274
LOSS_CLOCK_WEIGHT = 0.78
E2_LOSS_WEIGHT = 0.0000869
BASAL_HEAT_RATE = 49.44
CONSERVATION_GAIN = 9.76
HORMONAL_LOSS_GAIN = 15.29
SELF_LOSS_MAX = 125.0
SELF_LOSS_HALF_C = 39.08
INITIAL_TEMPERATURE_C = 37.0

# Hormone conditions
CYCLING_PROGESTERONE = 0.2
PREGNANCY_E2 = 215.0
PREGNANCY_PROGESTERONE = 25.0


def endocrine_clock(time_h: float) -> float:
    """Endocrine clock multiplier in [1, 1.7], peaking at 14:00."""
    unit_clock = 0.5 * (
        1.0
        + np.cos(
            2.0
            * np.pi
            * ((time_h % 24.0) - ENDOCRINE_CLOCK_PEAK_H)
            / 24.0
        )
    )
    return float(1.0 + CLOCK_GAIN * unit_clock)


def hill4(x: np.ndarray) -> np.ndarray:
    """Shared fourth-order Hill response: x^4 / (1 + x^4)."""
    return x**4 / (1.0 + x**4)


def simulate_endocrine(extra_output_h: float = 0.0) -> dict[str, np.ndarray]:
    """Integrate N -> GnRH -> LH -> E2 with the explicit Euler loop."""
    total_h = WARMUP_H + OUTPUT_H + extra_output_h
    steps = round(total_h / STEP_H)
    time_h = np.arange(steps + 1) * STEP_H
    n_state = np.empty(steps + 1)
    gnrh = np.empty(steps + 1)
    lh = np.empty(steps + 1)
    e2 = np.empty(steps + 1)

    n_state[0] = E2_FEEDBACK_HALF**5 / (
        E2_FEEDBACK_HALF**5 + 15.0**5
    )
    gnrh[0], lh[0], e2[0] = 0.0, 5.0, 15.0
    phase = 0.0

    for i in range(steps):
        n_value = n_state[i]
        amplitude_drive = n_value**4 / (N_PULSE_HALF**4 + n_value**4)
        frequency_drive = n_value**2 / (N_PULSE_HALF**2 + n_value**2)
        period_h = PULSE_PERIOD_LOW_H + (
            PULSE_PERIOD_HIGH_H - PULSE_PERIOD_LOW_H
        ) * frequency_drive

        new_phase = phase + STEP_H / period_h
        fires = np.floor(new_phase) > np.floor(phase)
        pulse = (
            PULSE_AMPLITUDE * amplitude_drive * endocrine_clock(time_h[i])
            if fires
            else 0.0
        )
        phase = new_phase

        n_state[i + 1] = n_value + STEP_H * (
            N_RATE
            * E2_FEEDBACK_HALF**5
            / (E2_FEEDBACK_HALF**5 + e2[i] ** 5)
            - N_RATE * n_value
        )
        gnrh[i + 1] = (
            gnrh[i] - STEP_H * K_GNRH * gnrh[i] + pulse
        )
        lh[i + 1] = lh[i] + STEP_H * (
            LH_GAIN * gnrh[i] - K_LH * lh[i]
        )
        e2[i + 1] = e2[i] + STEP_H * (
            E2_GAIN * lh[i] - K_E2 * e2[i]
        )

    return {
        "time_h": time_h,
        "N": n_state,
        "GnRH": gnrh,
        "LH": lh,
        "E2": e2,
    }


def temperature_rate(
    temperature_c: float, conservation: float, loss: float
) -> float:
    """Evaluate the final temperature ODE in degrees Celsius per hour."""
    return float(
        BASAL_HEAT_RATE
        + CONSERVATION_GAIN * conservation
        - HORMONAL_LOSS_GAIN * loss
        - self_loss(temperature_c)
    )


def integrate_temperature(
    conservation: np.ndarray, loss: np.ndarray
) -> np.ndarray:
    """Advance temperature with the same explicit Euler style as endocrine."""
    temperature = np.empty(len(conservation))
    temperature[0] = INITIAL_TEMPERATURE_C

    for i in range(len(temperature) - 1):
        temperature[i + 1] = temperature[i] + STEP_H * temperature_rate(
            temperature[i], conservation[i], loss[i]
        )

    return temperature

def thermal_clock(time_h: np.ndarray) -> np.ndarray:
    """Thermal pathway clock in [0, 1], peaking at the fitted time."""
    return 0.5 * (
        1.0
        + np.cos(
            2.0
            * np.pi
            * ((time_h % 24.0) - THERMAL_CLOCK_PEAK_H)
            / 24.0
        )
    )


def self_loss(temperature_c: float | np.ndarray) -> float | np.ndarray:
    """Complete temperature-dependent heat-loss rate in degrees C per hour."""
    return (
        SELF_LOSS_MAX * temperature_c**5
        / (SELF_LOSS_HALF_C**5 + temperature_c**5)
    )

def simulate_condition(
    time_h: np.ndarray,
    n_state: np.ndarray,
    e2: np.ndarray,
    progesterone: float,
) -> dict[str, np.ndarray]:
    """Calculate the thermal signals and temperature for one condition."""
    clock = thermal_clock(time_h)
    conservation = hill4(
        CONSERVATION_CLOCK_WEIGHT * clock
        + CONSERVATION_N_WEIGHT * n_state
        + PROGESTERONE_WEIGHT * progesterone
    )
    hormonal_loss = hill4(
        LOSS_CLOCK_WEIGHT * clock + E2_LOSS_WEIGHT * e2
    )
    temperature = integrate_temperature(conservation, hormonal_loss)
    return {
        "hour": time_h,
        "N": n_state,
        "E2_pg_ml": e2,
        "progesterone_ng_ml": np.full(len(time_h), progesterone),
        "clock": clock,
        "conservation": conservation,
        "hormonal_loss": hormonal_loss,
        "self_regulated_loss": np.asarray(self_loss(temperature)),
        "temperature_c": temperature,
    }


def keep_final_120_hours(condition: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Discard the 96-hour warm-up and reset the retained time to hour zero."""
    output_points = round(OUTPUT_H / STEP_H) + 1
    retained = {
        variable_name: values[-output_points:]
        for variable_name, values in condition.items()
    }
    retained["hour"] = retained["hour"] - retained["hour"][0]
    return retained


def run(
    output_start_clock_h: float = 0.0,
) -> dict[str, dict[str, np.ndarray]]:
    """Return 120 hours beginning at the requested clock hour.

    The default remains midnight.  A nonzero clock hour extends the simulation
    before taking its final 120 hours, so clock-anchored windows never wrap or
    reuse data from a different modeled day.
    """
    start_clock_h = output_start_clock_h % 24.0
    endocrine = simulate_endocrine(extra_output_h=start_clock_h)
    output_points = round(OUTPUT_H / STEP_H) + 1
    low_e2 = float(np.min(endocrine["E2"][-output_points:]))
    low_n = E2_FEEDBACK_HALF**5 / (E2_FEEDBACK_HALF**5 + low_e2**5)
    pregnancy_n = E2_FEEDBACK_HALF**5 / (
        E2_FEEDBACK_HALF**5 + PREGNANCY_E2**5
    )
    count = len(endocrine["time_h"])

    cycling = simulate_condition(
        endocrine["time_h"],
        endocrine["N"],
        endocrine["E2"],
        CYCLING_PROGESTERONE,
    )
    cycling["GnRH"] = endocrine["GnRH"]
    cycling["LH"] = endocrine["LH"]

    constant_low = simulate_condition(
        endocrine["time_h"],
        np.full(count, low_n),
        np.full(count, low_e2),
        CYCLING_PROGESTERONE,
    )
    pregnancy = simulate_condition(
        endocrine["time_h"],
        np.full(count, pregnancy_n),
        np.full(count, PREGNANCY_E2),
        PREGNANCY_PROGESTERONE,
    )

    return {
        "cycling": keep_final_120_hours(cycling),
        "constant_low": keep_final_120_hours(constant_low),
        "pregnancy": keep_final_120_hours(pregnancy),
    }
