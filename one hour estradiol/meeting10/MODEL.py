"""Meeting 10 deterministic endocrine and core-temperature model.

The state variables and equation forms are inherited from the Meeting 8
endocrine subsystem:

    E2 --| N -> pulse timing and amplitude -> GnRH -> LH -> E2

The separate direct pituitary E2-to-LH inhibition term has been removed.  In
particular, LH production is simply proportional to the current GnRH state:

    dLH/dt = G_LH * GnRH - k_LH * LH.

Literature-scale hormone half-lives and the 51--75 minute waking pulse-period
bounds are fixed. During the 23:00--07:00 Licinio sleep window, the current
first-edition hypothesis multiplies the instantaneous pulse period by 1.57.
The phenomenological calibration values are ordinary module-level constants.

The Meeting 8 heat-balance model is included directly in this module, with one
Meeting 9 hypothesis: E2 and a nonnegative SCN-associated BAT drive enter a
shared saturating thermogenic response,

    SCN_BAT = (1 - cos(2*pi*(ZT - 7)/24)) / 2
    X = E2 + k_SCN * SCN_BAT
    H(X) = X**5 / (30**5 + X**5).

The old independent additive SCN temperature term is not used. Temperature
does not feed back into the endocrine subsystem.
"""

import numpy as np


# Simulation grid.  Six-second steps preserve the Meeting 8 integration scale.
STEP_H = 1.0 / 600.0
ZT_PERIOD_H = 24.0
WARMUP_H = 96.0
OUTPUT_H = 72.0

# Literature-scale kinetic anchors retained from Meeting 8.
GNRH_HALF_LIFE_H = 1.0 / 15.0  # 4 min
LH_HALF_LIFE_H = 1.0 / 3.0  # 20 min
E2_HALF_LIFE_H = 22.0 / 60.0  # 22 min; lower literature-defensible bound
K_GNRH = np.log(2.0) / GNRH_HALF_LIFE_H
K_LH = np.log(2.0) / LH_HALF_LIFE_H
K_E2 = np.log(2.0) / E2_HALF_LIFE_H

# Pulse-period endpoints are held fixed to the literature-scale interval used
# in Meeting 8. The N turnover rate and both N response midpoints are fixed at 1.
PULSE_PERIOD_MIN_MIN = 51.0
PULSE_PERIOD_MAX_MIN = 75.0
HILL_COEFFICIENT = 5.0

# Phenomenological Meeting 9 calibration constants.
N_GAIN = 100.0
N_E2_HALF_PG_ML = 11.0
GNRH_PULSE_MAX_JUMP = 100.0
LH_GNRH_GAIN = 1.0
E2_GAIN = 15.0

# ZT0 is 07:00 in the shared project clock. Licinio et al. kept participants
# in darkness and asleep from 23:00 to 07:00, corresponding to ZT16--ZT24.
# Each row is (start ZT hour, end ZT hour, pulse-period multiplier).
WAKE_PULSE_PERIOD_MULTIPLIER = 1.0
SLEEP_PERIOD_BOUTS_ZT = (
    (16.0, 24.0, 1.57),
)

# First-edition endocrine SCN hypothesis. The smooth reproductive gain ranges
# from 1 to 2 and peaks at 08:00 (ZT1). Licinio fixes 08:00 at one hour after
# lights-on; Kerdelhue supports an early-morning LH-surge gate at 04:00--08:00
# clock time but did not report a ZT reference. This remains an inferred phase,
# not a direct measurement of peak human SCN-to-GnRH signaling.
SCN_ENDOCRINE_MIN_DRIVE = 1.0
SCN_ENDOCRINE_MAX_DRIVE = 2.0
SCN_ENDOCRINE_PEAK_ZT = 1.0


def _increasing_hill(value, half, exponent):
    value_power = value**exponent
    return value_power / (half**exponent + value_power)


def _decreasing_hill(value, half, exponent):
    half_power = half**exponent
    return half_power / (half_power + value**exponent)


def simulate_endocrine(
    output_h=OUTPUT_H,
    warmup_h=WARMUP_H,
    step_h=STEP_H,
    n_gain=N_GAIN,
    n_e2_half_pg_ml=N_E2_HALF_PG_ML,
    gnrh_pulse_max_jump=GNRH_PULSE_MAX_JUMP,
    lh_gnrh_gain=LH_GNRH_GAIN,
    e2_gain=E2_GAIN,
):
    """Integrate N, GnRH, LH, and E2 while retaining the pulse counter phi."""

    total_h = warmup_h + output_h
    steps = round(total_h / step_h)
    time_h = np.arange(steps + 1, dtype=float) * step_h
    scn_endocrine_drive = scn_drive(time_h)
    sleep_period_drive = sleep_period_multiplier(time_h)
    n_state = np.empty(steps + 1)
    gnrh = np.empty(steps + 1)
    lh = np.empty(steps + 1)
    e2 = np.empty(steps + 1)
    phi_trace = np.empty(steps + 1)

    gnrh[0], lh[0], e2[0], phi_trace[0] = 0.0, 5.0, 15.0, 0.0
    n_state[0] = n_gain * _decreasing_hill(
        e2[0], n_e2_half_pg_ml, HILL_COEFFICIENT
    )

    phi = 0.0
    for i in range(steps):
        n_value = n_state[i]
        period_activation = _increasing_hill(
            n_value, 1.0, HILL_COEFFICIENT
        )
        waking_period_min = PULSE_PERIOD_MAX_MIN - (
            PULSE_PERIOD_MAX_MIN - PULSE_PERIOD_MIN_MIN
        ) * period_activation
        period_min = waking_period_min * sleep_period_drive[i]

        phi_new = phi + 60.0 * step_h / period_min
        fires = np.floor(phi_new) > np.floor(phi)
        pulse = 0.0
        if fires:
            pulse = (
                gnrh_pulse_max_jump
                * _increasing_hill(
                    n_value,
                    1.0,
                    HILL_COEFFICIENT,
                )
                * scn_endocrine_drive[i]
            )
        phi = phi_new

        n_target = n_gain * _decreasing_hill(
            e2[i], n_e2_half_pg_ml, HILL_COEFFICIENT
        )
        n_state[i + 1] = n_value + step_h * (n_target - n_value)
        gnrh[i + 1] = gnrh[i] + step_h * (-K_GNRH * gnrh[i]) + pulse

        # No direct E2 -> LH inhibition: this is the Meeting 9 model decision.
        lh[i + 1] = lh[i] + step_h * (
            lh_gnrh_gain * gnrh[i] - K_LH * lh[i]
        )
        e2[i + 1] = e2[i] + step_h * (
            e2_gain * lh[i] - K_E2 * e2[i]
        )
        phi_trace[i + 1] = phi

    return {
        "time_h": time_h,
        "SCN_drive": scn_endocrine_drive,
        "sleep_period_multiplier": sleep_period_drive,
        "N": n_state,
        "phi": phi_trace,
        "GnRH": gnrh,
        "LH": lh,
        "E2": e2,
    }


# Temperature model imported from Meeting 8. The SCN-associated BAT drive is
# smallest at the ZT7 SCN neural peak and largest 12 hours later, at ZT19.
SCN_NEURAL_PEAK_ZT = 7.0

# Thermal measurements and literature-informed scales.
AMBIENT_TEMPERATURE_C = 21.25
INITIAL_TEMPERATURE_C = 37.0

ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G = 22.36
ZENGIN_SHORT_TERM_WT_OVX_FAT_MASS_G = 3.98
ZENGIN_SHORT_TERM_WT_OVX_FAT_FRACTION = (
    ZENGIN_SHORT_TERM_WT_OVX_FAT_MASS_G
    / ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
)

MOUSE_SPECIFIC_HEAT_CAL_G_C = (
    0.9100 - 0.5231 * ZENGIN_SHORT_TERM_WT_OVX_FAT_FRACTION
)
BODY_HEAT_CAPACITY_KCAL_C = (
    MOUSE_SPECIFIC_HEAT_CAL_G_C
    * ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
    / 1000.0
)

# The digitized OVX resting-heat estimate is approximately 0.32 kcal/h;
# Meeting 9 now uses that central estimate rather than the former 0.31 lower
# calibration value.
OVX_BACKGROUND_HEAT_KCAL_H = 0.32
OVX_MEAN_RESTING_HEAT_KCAL_H = OVX_BACKGROUND_HEAT_KCAL_H

ABREU_VIEIRA_REFERENCE_BODY_MASS_G = 26.97
ABREU_VIEIRA_REFERENCE_CONDUCTANCE_KCAL_H_C = 0.0255
WHOLE_BODY_CONDUCTANCE_KCAL_H_C = (
    ABREU_VIEIRA_REFERENCE_CONDUCTANCE_KCAL_H_C
    * (
        ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
        / ABREU_VIEIRA_REFERENCE_BODY_MASS_G
    ) ** (2.0 / 3.0)
)

# E2 and SCN share one whole-body thermogenic Hill response. The 14% ceiling
# is the upper Camporez-supported whole-body thermogenic fraction. K=30 and
# n=5 are model-shape parameters, not a measured receptor EC50 or
# cooperativity coefficient. At E2=39 pg/mL with zero added SCN drive, this
# calibration produces approximately 11% extra heat relative to background.
E2_SCN_WHOLE_BODY_THERMOGENIC_FRACTION = 0.14
E2_HEAT_HALF_PG_ML = 30.0
E2_HEAT_HILL_COEFFICIENT = 5.0

# Convert the dimensionless 0--1 SCN-associated BAT drive into E2-equivalent
# pg/mL before summing the inputs. This remains the parameter intended for a
# sweep; its current maximum contribution is 11 pg/mL.
SCN_E2_EQUIVALENT_SCALE_PG_ML = 11.0

# Activity remains a separate temperature-equivalent input.
ACTIVITY_TEMPERATURE_DIFFERENCE_C = 1.9
ACTIVITY_BOUTS_ZT = (
    (12.0, 17.0, 1.0),
)

THERMAL_RATE_CONSTANT_H = (
    WHOLE_BODY_CONDUCTANCE_KCAL_H_C / BODY_HEAT_CAPACITY_KCAL_C
)
BASAL_TEMPERATURE_RATE_C_H = (
    OVX_BACKGROUND_HEAT_KCAL_H / BODY_HEAT_CAPACITY_KCAL_C
)
E2_SCN_WHOLE_BODY_THERMOGENIC_GAIN_C_H = (
    E2_SCN_WHOLE_BODY_THERMOGENIC_FRACTION
    * OVX_BACKGROUND_HEAT_KCAL_H
    / BODY_HEAT_CAPACITY_KCAL_C
)


def zeitgeber_hour(time_h):
    """Map simulation time to the repeating ZT0--ZT24 clock."""

    return np.mod(time_h, ZT_PERIOD_H)


def sleep_period_multiplier(zt_h):
    """Return 1 awake and 1.57 during the Licinio sleep window."""

    phase_h = zeitgeber_hour(zt_h)
    multiplier = np.full_like(
        phase_h,
        WAKE_PULSE_PERIOD_MULTIPLIER,
        dtype=float,
    )

    for start_h, end_h, height in SLEEP_PERIOD_BOUTS_ZT:
        if start_h <= end_h:
            in_bout = (phase_h >= start_h) & (phase_h < end_h)
        else:
            in_bout = (phase_h >= start_h) | (phase_h < end_h)
        multiplier = np.where(in_bout, height, multiplier)

    return multiplier


def scn_drive(zt_h):
    """Return the smooth 1--2 reproductive gain, peaking at ZT1/08:00."""

    midpoint = (SCN_ENDOCRINE_MAX_DRIVE + SCN_ENDOCRINE_MIN_DRIVE) / 2.0
    amplitude = (SCN_ENDOCRINE_MAX_DRIVE - SCN_ENDOCRINE_MIN_DRIVE) / 2.0
    drive = midpoint + amplitude * np.cos(
        2.0
        * np.pi
        * (zeitgeber_hour(zt_h) - SCN_ENDOCRINE_PEAK_ZT)
        / ZT_PERIOD_H
    )

    return drive


def e2_scn_heat_hill(
    e2_pg_ml,
    scn,
    scn_e2_scale_pg_ml=SCN_E2_EQUIVALENT_SCALE_PG_ML,
):
    """Return the K=30, n=5 shared E2-plus-SCN heat response."""

    drive = e2_pg_ml + scn_e2_scale_pg_ml * scn
    numerator = drive**E2_HEAT_HILL_COEFFICIENT
    return numerator / (
        E2_HEAT_HALF_PG_ML**E2_HEAT_HILL_COEFFICIENT + numerator
    )


def scn_bat_drive(zt_h):
    """Return the smooth SCN-associated BAT drive in [0, 1]."""

    return (
        1.0
        - np.cos(
            2.0
            * np.pi
            * (zeitgeber_hour(zt_h) - SCN_NEURAL_PEAK_ZT)
            / ZT_PERIOD_H
        )
    ) / 2.0


def activity_drive(zt_h):
    """Return the editable piecewise-constant activity schedule."""

    phase_h = zeitgeber_hour(zt_h)
    activity = np.zeros_like(phase_h, dtype=float)

    for start_h, end_h, height in ACTIVITY_BOUTS_ZT:
        if start_h <= end_h:
            in_bout = (phase_h >= start_h) & (phase_h < end_h)
        else:
            in_bout = (phase_h >= start_h) | (phase_h < end_h)
        activity = np.where(in_bout, height, activity)

    return activity


def temperature_rate(
    temperature_c,
    e2_pg_ml,
    scn,
    activity,
    activity_multiplier=1.0,
    scn_e2_scale_pg_ml=SCN_E2_EQUIVALENT_SCALE_PG_ML,
):
    """Evaluate the combined E2-SCN heat-balance ODE in degrees C/hour."""

    combined_e2_scn_thermogenesis = (
        E2_SCN_WHOLE_BODY_THERMOGENIC_GAIN_C_H
        * e2_scn_heat_hill(e2_pg_ml, scn, scn_e2_scale_pg_ml)
    )
    activity_thermal_effect = (
        THERMAL_RATE_CONSTANT_H
        * activity_multiplier
        * ACTIVITY_TEMPERATURE_DIFFERENCE_C
        * activity
    )
    heat_loss = THERMAL_RATE_CONSTANT_H * (
        temperature_c - AMBIENT_TEMPERATURE_C
    )
    return (
        BASAL_TEMPERATURE_RATE_C_H
        + combined_e2_scn_thermogenesis
        + activity_thermal_effect
        - heat_loss
    )


def simulate_temperature(
    time_h,
    e2,
    activity_multiplier=1.0,
    scn_e2_scale_pg_ml=SCN_E2_EQUIVALENT_SCALE_PG_ML,
):
    """Return the temperature trajectory on the endocrine time grid."""

    if len(time_h) != len(e2):
        raise ValueError("time_h and e2 must have the same length")

    temperature = np.empty(len(time_h))
    temperature[0] = INITIAL_TEMPERATURE_C
    zt_h = zeitgeber_hour(time_h)
    scn = scn_bat_drive(zt_h)
    activity = activity_drive(zt_h)

    for i in range(len(time_h) - 1):
        step_h = time_h[i + 1] - time_h[i]
        temperature[i + 1] = temperature[i] + step_h * temperature_rate(
            temperature[i],
            e2[i],
            scn[i],
            activity[i],
            activity_multiplier,
            scn_e2_scale_pg_ml,
        )

    return temperature


def run(
    output_h=OUTPUT_H,
    scn_e2_scale_pg_ml=SCN_E2_EQUIVALENT_SCALE_PG_ML,
):
    """Return the retained endocrine and temperature trajectories."""

    result = simulate_endocrine(output_h=output_h)
    thermal = simulate_temperature(
        result["time_h"],
        result["E2"],
        scn_e2_scale_pg_ml=scn_e2_scale_pg_ml,
    )
    output_points = round(output_h / STEP_H) + 1
    trajectory = {
        "hour": result["time_h"][-output_points:] - WARMUP_H,
        "zt_h": zeitgeber_hour(result["time_h"][-output_points:]),
        "N": result["N"][-output_points:],
        "phi": result["phi"][-output_points:],
        "GnRH": result["GnRH"][-output_points:],
        "LH": result["LH"][-output_points:],
        "LH_IU_L": result["LH"][-output_points:],
        "E2_pg_ml": result["E2"][-output_points:],
        "SCN_drive": result["SCN_drive"][-output_points:],
        "sleep_period_multiplier": result["sleep_period_multiplier"][
            -output_points:
        ],
        "scn_e2_scale_pg_ml": scn_e2_scale_pg_ml,
    }
    trajectory["temperature_c"] = thermal[-output_points:]
    return trajectory


if __name__ == "__main__":
    trajectory = run()
    e2 = trajectory["E2_pg_ml"]
    lh = trajectory["LH_IU_L"]
    temperature = trajectory["temperature_c"]
    print(
        "E2 min/mean/max: "
        f"{e2.min():.3f}/{e2.mean():.3f}/{e2.max():.3f} pg/mL"
    )
    print(
        "LH min/mean/max: "
        f"{lh.min():.3f}/{lh.mean():.3f}/{lh.max():.3f} IU/L"
    )
    print(
        "Temperature min/mean/max: "
        f"{temperature.min():.3f}/{temperature.mean():.3f}/"
        f"{temperature.max():.3f} degrees C"
    )
