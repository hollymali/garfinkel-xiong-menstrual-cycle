"""Reduced E2-SCN-activity mouse core-temperature model.

``run()`` returns the final 72-hour time series after a warm-up period.  The
module intentionally contains no plotting or file-export code, matching the
style of the project's top-level ``MODEL.py``.  It returns the cycling-E2
non-pregnant condition, a constant-high-E2-only counterfactual, and a
constant-high-E2 pregnancy comparison with progesterone and reduced SCN-BAT
coupling.
"""

import numpy as np


# Shared simulation grid (hours): one 6-second step for both subsystems.
STEP_H = 1.0 / 600.0
ZT_PERIOD_H = 24.0
LIGHTS_ON_ZT = 0.0
LIGHTS_OFF_ZT = 12.0
WARMUP_H = 4.0 * ZT_PERIOD_H
OUTPUT_H = 72.0

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
E2_RESPONSE_HALF = 20.1
PULSE_AMPLITUDE = 6.8
N_PULSE_HALF = 0.28
PULSE_PERIOD_LOW_H = 1.05
PULSE_PERIOD_HIGH_H = 1.90
# All circadian phases are expressed in Zeitgeber time. Thirouin et al. (2023)
# used 07:00 lights-on and found peak female SCN-VP firing at CT7 (14:00).
ENDOCRINE_CLOCK_PEAK_ZT = 7.0
ENDOCRINE_CLOCK_GAIN = 0.7

# The thermal SCN term uses the same ZT7-centered clock as the endocrine
# model. The clock is set to zero before ZT12 and preserved unchanged after
# ZT12. Its negative values are suppression of heat production, not negative
# SCN firing. The retained signal is not mean-centered.
SCN_RHYTHM_PEAK_ZT = ENDOCRINE_CLOCK_PEAK_ZT

# Thermal measurements and literature-derived scales
AMBIENT_TEMPERATURE_C = 22.0
INITIAL_TEMPERATURE_C = 37.0
BODY_HEAT_CAPACITY_KCAL_C = 0.018379
# Approximate short-term (14-18 d post-OVX) WT OVX resting metabolic rate
# read from Zengin et al. (2013), Figure 5A. This is resting heat at 22 C,
# not classical thermoneutral, postabsorptive BMR.
OVX_MEAN_RESTING_HEAT_KCAL_H = 0.32
# Provisional mass-adjusted conductance anchored to the 0.0255 kcal/(h*C)
# light-phase TEE-PAEE slope measured in 26.97-g chow-fed male C57BL/6J mice
# (Abreu-Vieira et al., 2015, Table 1 and Table S1). A value of 0.02425 is the
# approximate M^(2/3) scaling for a 25-g mouse; replace it when the matched
# female body mass is available.
WHOLE_BODY_CONDUCTANCE_KCAL_H_C = 0.02425

# Reported effect sizes used by the reduced thermal model. This E2 term is
# whole-body thermogenesis; it is not assigned exclusively to BAT.
E2_WHOLE_BODY_THERMOGENIC_FRACTION = 0.11
# Shape-calibrated thermogenesis gate. The midpoint is constrained to the
# 8.3-39.1 pg/mL interval bracketed by the two serum-E2 groups in Camporez et
# al.; neither this midpoint nor the exponent is a measured thermal EC50.
E2_HEAT_HALF_PG_ML = 9.0
E2_HEAT_HILL_COEFFICIENT = 1.41
E2_LOSS_HALF_PG_ML = 7.49
E2_LOSS_HILL_COEFFICIENT = 1.83
# Half of the approximately 0.6 C figure-informed light/dark span is assigned
# as a one-directional suppression from the basal reference. Because the raw
# retained cosine reaches -1, this gives a maximum -0.3 C-equivalent effect.
SCN_SUPPRESSION_SCALE_C = 0.3
ACTIVITY_TEMPERATURE_DIFFERENCE_C = 1.9
# Maximum whole-body heat-loss fraction conserved by E2 through the tail:
# 0.55 figure-estimated suppression of rat tail HLI in Dacks and Rance,
# multiplied by the approximately 0.05 mouse-tail share from Skop et al.
E2_HEAT_CONSERVATION_FRACTION = 0.0275

# Conservative defended-temperature increase from 10 mg/21 d progesterone
# replacement in gonadectomized female C57BL/6 mice (Cintron-Colon et al.,
# 2019).  This term is enabled only in the pregnancy comparison; it is not
# included in the cycling-E2 or constant-high-E2-only equations.
P4_TEMPERATURE_SHIFT_C = 0.3

# Edit this tuple to change the rectangular activity bouts used by the model.
# Each row is: (start ZT hour, end ZT hour, relative activity height).
ACTIVITY_BOUTS_ZT = (
    (12.0, 17.0, 1.0),
    # (21.5, 23.5, 0.8),  # Example: uncomment for a second dark-phase bout.
)

# Pregnancy comparison. Pregnancy retains 51% of inducible BAT capacity
# (49% reduction at GD14). The activity-associated heat drive is calibrated
# to 73.22% of the nonpregnant drive, a conservative change relative to the
# larger pregnancy-associated reduction reported for measured locomotion.
PREGNANT_E2_PG_ML = 50.88
PREGNANT_SCN_BAT_MULTIPLIER = 0.51
PREGNANT_ACTIVITY_MULTIPLIER = 0.7322
PREGNANT_P4_MULTIPLIER = 1.0

# Derived ODE coefficients (degrees Celsius per hour where applicable)
THERMAL_FLOW_RATE_H = (
    WHOLE_BODY_CONDUCTANCE_KCAL_H_C / BODY_HEAT_CAPACITY_KCAL_C
)
BASAL_TEMPERATURE_RATE_C_H = (
    OVX_MEAN_RESTING_HEAT_KCAL_H / BODY_HEAT_CAPACITY_KCAL_C
)
E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H = (
    E2_WHOLE_BODY_THERMOGENIC_FRACTION
    * OVX_MEAN_RESTING_HEAT_KCAL_H
    / BODY_HEAT_CAPACITY_KCAL_C
)
SCN_THERMAL_GAIN_C_H = (
    THERMAL_FLOW_RATE_H * SCN_SUPPRESSION_SCALE_C
)
ACTIVITY_THERMOGENIC_GAIN_C_H = (
    THERMAL_FLOW_RATE_H * ACTIVITY_TEMPERATURE_DIFFERENCE_C
)
P4_THERMAL_GAIN_C_H = (
    THERMAL_FLOW_RATE_H * P4_TEMPERATURE_SHIFT_C
)


def zeitgeber_hour(time_h):
    """Map simulation time to the shared repeating ZT0-ZT24 clock."""
    return np.mod(time_h, ZT_PERIOD_H)


def scn_circadian_wave(zt_h):
    """Return the endocrine ZT7-centered cosine in [-1, 1]."""
    phase_h = zeitgeber_hour(zt_h)
    return np.cos(
        2.0
        * np.pi
        * (phase_h - ENDOCRINE_CLOCK_PEAK_ZT)
        / ZT_PERIOD_H
    )


def endocrine_clock(zt_h):
    """Return the SCN-linked endocrine multiplier, peaking at ZT7."""
    unit_clock = 0.5 * (1.0 + scn_circadian_wave(zt_h))
    return 1.0 + ENDOCRINE_CLOCK_GAIN * unit_clock


def e2_heat_hill(e2_pg_ml):
    """Return the calibrated increasing whole-body E2 heat response."""
    numerator = e2_pg_ml**E2_HEAT_HILL_COEFFICIENT
    return numerator / (
        E2_HEAT_HALF_PG_ML**E2_HEAT_HILL_COEFFICIENT
        + numerator
    )


def e2_loss_hill(e2_pg_ml):
    """Return the increasing E2 tail heat-conservation response."""
    numerator = e2_pg_ml**E2_LOSS_HILL_COEFFICIENT
    return numerator / (
        E2_LOSS_HALF_PG_ML**E2_LOSS_HILL_COEFFICIENT
        + numerator
    )


def scn_retained_wave(zt_h):
    """Return zero before ZT12 and the unchanged ZT7 cosine afterward."""
    phase_h = zeitgeber_hour(zt_h)
    return np.where(
        phase_h >= LIGHTS_OFF_ZT,
        np.cos(
            2.0 * np.pi * (phase_h - SCN_RHYTHM_PEAK_ZT) / ZT_PERIOD_H
        ),
        0.0,
    )


def scn_drive(zt_h):
    """Return the uncentered, time-cut SCN cosine thermal adjustment."""
    return scn_retained_wave(zt_h)


def activity_drive(zt_h):
    """Return the editable piecewise-constant activity schedule."""
    phase_h = zeitgeber_hour(zt_h)
    activity = np.zeros_like(phase_h, dtype=float)

    for start_h, end_h, height in ACTIVITY_BOUTS_ZT:
        if start_h <= end_h:
            in_bout = (phase_h >= start_h) & (phase_h < end_h)
        else:
            # This permits a bout such as (22.0, 2.0, 0.5) across ZT0.
            in_bout = (phase_h >= start_h) | (phase_h < end_h)
        activity = np.where(in_bout, height, activity)

    return activity


def simulate_endocrine():
    """Integrate N -> GnRH -> LH -> E2 with an explicit Euler loop."""
    total_h = WARMUP_H + OUTPUT_H
    steps = round(total_h / STEP_H)
    time_h = np.arange(steps + 1) * STEP_H
    n_state = np.empty(steps + 1)
    gnrh = np.empty(steps + 1)
    lh = np.empty(steps + 1)
    e2 = np.empty(steps + 1)

    n_state[0] = E2_RESPONSE_HALF**5 / (
        E2_RESPONSE_HALF**5 + 15.0**5
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
        zt_h = zeitgeber_hour(time_h[i])
        pulse = (
            PULSE_AMPLITUDE * amplitude_drive * endocrine_clock(zt_h)
            if fires
            else 0.0
        )
        phase = new_phase

        n_state[i + 1] = n_value + STEP_H * (
            N_RATE
            * E2_RESPONSE_HALF**5
            / (E2_RESPONSE_HALF**5 + e2[i] ** 5)
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
    temperature_c,
    e2_pg_ml,
    scn,
    activity,
    scn_bat_multiplier=1.0,
    activity_multiplier=1.0,
    p4_multiplier=0.0,
):
    """Evaluate the reduced heat-balance ODE in degrees Celsius per hour."""
    e2_whole_body_thermogenesis = (
        E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H * e2_heat_hill(e2_pg_ml)
    )
    scn_thermal_effect = (
        scn_bat_multiplier * SCN_THERMAL_GAIN_C_H * scn
    )
    activity_thermogenesis = (
        activity_multiplier * ACTIVITY_THERMOGENIC_GAIN_C_H * activity
    )
    p4_thermal_effect = p4_multiplier * P4_THERMAL_GAIN_C_H
    baseline_heat_loss = (
        THERMAL_FLOW_RATE_H
        * (temperature_c - AMBIENT_TEMPERATURE_C)
    )
    e2_heat_conservation = (
        THERMAL_FLOW_RATE_H
        * E2_HEAT_CONSERVATION_FRACTION
        * e2_loss_hill(e2_pg_ml)
        * (temperature_c - AMBIENT_TEMPERATURE_C)
    )
    return (
        BASAL_TEMPERATURE_RATE_C_H
        + e2_whole_body_thermogenesis
        + scn_thermal_effect
        + activity_thermogenesis
        + p4_thermal_effect
        - baseline_heat_loss
        + e2_heat_conservation
    )


def simulate_temperature(
    time_h,
    e2,
    scn_bat_multiplier=1.0,
    activity_multiplier=1.0,
    p4_multiplier=0.0,
):
    """Integrate the temperature ODE on the shared endocrine time grid."""
    temperature = np.empty(len(time_h))
    temperature[0] = INITIAL_TEMPERATURE_C
    zt_h = zeitgeber_hour(time_h)
    scn = scn_drive(zt_h)
    activity = activity_drive(zt_h)

    for i in range(len(time_h) - 1):
        temperature[i + 1] = temperature[i] + STEP_H * temperature_rate(
            temperature[i],
            e2[i],
            scn[i],
            activity[i],
            scn_bat_multiplier,
            activity_multiplier,
            p4_multiplier,
        )

    e2_whole_body_thermogenesis = (
        E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H * e2_heat_hill(e2)
    )
    scn_thermal_effect = (
        scn_bat_multiplier * SCN_THERMAL_GAIN_C_H * scn
    )
    activity_thermogenesis = (
        activity_multiplier * ACTIVITY_THERMOGENIC_GAIN_C_H * activity
    )
    p4_thermal_effect = np.full(
        len(time_h), p4_multiplier * P4_THERMAL_GAIN_C_H
    )
    baseline_heat_loss = (
        THERMAL_FLOW_RATE_H
        * (temperature - AMBIENT_TEMPERATURE_C)
    )
    e2_heat_conservation = (
        THERMAL_FLOW_RATE_H
        * E2_HEAT_CONSERVATION_FRACTION
        * e2_loss_hill(e2)
        * (temperature - AMBIENT_TEMPERATURE_C)
    )
    net_heat_loss = baseline_heat_loss - e2_heat_conservation

    return {
        "temperature_c": temperature,
        "scn_drive": scn,
        "activity_drive": activity,
        "basal_rate_c_h": np.full(len(time_h), BASAL_TEMPERATURE_RATE_C_H),
        "e2_whole_body_thermogenesis_c_h": e2_whole_body_thermogenesis,
        "scn_thermal_effect_c_h": scn_thermal_effect,
        "activity_thermogenesis_c_h": activity_thermogenesis,
        "p4_thermal_effect_c_h": p4_thermal_effect,
        "baseline_heat_loss_c_h": baseline_heat_loss,
        "e2_heat_conservation_c_h": e2_heat_conservation,
        "heat_loss_c_h": net_heat_loss,
    }


def run():
    """Return 72 hours of non-pregnant and pregnancy-comparison variables."""
    endocrine = simulate_endocrine()
    nonpregnant_thermal = simulate_temperature(
        endocrine["time_h"], endocrine["E2"]
    )
    pregnant_e2 = np.full_like(endocrine["E2"], PREGNANT_E2_PG_ML)
    constant_high_thermal = simulate_temperature(
        endocrine["time_h"], pregnant_e2
    )
    pregnant_thermal = simulate_temperature(
        endocrine["time_h"],
        pregnant_e2,
        scn_bat_multiplier=PREGNANT_SCN_BAT_MULTIPLIER,
        activity_multiplier=PREGNANT_ACTIVITY_MULTIPLIER,
        p4_multiplier=PREGNANT_P4_MULTIPLIER,
    )
    output_points = round(OUTPUT_H / STEP_H) + 1
    output_slice = slice(-output_points, None)

    result = {
        "hour": endocrine["time_h"][output_slice],
        "zt_h": zeitgeber_hour(endocrine["time_h"][output_slice]),
        "N": endocrine["N"][output_slice],
        "GnRH": endocrine["GnRH"][output_slice],
        "LH": endocrine["LH"][output_slice],
        "E2_pg_ml": endocrine["E2"][output_slice],
    }
    # Preserve the original unprefixed fields as the non-pregnant outputs.
    result.update({
        name: values[output_slice]
        for name, values in nonpregnant_thermal.items()
    })
    result["nonpregnant_temperature_c"] = result["temperature_c"]
    result["constant_high_E2_pg_ml"] = pregnant_e2[output_slice]
    result.update({
        f"constant_high_{name}": values[output_slice]
        for name, values in constant_high_thermal.items()
    })
    result["pregnant_E2_pg_ml"] = pregnant_e2[output_slice]
    result.update({
        f"pregnant_{name}": values[output_slice]
        for name, values in pregnant_thermal.items()
    })
    result["hour"] = result["hour"] - result["hour"][0]
    return result
