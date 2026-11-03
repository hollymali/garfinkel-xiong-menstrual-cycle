"""Final meeting-7 E2-SCN-activity mouse core-temperature model.

``run()`` returns the final 72-hour time series after a warm-up period.  The
module intentionally contains no plotting or file-export code, matching the
style of the project's top-level ``MODEL.py``.  It returns the cycling-E2
non-pregnant condition, a fixed 13.3-pg/mL aged-female comparison, a
constant-high-E2 counterfactual, and a constant-high-E2 pregnancy comparison
with progesterone and reduced activity.

The implemented temperature balance is

    dT/dt = B/Cb + (0.11 B/Cb) H_E2(E2)
           + kappa [0.3 SCN_T + 1.9 m_A Activity + m_P4 DeltaT_P4]
           - kappa (T - Ta),

where kappa = G/Cb.  Estradiol does not modify whole-body conductance.

With the thermal Hill gate n=5 and K=9.2 pg/mL, the aged-female comparison
trough is 0.375 C above the cycling trough while its peak is 0.113 C lower.
The constant-high-E2 trough is 0.584 C above cycling while its peak is 0.096 C
higher. No pregnancy-specific basal-heat term is used.
"""

import numpy as np


# Shared simulation grid (hours): one 6-second step for both subsystems.
STEP_H = 1.0 / 600.0
ZT_PERIOD_H = 24.0
LIGHTS_ON_ZT = 0.0
LIGHTS_OFF_ZT = 12.0
WARMUP_H = 4.0 * ZT_PERIOD_H
OUTPUT_H = 72.0

# Endocrine model
# Constants use hours and pg/mL unless otherwise noted.
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


# Temperature model
# SCN-VP firing peaks near ZT7, but the AVP temperature output is hypothermic.
# The thermal drive is therefore the negative of the ZT7-centered neural
# waveform: -1 at ZT7 (mid-light), +1 at ZT19 (mid-dark), and zero mean.
SCN_NEURAL_PEAK_ZT = ENDOCRINE_CLOCK_PEAK_ZT

# Thermal measurements and literature-derived scales
AMBIENT_TEMPERATURE_C = 22.0
INITIAL_TEMPERATURE_C = 37.0
BODY_HEAT_CAPACITY_KCAL_C = 0.018463240021
# Approximate short-term (14-18 d post-OVX) WT OVX resting metabolic rate
# read from Zengin et al. (2013), Figure 5A. This is resting heat at 22 C,
# not classical thermoneutral, postabsorptive BMR.
OVX_MEAN_RESTING_HEAT_KCAL_H = 0.3226854877
# Female mass-scaled value informed by the 0.0255 kcal/(h*C) light-phase
# TEE-PAEE slope measured in 26.97-g male C57BL/6J mice (Abreu-Vieira et al.,
# 2015, Table 1 and Table S1). This is a calibrated literature-informed value,
# not a direct female measurement.
WHOLE_BODY_CONDUCTANCE_KCAL_H_C = 0.0232

# Reported effect sizes used by the reduced thermal model. This E2 term is
# whole-body thermogenesis; it is not assigned exclusively to BAT.
E2_WHOLE_BODY_THERMOGENIC_FRACTION = 0.11
# Shape-calibrated thermogenesis gate. K lies near the low end of the 8.3-39.1
# pg/mL serum-E2 interval bracketed by the two groups in Camporez et al. The
# steep gate raises the 13.3 pg/mL aged-female comparison trough without
# materially separating the nonpregnant peaks. Neither K nor n is a measured
# thermal EC50/cooperativity coefficient.
E2_HEAT_HALF_PG_ML = 9.2
E2_HEAT_HILL_COEFFICIENT = 5.0
# Half of the approximately 0.6 C figure-informed light/dark span. Because B
# is a 24-h mean, the SCN temperature-equivalent drive is centered on zero.
SCN_TEMPERATURE_AMPLITUDE_C = 0.3
ACTIVITY_TEMPERATURE_DIFFERENCE_C = 1.9

# Calibrated defended-temperature increase within the approximately 0.3-0.4 C
# active-phase effect of 10 mg/21 d progesterone replacement in gonadectomized
# female C57BL/6 mice (Cintron-Colon et al., 2019). This term is enabled only
# in pregnancy, not in the cycling-E2 or constant-high-E2-only equations.
P4_TEMPERATURE_SHIFT_C = 0.343987915

# Edit this tuple to change the rectangular activity bouts used by the model.
# Each row is: (start ZT hour, end ZT hour, relative activity height).
ACTIVITY_BOUTS_ZT = (
    (12.0, 17.0, 1.0),
    # (21.5, 23.5, 0.8),  # Example: uncomment for a second dark-phase bout.
)

# Pregnancy comparison. The activity-associated temperature adjustment retains
# 70% of the nonpregnant input. This is deliberately conservative: Ladyman et
# al. report an approximately 50% early-pregnancy reduction in voluntary wheel
# running, while Clarke et al. report an approximately 70% reduction in total
# activity, especially in the dark phase. Those measurements do not translate
# one-to-one into activity heat, so 0.70 remains an explicit modeling choice.
# A formerly used 0.51 multiplier on the entire signed SCN term has been
# removed: inducible BAT capacity does not measure the full circadian
# temperature influence represented here.
PREGNANT_E2_PG_ML = 50.88
# Glassberg et al. (2014): intact 24-month-old female C57BL/6J mice,
# 13.3 +/- 0.6 pg/mL (group mean, n=9). This is a fixed comparison input,
# not a continuously measured time series.
AGED_CONSTANT_E2_PG_ML = 13.3
PREGNANT_ACTIVITY_MULTIPLIER = 0.7
PREGNANT_P4_MULTIPLIER = 1.0

# Derived ODE coefficients (degrees Celsius per hour where applicable)
# Thermal rate constant kappa = G / C_b, with units h^-1. It converts the
# temperature-equivalent SCN, activity, P4, and loss terms to degrees C/hour.
THERMAL_RATE_CONSTANT_H = (
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


def e2_heat_hill(e2_pg_ml):
    """Return the calibrated increasing whole-body E2 heat response."""
    numerator = e2_pg_ml**E2_HEAT_HILL_COEFFICIENT
    return numerator / (
        E2_HEAT_HALF_PG_ML**E2_HEAT_HILL_COEFFICIENT
        + numerator
    )


def scn_drive(zt_h):
    """Return the zero-mean SCN temperature drive in [-1, 1].

    SCN-VP neural firing peaks at ZT7, but its AVP output lowers temperature;
    hence the minus sign. The thermal drive troughs at ZT7 and peaks at ZT19.
    """
    return -np.cos(
        2.0
        * np.pi
        * (zeitgeber_hour(zt_h) - SCN_NEURAL_PEAK_ZT)
        / ZT_PERIOD_H
    )


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


def temperature_rate(
    temperature_c,
    e2_pg_ml,
    scn,
    activity,
    activity_multiplier=1.0,
    p4_multiplier=0.0,
):
    """Evaluate the reduced heat-balance ODE in degrees Celsius per hour."""
    e2_whole_body_thermogenesis = (
        E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H * e2_heat_hill(e2_pg_ml)
    )
    scn_thermal_effect = (
        THERMAL_RATE_CONSTANT_H * SCN_TEMPERATURE_AMPLITUDE_C * scn
    )
    activity_thermal_effect = (
        THERMAL_RATE_CONSTANT_H
        * activity_multiplier
        * ACTIVITY_TEMPERATURE_DIFFERENCE_C
        * activity
    )
    p4_thermal_effect = (
        THERMAL_RATE_CONSTANT_H
        * p4_multiplier
        * P4_TEMPERATURE_SHIFT_C
    )
    heat_loss = THERMAL_RATE_CONSTANT_H * (
        temperature_c - AMBIENT_TEMPERATURE_C
    )
    return (
        BASAL_TEMPERATURE_RATE_C_H
        + e2_whole_body_thermogenesis
        + scn_thermal_effect
        + activity_thermal_effect
        + p4_thermal_effect
        - heat_loss
    )


def simulate_temperature(
    time_h,
    e2,
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
            activity_multiplier,
            p4_multiplier,
        )

    e2_whole_body_thermogenesis = (
        E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H * e2_heat_hill(e2)
    )
    scn_thermal_effect = (
        THERMAL_RATE_CONSTANT_H * SCN_TEMPERATURE_AMPLITUDE_C * scn
    )
    activity_thermal_effect = (
        THERMAL_RATE_CONSTANT_H
        * activity_multiplier
        * ACTIVITY_TEMPERATURE_DIFFERENCE_C
        * activity
    )
    p4_thermal_effect = np.full(
        len(time_h),
        THERMAL_RATE_CONSTANT_H * p4_multiplier * P4_TEMPERATURE_SHIFT_C,
    )
    heat_loss = (
        THERMAL_RATE_CONSTANT_H
        * (temperature - AMBIENT_TEMPERATURE_C)
    )

    return {
        "temperature_c": temperature,
        "scn_drive": scn,
        "activity_drive": activity,
        "basal_rate_c_h": np.full(len(time_h), BASAL_TEMPERATURE_RATE_C_H),
        "e2_whole_body_thermogenesis_c_h": e2_whole_body_thermogenesis,
        "scn_thermal_effect_c_h": scn_thermal_effect,
        "activity_thermal_effect_c_h": activity_thermal_effect,
        "p4_thermal_effect_c_h": p4_thermal_effect,
        "heat_loss_c_h": heat_loss,
    }


def run():
    """Return 72 hours of all four final temperature comparisons."""
    endocrine = simulate_endocrine()
    output_points = round(OUTPUT_H / STEP_H) + 1
    output_slice = slice(-output_points, None)
    nonpregnant_thermal = simulate_temperature(
        endocrine["time_h"], endocrine["E2"]
    )
    constant_low_e2_value = AGED_CONSTANT_E2_PG_ML
    constant_low_e2 = np.full_like(
        endocrine["E2"], constant_low_e2_value
    )
    constant_low_thermal = simulate_temperature(
        endocrine["time_h"], constant_low_e2
    )
    pregnant_e2 = np.full_like(endocrine["E2"], PREGNANT_E2_PG_ML)
    constant_high_thermal = simulate_temperature(
        endocrine["time_h"], pregnant_e2
    )
    pregnant_thermal = simulate_temperature(
        endocrine["time_h"],
        pregnant_e2,
        activity_multiplier=PREGNANT_ACTIVITY_MULTIPLIER,
        p4_multiplier=PREGNANT_P4_MULTIPLIER,
    )
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
    result["constant_low_E2_pg_ml"] = constant_low_e2[output_slice]
    result.update({
        f"constant_low_{name}": values[output_slice]
        for name, values in constant_low_thermal.items()
    })
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
