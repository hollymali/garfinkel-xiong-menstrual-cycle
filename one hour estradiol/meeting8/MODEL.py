"""Final Meeting 8 E2-SCN-activity mouse core-temperature model.

``run()`` returns the final 72-hour time series after a 96-hour warm-up. The
model has four conditions: cycling-E2/nonpregnant, fixed 13.3-pg/mL aged E2,
constant-high E2, and constant-high E2 with reduced activity.

The implemented temperature balance is

    dT/dt = B/Cb + (0.14 B/Cb) H_E2(E2)
           + kappa [0.19 SCN_T + 1.9 m_A Activity]
           - kappa (T - Ta),

where kappa = G/Cb and

    H_E2(E2) = E2**5 / (11**5 + E2**5).

Estradiol does not modify whole-body conductance. The cycling E2 driver uses
one E2 state, literature-scale hormone half-lives, and a pulse generator
retuned to the frequently sampled Licinio mid-to-late-follicular waveform.
The pulse generator is fully deterministic: aggregate KNDy/kisspeptin drive
sets the period directly and the phi counter schedules events. Pulse amplitude
has no direct clock term and no random component. Instantaneous E2 applies a
standard inhibitory Hill modulation to pituitary GnRH-to-LH gain, with no
delayed or averaged E2 state. Its 21-pg/mL half-saturation lies inside the
21-58-pg/mL Shaw calibration interval. The remaining unanchored HPG constants
are calibrated against the digitized Licinio waveform after removal of the
endocrine clock envelope.
No pregnancy-specific basal-heat or SCN multiplier is used in the
reduced-activity comparison.
"""

import numpy as np


# Shared simulation grid (hours): one 6-second step for both subsystems.
STEP_H = 1.0 / 600.0
ZT_PERIOD_H = 24.0
LIGHTS_ON_ZT = 0.0
LIGHTS_OFF_ZT = 12.0
WARMUP_H = 4.0 * ZT_PERIOD_H
OUTPUT_H = 72.0


# Endocrine model. Constants use hours and pg/mL unless otherwise noted.
GNRH_HALF_LIFE_H = 1.0 / 15.0 # 4 min
LH_HALF_LIFE_H = 1.0 / 3.0 # 20 min
E2_HALF_LIFE_H = 1.0 / 2.0 # 30 min
K_GNRH = np.log(2.0) / GNRH_HALF_LIFE_H
K_LH = np.log(2.0) / LH_HALF_LIFE_H
K_E2 = np.log(2.0) / E2_HALF_LIFE_H
E2_GAIN = 2.5
N_TURNOVER_RATE_H = 1.85
N_MAX = 100.0
N_E2_HALF_PG_ML = 5.6
N_E2_HILL_COEFFICIENT = 4.6

# Unanchored shape constants for the direct N-to-period and N-to-amplitude
# responses. The pulse-period limits themselves remain fixed at 51 and 75 min.
N_PERIOD_HALF = 45.0
N_PERIOD_HILL_COEFFICIENT = 1.1
N_AMPLITUDE_HALF = 11.3
N_AMPLITUDE_HILL_COEFFICIENT = 1.0

# Direct pituitary E2 feedback is a decreasing sensitivity Hill response.
# Its half-saturation is constrained to the 21-58-pg/mL interval bracketed by
# Shaw et al.'s measured E2 values, and its fitted exponent remains below 5.
# Sensitivity falls from 1 at E2=0 toward 1 - PITUITARY_E2_MAX_INHIBITION at
# saturating E2. The rounded 20% span is a
# conservative approximation to the 23-26% reduction in LH responses to
# controlled GnRH challenges reported by Shaw et al. after one month of
# estrogen. It is a model ceiling, not a measured concentration-response fit.
PITUITARY_E2_MAX_INHIBITION = 0.2
PITUITARY_E2_HALF_PG_ML = 21.0
PITUITARY_E2_HILL_COEFFICIENT = 4.9

# Maximum (zero-E2) GnRH-to-LH gain, calibrated against the digitized Licinio
# waveform after removing the endocrine circadian-amplitude envelope.
LH_GNRH_MAX_GAIN = 13.4

# Constant GnRH pulse scale. This is deliberately not circadian.
GNRH_PULSE_AMPLITUDE_SCALE = 102.4


def zeitgeber_hour(time_h):
    """Map simulation time to the shared repeating ZT0-ZT24 clock."""
    return np.mod(time_h, ZT_PERIOD_H)


def simulate_endocrine(
    output_h=OUTPUT_H,
    pituitary_e2_max_inhibition=PITUITARY_E2_MAX_INHIBITION,
):
    """Integrate N -> GnRH -> LH -> E2 with an explicit Euler loop."""
    total_h = WARMUP_H + output_h
    steps = round(total_h / STEP_H)
    time_h = np.arange(steps + 1) * STEP_H # step + 1 adds ic + steps
    n_state = np.empty(steps + 1)
    gnrh = np.empty(steps + 1)
    lh = np.empty(steps + 1)
    e2 = np.empty(steps + 1)
    feedback_half_hill = N_E2_HALF_PG_ML**N_E2_HILL_COEFFICIENT

    # N is expressed on a 0-100 scale.
    gnrh[0], lh[0], e2[0] = 0.0, 5.0, 15.0
    n_state[0] = N_MAX * feedback_half_hill / (
        feedback_half_hill + e2[0]**N_E2_HILL_COEFFICIENT
    )
    # phi is the accumulated pulse-clock cycle count. A pulse fires whenever
    # phi crosses the next integer, exactly as in the earlier model.
    phi = 0.0

    for i in range(steps):
        n_value = n_state[i]
        # Aggregate KNDy/kisspeptin drive acts directly on period in minutes.
        # Stronger N shortens the interval from 75 toward 51 minutes.
        period_min = (
            75.0
            - 24.0
            * n_value**N_PERIOD_HILL_COEFFICIENT
            / (
                N_PERIOD_HALF**N_PERIOD_HILL_COEFFICIENT
                + n_value**N_PERIOD_HILL_COEFFICIENT
            )
        )

        phi_new = phi + 60 * STEP_H / period_min # times 60 converts period to hrs
        fires = np.floor(phi_new) > np.floor(phi)
        pulse = 0.0
        if fires:
            # Neuropeptide Hill response with a constant (non-circadian) scale.
            pulse = (
                GNRH_PULSE_AMPLITUDE_SCALE
                * 3.5
                * n_value**N_AMPLITUDE_HILL_COEFFICIENT
                / (
                    N_AMPLITUDE_HALF**N_AMPLITUDE_HILL_COEFFICIENT
                    + n_value**N_AMPLITUDE_HILL_COEFFICIENT
                )
            )
        phi = phi_new

        e2_feedback = feedback_half_hill / (
            feedback_half_hill + e2[i]**N_E2_HILL_COEFFICIENT
        )

        # Simple Euler method. The neuropeptide production and decay terms are
        # written separately to keep the half-saturation and Hill coefficient
        # explicit.
        n_state[i + 1] = n_value + STEP_H * (
            N_TURNOVER_RATE_H
            * N_MAX
            * e2_feedback
            - N_TURNOVER_RATE_H * n_value
        )
        gnrh[i + 1] = (
            gnrh[i] - STEP_H * K_GNRH * gnrh[i] + pulse
        )
        pituitary_gain = 1.0 - pituitary_e2_max_inhibition * (
            e2[i] ** PITUITARY_E2_HILL_COEFFICIENT
            / (
                PITUITARY_E2_HALF_PG_ML
                ** PITUITARY_E2_HILL_COEFFICIENT
                + e2[i] ** PITUITARY_E2_HILL_COEFFICIENT
            )
        )
        lh[i + 1] = lh[i] + STEP_H * (
            LH_GNRH_MAX_GAIN * pituitary_gain * gnrh[i]
            - K_LH * lh[i]
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


# Temperature model. SCN-VP firing peaks near ZT7, but the AVP temperature
# output is hypothermic. This thermal timing remains literature-anchored; the
# endocrine pulse amplitude has no circadian term.
SCN_NEURAL_PEAK_ZT = 7.0

# Thermal measurements and literature-informed scales.
# Jin et al. report 22 +/- 1 C housing. Use a value inside that measured
# interval to align the absolute level with the digitized nonpregnant trace.
AMBIENT_TEMPERATURE_C = 21.25
INITIAL_TEMPERATURE_C = 37.0

# Zengin et al. Figure 1 short-term WT-OVX group. These group means were
# digitized from panels b and c because the paper does not tabulate them.
# Their limited precision is intentional: they are approximately 22.36 g
# total body mass and 3.98 g whole-body fat mass.
ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G = 22.36
ZENGIN_SHORT_TERM_WT_OVX_FAT_MASS_G = 3.98
ZENGIN_SHORT_TERM_WT_OVX_FAT_FRACTION = (
    ZENGIN_SHORT_TERM_WT_OVX_FAT_MASS_G
    / ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
)

# Abreu-Vieira et al. use cp = 0.9100 - 0.5231*f cal/(g degrees C).
MOUSE_SPECIFIC_HEAT_CAL_G_C = (
    0.9100 - 0.5231 * ZENGIN_SHORT_TERM_WT_OVX_FAT_FRACTION
)
BODY_HEAT_CAPACITY_KCAL_C = (
    MOUSE_SPECIFIC_HEAT_CAL_G_C
    * ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
    / 1000.0
)

# Lower-end calibration within the uncertainty of the approximately
# 0.32-kcal/h short-term WT-OVX bar in Zengin et al. Figure 5A. This is not
# asserted to be the reported group mean.
OVX_BACKGROUND_HEAT_KCAL_H = 0.31

# Backward-compatible alias for earlier Meeting 8 imports. New code should
# use OVX_BACKGROUND_HEAT_KCAL_H to avoid implying a directly tabulated mean.
OVX_MEAN_RESTING_HEAT_KCAL_H = OVX_BACKGROUND_HEAT_KCAL_H

# Zengin et al. did not measure whole-body conductance. Scale the
# Abreu-Vieira 26.97-g male estimate by the two-thirds power of body mass,
# an approximate body-surface-area transfer, using the same Zengin mass.
ABREU_VIEIRA_REFERENCE_BODY_MASS_G = 26.97
ABREU_VIEIRA_REFERENCE_CONDUCTANCE_KCAL_H_C = 0.0255
WHOLE_BODY_CONDUCTANCE_KCAL_H_C = (
    ABREU_VIEIRA_REFERENCE_CONDUCTANCE_KCAL_H_C
    * (
        ZENGIN_SHORT_TERM_WT_OVX_BODY_MASS_G
        / ABREU_VIEIRA_REFERENCE_BODY_MASS_G
    ) ** (2.0 / 3.0)
)

# E2 adds whole-body thermogenesis; it is not assigned exclusively to BAT.
# Upper end of the approximately 10-14% mouse whole-body heat difference
# supported by the Camporez et al. cross-check. Martinez de Morentin et al.'s
# rat central-E2 comparison gives an approximately 11% lower anchor.
E2_WHOLE_BODY_THERMOGENIC_FRACTION = 0.14

# Shape-calibrated thermal gate. K=11 pg/mL and n=5 are modeling choices, not
# an experimental EC50 or measured receptor cooperativity. K lies between the
# cycling trough and the 13.3-pg/mL aged input and near the low end of the
# 8.3-39.1-pg/mL serum-E2 bracket reported by Camporez et al.
E2_HEAT_HALF_PG_ML = 11.0
E2_HEAT_HILL_COEFFICIENT = 5.0

# Temperature-equivalent inputs. B is a 24-hour mean, so the SCN term is
# centered on zero rather than added as a one-directional heat input. The
# +/-0.19 C allocation is half of the 0.38 C activity-independent subjective
# day-night resting-temperature contrast reported in male rats by Scheer et
# al. (2005; doi:10.1016/j.neuroscience.2004.12.012). Scheer reported phase
# averages rather than extrema, so 0.19 C is an explicit modeling allocation.
SCN_TEMPERATURE_AMPLITUDE_C = 0.19
ACTIVITY_TEMPERATURE_DIFFERENCE_C = 1.9

# Piecewise-constant activity schedule: (start ZT, end ZT, relative height).
ACTIVITY_BOUTS_ZT = (
    (12.0, 17.0, 1.0),
)

# Comparison-condition inputs.
PREGNANT_E2_PG_ML = 50.88
AGED_CONSTANT_E2_PG_ML = 13.3
PREGNANT_ACTIVITY_MULTIPLIER = 0.7


# Derived ODE coefficients (degrees Celsius per hour where applicable).
THERMAL_RATE_CONSTANT_H = (
    WHOLE_BODY_CONDUCTANCE_KCAL_H_C / BODY_HEAT_CAPACITY_KCAL_C
)
BASAL_TEMPERATURE_RATE_C_H = (
    OVX_BACKGROUND_HEAT_KCAL_H / BODY_HEAT_CAPACITY_KCAL_C
)
E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H = (
    E2_WHOLE_BODY_THERMOGENIC_FRACTION
    * OVX_BACKGROUND_HEAT_KCAL_H
    / BODY_HEAT_CAPACITY_KCAL_C
)


def e2_heat_hill(e2_pg_ml):
    """Return the increasing K=11, n=5 whole-body E2 heat response."""
    numerator = e2_pg_ml**E2_HEAT_HILL_COEFFICIENT
    return numerator / (
        E2_HEAT_HALF_PG_ML**E2_HEAT_HILL_COEFFICIENT
        + numerator
    )


def scn_drive(zt_h):
    """Return the zero-mean SCN temperature drive in [-1, 1]."""
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
            in_bout = (phase_h >= start_h) | (phase_h < end_h)
        activity = np.where(in_bout, height, activity)

    # returns an array of [0,1,1...] at time points with activity
    return activity


def temperature_rate(
    temperature_c,
    e2_pg_ml,
    scn,
    activity,
    activity_multiplier=1.0,
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
    heat_loss = THERMAL_RATE_CONSTANT_H * (
        temperature_c - AMBIENT_TEMPERATURE_C
    )
    return (
        BASAL_TEMPERATURE_RATE_C_H
        + e2_whole_body_thermogenesis
        + scn_thermal_effect
        + activity_thermal_effect
        - heat_loss
    )


def simulate_temperature(
    time_h,
    e2,
    activity_multiplier=1.0,
):
    """Integrate the temperature ODE on the shared endocrine time grid."""
    temperature = np.empty(len(time_h))
    temperature[0] = INITIAL_TEMPERATURE_C
    zt_h = zeitgeber_hour(time_h) # array of time in ZT hour
    scn = scn_drive(zt_h) # array of scn values
    activity = activity_drive(zt_h) # get array of activity values

    # pull values out of the arrays to calculate temperature rate and integrate
    for i in range(len(time_h) - 1):
        temperature[i + 1] = temperature[i] + STEP_H * temperature_rate(
            temperature[i],
            e2[i],
            scn[i],
            activity[i],
            activity_multiplier,
        )
        
    # these return arrays of values
    e2_thermogenesis = (
        E2_WHOLE_BODY_THERMOGENIC_GAIN_C_H * e2_heat_hill(e2)
    )
    scn_effect = (
        THERMAL_RATE_CONSTANT_H * SCN_TEMPERATURE_AMPLITUDE_C * scn
    )
    activity_effect = (
        THERMAL_RATE_CONSTANT_H
        * activity_multiplier
        * ACTIVITY_TEMPERATURE_DIFFERENCE_C
        * activity
    )
    heat_loss = THERMAL_RATE_CONSTANT_H * (
        temperature - AMBIENT_TEMPERATURE_C
    )

    return {
        "temperature_c": temperature,
        "scn_drive": scn,
        "activity_drive": activity,
        "basal_rate_c_h": np.full(len(time_h), BASAL_TEMPERATURE_RATE_C_H),
        "e2_whole_body_thermogenesis_c_h": e2_thermogenesis,
        "scn_thermal_effect_c_h": scn_effect,
        "activity_thermal_effect_c_h": activity_effect,
        "heat_loss_c_h": heat_loss,
    }


def run():
    """Return 72 hours of endocrine output and four temperature conditions."""
    endocrine = simulate_endocrine()
    output_points = round(OUTPUT_H / STEP_H) + 1

    nonpregnant_thermal = simulate_temperature(
        endocrine["time_h"], endocrine["E2"]
    )

    constant_low_e2 = np.full_like(
        endocrine["E2"], AGED_CONSTANT_E2_PG_ML
    )
    constant_low_thermal = simulate_temperature(
        endocrine["time_h"], constant_low_e2
    )

    constant_high_e2 = np.full_like(
        endocrine["E2"], PREGNANT_E2_PG_ML
    )
    constant_high_thermal = simulate_temperature(
        endocrine["time_h"], constant_high_e2
    )
    pregnant_thermal = simulate_temperature(
        endocrine["time_h"],
        constant_high_e2,
        activity_multiplier=PREGNANT_ACTIVITY_MULTIPLIER,
    )

    result = {
        "hour": endocrine["time_h"][-output_points:],
        "zt_h": zeitgeber_hour(endocrine["time_h"][-output_points:]),
        "N": endocrine["N"][-output_points:],
        "GnRH": endocrine["GnRH"][-output_points:],
        "LH": endocrine["LH"][-output_points:],
        "E2_pg_ml": endocrine["E2"][-output_points:],
    }

    result.update({
        f"nonpregnant_{name}": values[-output_points:]
        for name, values in nonpregnant_thermal.items()
    })

    result["constant_low_E2_pg_ml"] = constant_low_e2[-output_points:]
    result.update({
        f"constant_low_{name}": values[-output_points:]
        for name, values in constant_low_thermal.items()
    })

    result["constant_high_E2_pg_ml"] = constant_high_e2[-output_points:]
    result.update({
        f"constant_high_{name}": values[-output_points:]
        for name, values in constant_high_thermal.items()
    })

    result["pregnant_E2_pg_ml"] = constant_high_e2[-output_points:]
    result.update({
        f"pregnant_{name}": values[-output_points:]
        for name, values in pregnant_thermal.items()
    })

    result["hour"] = result["hour"] - result["hour"][0]
    return result
