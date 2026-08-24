"""Dynamic Meeting 6 body-temperature model with self-regulated heat loss.

The endocrine trace comes from the canonical ``MODEL.py``. Progesterone and
E2 enter the outer conservation and hormonal-loss pathways through one linear
coefficient each; neither hormone has its own Hill transform:

    dT/dt = b + a_C Conservation - a_L HormonalLoss
            - a_T T^n / (K_T^n + T^n)

The parameters below are the rounded, age-trace fit from
``fit_age_temperature.js``.  Run this file to export the fitted five-day
traces and their summary metrics.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from types import ModuleType

HERE = Path(__file__).resolve().parent
ENDOCRINE_MODEL_PATH = HERE.parent / "MODEL.py"

MODEL_STEP_MINUTES = 0.1
MODEL_TOTAL_DAYS = 9
SAMPLE_MINUTES = 5
ANALYSIS_DAYS = 5
CONSERVATION_CLOCK_WEIGHT = 2.27
N_E2_HALF_INHIBITION_PG_ML = 20.1
N_E2_HILL_EXPONENT = 5
CYCLING_PROGESTERONE_NG_ML = 0.2
MENOPAUSE_PROGESTERONE_NG_ML = 0.2
PREGNANCY_E2_PG_ML = 215.0
PREGNANCY_PROGESTERONE_NG_ML = 25.0

# Rounded age-trace fit.  Only simple pathway coefficients and the thermal
# clock center changed; the self-loss Hill function is unchanged.  Every rate
# has units degC/hour after division by whole-body heat capacity.
BASAL_HEAT_RATE_C_PER_H = 49.44
CONSERVATION_GAIN_C_PER_H = 9.76
HORMONAL_LOSS_GAIN_C_PER_H = 15.29
SELF_LOSS_MAX_C_PER_H = 125.0
SELF_LOSS_HALF_SATURATION_C = 39.08
SELF_LOSS_HILL_EXPONENT = 5
INITIAL_TEMPERATURE_C = 37.0
HORMONAL_LOSS_CLOCK_WEIGHT = 0.78
THERMAL_CLOCK_PEAK_HOURS = 18.16
FITTED_CONSERVATION_N_WEIGHT = 0.60
PROGESTERONE_COEFFICIENT_PER_NG_ML = 0.0274
E2_COEFFICIENT_PER_PG_ML = 0.0000869


def load_endocrine_model() -> ModuleType:
    """Load the canonical endocrine model from the parent directory."""
    spec = importlib.util.spec_from_file_location(
        "final_endocrine_model", ENDOCRINE_MODEL_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {ENDOCRINE_MODEL_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def thermal_hill(value: float) -> float:
    """Return the outer fourth-order pathway readout."""
    nonnegative = max(value, 0.0)
    power = nonnegative**4
    return power / (1.0 + power)


def temperature_metrics(values: list[float]) -> dict[str, float]:
    """Return the peak, trough, and range of a temperature trace."""
    peak = max(values)
    trough = min(values)
    return {"peak_c": peak, "trough_c": trough, "range_c": peak - trough}


def thermal_clock(hour: float) -> float:
    """Return the fitted unit-range thermal clock."""
    hour_of_day = hour % 24.0
    cosine = math.cos(
        2.0 * math.pi * (hour_of_day - THERMAL_CLOCK_PEAK_HOURS) / 24.0
    )
    return 0.5 * (1.0 + cosine)


def self_regulated_loss(temperature_c: float) -> float:
    """Return the direct temperature-driven Hill activation (0 to 1)."""
    ratio = max(temperature_c, 0.0) / SELF_LOSS_HALF_SATURATION_C
    power = ratio**SELF_LOSS_HILL_EXPONENT
    return power / (1.0 + power)


def temperature_rate(
    temperature_c: float,
    conservation: float,
    hormonal_loss: float,
) -> float:
    """Calculate dT/dt in degC/hour."""
    return (
        BASAL_HEAT_RATE_C_PER_H
        + CONSERVATION_GAIN_C_PER_H * conservation
        - HORMONAL_LOSS_GAIN_C_PER_H * hormonal_loss
        - SELF_LOSS_MAX_C_PER_H * self_regulated_loss(temperature_c)
    )


def integrate_temperature(rows: list[dict[str, float]]) -> list[float]:
    """Integrate the temperature ODE with fourth-order Runge-Kutta."""
    if not rows:
        return []
    step_hours = SAMPLE_MINUTES / 60.0
    values = [INITIAL_TEMPERATURE_C]
    temperature = INITIAL_TEMPERATURE_C
    for current, following in zip(rows, rows[1:]):
        c0 = current["conservation"]
        c1 = following["conservation"]
        l0 = current["loss"]
        l1 = following["loss"]
        cm = 0.5 * (c0 + c1)
        lm = 0.5 * (l0 + l1)

        k1 = temperature_rate(temperature, c0, l0)
        k2 = temperature_rate(temperature + 0.5 * step_hours * k1, cm, lm)
        k3 = temperature_rate(temperature + 0.5 * step_hours * k2, cm, lm)
        k4 = temperature_rate(temperature + step_hours * k3, c1, l1)
        temperature += step_hours * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        values.append(temperature)
    return values


def _full_endocrine_trace() -> tuple[list[float], list[float]]:
    """Return the full nine-day N and E2 traces at five-minute spacing."""
    model = load_endocrine_model()
    time, n_state, _gnrh, _lh, e2 = model.run(
        t_end=MODEL_TOTAL_DAYS * 1440,
        dt=MODEL_STEP_MINUTES,
        circ=True,
        return_n=True,
    )
    stride = round(SAMPLE_MINUTES / MODEL_STEP_MINUTES)
    return n_state[::stride].tolist(), e2[::stride].tolist()


def _condition_rows(
    n_state: list[float],
    e2: list[float],
    progesterone: float,
) -> list[dict[str, float]]:
    """Build upstream signals, integrate T, and retain the final five days."""
    rows: list[dict[str, float]] = []
    for index, (n_value, e2_value) in enumerate(zip(n_state, e2)):
        hour = index * SAMPLE_MINUTES / 60.0
        clock = thermal_clock(hour)
        conservation = thermal_hill(
            CONSERVATION_CLOCK_WEIGHT * clock
            + FITTED_CONSERVATION_N_WEIGHT * n_value
            + PROGESTERONE_COEFFICIENT_PER_NG_ML * progesterone
        )
        hormonal_loss = thermal_hill(
            HORMONAL_LOSS_CLOCK_WEIGHT * clock
            + E2_COEFFICIENT_PER_PG_ML * e2_value
        )
        rows.append(
            {
                "hour": hour,
                "N": n_value,
                "e2_pg_ml": e2_value,
                "progesterone_ng_ml": progesterone,
                "clock": clock,
                "conservation": conservation,
                "loss": hormonal_loss,
            }
        )
    dynamic_temperature = integrate_temperature(rows)
    for row, value in zip(rows, dynamic_temperature):
        row["temperature_c"] = value
        row["self_regulated_loss"] = self_regulated_loss(value)

    retained_intervals = round(ANALYSIS_DAYS * 24 * 60 / SAMPLE_MINUTES)
    retained = rows[-(retained_intervals + 1):]
    first_hour = retained[0]["hour"]
    for row in retained:
        row["hour"] -= first_hour
    return retained


def calculate_dynamic_timeseries() -> tuple[
    list[dict[str, float]],
    list[dict[str, float]],
    list[dict[str, float]],
]:
    """Calculate cycling, constant-low-E2, and pregnancy traces."""
    cycling_n, cycling_e2 = _full_endocrine_trace()
    count = len(cycling_e2)
    retained_intervals = round(ANALYSIS_DAYS * 24 * 60 / SAMPLE_MINUTES)
    constant_low_e2 = min(cycling_e2[-(retained_intervals + 1):])

    low_n_value = 1.0 / (
        1.0
        + (constant_low_e2 / N_E2_HALF_INHIBITION_PG_ML)
        ** N_E2_HILL_EXPONENT
    )
    pregnancy_n_value = 1.0 / (
        1.0
        + (PREGNANCY_E2_PG_ML / N_E2_HALF_INHIBITION_PG_ML)
        ** N_E2_HILL_EXPONENT
    )

    cycling = _condition_rows(
        cycling_n,
        cycling_e2,
        CYCLING_PROGESTERONE_NG_ML,
    )
    low = _condition_rows(
        [low_n_value] * count,
        [constant_low_e2] * count,
        MENOPAUSE_PROGESTERONE_NG_ML,
    )
    pregnancy = _condition_rows(
        [pregnancy_n_value] * count,
        [PREGNANCY_E2_PG_ML] * count,
        PREGNANCY_PROGESTERONE_NG_ML,
    )
    return cycling, low, pregnancy


def _summary(
    cycling: list[dict[str, float]],
    low: list[dict[str, float]],
    pregnancy: list[dict[str, float]],
) -> dict[str, object]:
    def metrics(rows: list[dict[str, float]]) -> dict[str, object]:
        return temperature_metrics([row["temperature_c"] for row in rows])

    return {
        "equation": (
            "dT/dt = b + aC*Conservation - aL*HormonalLoss "
            "- aT*T^n/(KT^n + T^n)"
        ),
        "parameters": {
            "b_c_per_h": BASAL_HEAT_RATE_C_PER_H,
            "aC_c_per_h": CONSERVATION_GAIN_C_PER_H,
            "aL_c_per_h": HORMONAL_LOSS_GAIN_C_PER_H,
            "aT_c_per_h": SELF_LOSS_MAX_C_PER_H,
            "KT_c": SELF_LOSS_HALF_SATURATION_C,
            "n": SELF_LOSS_HILL_EXPONENT,
            "initial_temperature_c": INITIAL_TEMPERATURE_C,
            "thermal_clock_peak_hours": THERMAL_CLOCK_PEAK_HOURS,
            "conservation_n_weight": FITTED_CONSERVATION_N_WEIGHT,
            "progesterone_coefficient_per_ng_ml": (
                PROGESTERONE_COEFFICIENT_PER_NG_ML
            ),
            "e2_coefficient_per_pg_ml": E2_COEFFICIENT_PER_PG_ML,
            "hormonal_loss_clock_weight": HORMONAL_LOSS_CLOCK_WEIGHT,
            "constant_low_e2_pg_ml": low[0]["e2_pg_ml"],
        },
        "temperature": {
            "cycling_e2": metrics(cycling),
            "constant_low_e2": metrics(low),
            "pregnancy": metrics(pregnancy),
        },
    }


def _write_csv(
    cycling: list[dict[str, float]],
    low: list[dict[str, float]],
    pregnancy: list[dict[str, float]],
    path: Path,
) -> None:
    """Write the three final dynamic temperature traces to one compact CSV."""
    fieldnames = [
        "hour",
        "cycling_temperature_c",
        "constant_low_temperature_c",
        "pregnancy_temperature_c",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for cycling_row, low_row, pregnancy_row in zip(cycling, low, pregnancy):
            writer.writerow(
                {
                    "hour": cycling_row["hour"],
                    "cycling_temperature_c": cycling_row["temperature_c"],
                    "constant_low_temperature_c": low_row["temperature_c"],
                    "pregnancy_temperature_c": pregnancy_row["temperature_c"],
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cycling, low, pregnancy = calculate_dynamic_timeseries()
    summary = _summary(cycling, low, pregnancy)
    json_path = args.output_dir / "dynamic_temp_model_results.json"
    csv_path = args.output_dir / "dynamic_temp_model_timeseries.csv"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_csv(cycling, low, pregnancy, csv_path)
    print(json.dumps(summary, indent=2))
    print(f"\nWrote: {json_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
