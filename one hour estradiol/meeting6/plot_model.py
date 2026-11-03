"""Plot the final model outputs; MODEL.py contains no plotting code.

Run ``python plot_model.py`` to regenerate all final PNG files.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from MODEL import run


HERE = Path(__file__).resolve().parent
AGE_REFERENCE_DATA = HERE / "references" / "mice_age_temperature_digitized.csv"
MOUSEBODY_REFERENCE_DATA = (
    HERE / "references" / "mousebodytemp_digitized.csv"
)
COLORS = {
    "cycling": "#5046e5",
    "constant_low": "#16803b",
    "pregnancy": "#c2410c",
}


def read_mouse_reference() -> tuple[list[float], list[float], list[float]]:
    """Read the digitized 3-month and 24-month female reference traces."""
    time_h: list[float] = []
    young: list[float] = []
    aged: list[float] = []
    with AGE_REFERENCE_DATA.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            time_h.append(float(row["hour"]))
            young.append(float(row["three_month_f_c"]))
            aged.append(float(row["twenty_four_month_f_c"]))
    return time_h, young, aged


def read_mousebody_reference() -> tuple[list[float], list[float], list[float]]:
    """Read non-pregnant and pregnant-control traces from mousebodytemp.png."""
    time_h: list[float] = []
    nonpregnant: list[float] = []
    pregnant_control: list[float] = []
    with MOUSEBODY_REFERENCE_DATA.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            time_h.append(24.0 * float(row["relative_day"]))
            nonpregnant.append(float(row["nonpregnant_c"]))
            pregnant_control.append(float(row["pregnant_control_c"]))
    return time_h, nonpregnant, pregnant_control


def phase_align_model(
    reference_time_h: list[float],
    reference_values: list[float],
    model_time_h: np.ndarray,
    model_values: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Circularly align a model trace to a reference by waveform correlation.

    The search is restricted to one 24-hour cycle.  Both traces are standardized
    before correlation, so the alignment uses waveform timing rather than their
    different means or amplitudes.  This is for a shape-comparison figure; it is
    deliberately not presented as an inferred clock-time alignment.
    """
    reference_time = np.asarray(reference_time_h, dtype=float)
    reference = np.asarray(reference_values, dtype=float)
    model_time = np.asarray(model_time_h, dtype=float)
    model = np.asarray(model_values, dtype=float)
    reference_standard = (reference - np.mean(reference)) / np.std(reference)
    period_h = float(model_time[-1] - model_time[0])

    best_score = -np.inf
    best_shift_h = 0.0
    best_values = np.interp(reference_time, model_time, model)
    for shift_h in np.arange(0.0, 24.0, 5.0 / 60.0):
        sample_time_h = (reference_time + shift_h) % period_h
        candidate = np.interp(sample_time_h, model_time, model)
        candidate_standard = (candidate - np.mean(candidate)) / np.std(candidate)
        score = float(np.mean(reference_standard * candidate_standard))
        if score > best_score:
            best_score = score
            best_shift_h = float(shift_h)
            best_values = candidate
    return best_values, best_shift_h


def plot_timeseries(series: dict, output_path: Path) -> None:
    """Plot the three final 120-hour temperature conditions."""
    time_h = series["cycling"]["hour"]
    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    for name, label, style in [
        ("cycling", "Cycling E2", "-"),
        ("constant_low", "Constant low E2", "--"),
        ("pregnancy", "Pregnancy", ":"),
    ]:
        values = series[name]["temperature_c"]
        axis.plot(
            time_h,
            values,
            style,
            color=COLORS[name],
            linewidth=1.5,
            label=f"{label} (trough {min(values):.2f} °C)",
        )
    axis.set(
        xlabel="Time (h)",
        ylabel="Body temperature (°C)",
        title="Final dynamic temperature model — Euler integration",
    )
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=3)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_e2_timeseries(series: dict, output_path: Path) -> None:
    """Plot one modeled day of cycling estradiol, from 08:00 to 08:00."""
    time_h = series["cycling"]["hour"]
    e2_pg_ml = series["cycling"]["E2_pg_ml"]
    one_day = (time_h >= 8.0) & (time_h <= 32.0)
    figure, axis = plt.subplots(figsize=(9.0, 3.6))
    axis.plot(
        time_h[one_day],
        e2_pg_ml[one_day],
        color="#FC05E3",
        linewidth=1.5,
    )
    axis.set(
        xlabel="Time of day",
        ylabel="Estradiol, E2 (pg/mL)",
        xlim=(8.0, 32.0),
        ylim=(0.0, 55.0),
    )
    axis.set_xticks(
        [8.0, 14.0, 20.0, 26.0, 32.0],
        ["08:00", "14:00", "20:00", "02:00", "08:00"],
    )
    axis.set_yticks([0.0, 20.0, 40.0])
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_reference_comparison(
    series: dict,
    clock_anchored_series: dict,
    output_path: Path,
) -> None:
    """Plot two waveform alignments and one clock-anchored comparison."""
    age_time_h, _, aged = read_mouse_reference()
    mousebody_time_h, nonpregnant, pregnant_control = (
        read_mousebody_reference()
    )
    model_time_h = series["cycling"]["hour"]
    mousebody_time_days = np.asarray(mousebody_time_h) / 24.0
    age_time_days = np.asarray(age_time_h) / 24.0
    figure, axes = plt.subplots(3, 1, figsize=(8.0, 12.0), sharex=True)
    figure.suptitle(
        "Temperature-model comparisons: phase aligned or clock anchored",
        fontsize=15,
    )

    cycling_aligned, cycling_shift_h = phase_align_model(
        mousebody_time_h,
        nonpregnant,
        model_time_h,
        series["cycling"]["temperature_c"],
    )
    axes[0].plot(
        mousebody_time_days,
        nonpregnant,
        color="#f912f2",
        linewidth=1.3,
        label="Non-pregnant mousebodytemp",
    )
    axes[0].plot(
        mousebody_time_days,
        cycling_aligned,
        color=COLORS["cycling"],
        linewidth=1.5,
        label=f"Cycling E2 model (phase shift {cycling_shift_h:.1f} h)",
    )
    axes[0].set_title("Cycling E2 — waveform aligned")

    aged_clock_anchored = np.interp(
        age_time_h,
        clock_anchored_series["constant_low"]["hour"],
        clock_anchored_series["constant_low"]["temperature_c"],
    )
    for day in range(5):
        axes[1].axvspan(
            day + 0.5,
            day + 1.0,
            color="#777777",
            alpha=0.10,
            linewidth=0,
        )
    axes[1].plot(
        age_time_days,
        aged,
        color="#f912f2",
        linewidth=1.3,
        label="24-month F mouse (2011)",
    )
    axes[1].plot(
        age_time_days,
        aged_clock_anchored,
        color=COLORS["constant_low"],
        linewidth=1.5,
        label="Constant low E2 model (clock anchored)",
    )
    axes[1].set_title(
        "Constant low E2 — 2011 clock anchor: day starts 06:00 (ZT0)"
    )

    pregnancy_aligned, pregnancy_shift_h = phase_align_model(
        mousebody_time_h,
        pregnant_control,
        model_time_h,
        series["pregnancy"]["temperature_c"],
    )
    axes[2].plot(
        mousebody_time_days,
        pregnant_control,
        color="#f912f2",
        linewidth=1.3,
        label="Pregnant-control mousebodytemp",
    )
    axes[2].plot(
        mousebody_time_days,
        pregnancy_aligned,
        color=COLORS["pregnancy"],
        linewidth=1.5,
        label=f"Pregnancy/high E2 model (phase shift {pregnancy_shift_h:.1f} h)",
    )
    axes[2].set_title("Pregnancy/high E2 — waveform aligned")

    for axis in axes:
        axis.set_ylabel("Body temperature (°C)")
        axis.set_xlim(0.0, 5.0)
        axis.set_xticks(np.arange(0.0, 5.1, 1.0))
        axis.grid(alpha=0.22)
        axis.legend(frameon=False, ncol=2)
    axes[2].set_xlabel("Relative time (days)")
    figure.tight_layout(rect=(0, 0, 1, 0.975))
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_cycling_2011_comparison(
    clock_anchored_series: dict,
    output_path: Path,
) -> None:
    """Compare cycling E2 with the young 2011 trace from 06:00/ZT0."""
    reference_time_h, young, _ = read_mouse_reference()
    reference_time_days = np.asarray(reference_time_h) / 24.0
    cycling = clock_anchored_series["cycling"]
    cycling_temperature = np.interp(
        reference_time_h,
        cycling["hour"],
        cycling["temperature_c"],
    )

    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    for day in range(5):
        axis.axvspan(
            day + 0.5,
            day + 1.0,
            color="#777777",
            alpha=0.10,
            linewidth=0,
        )
    axis.plot(
        reference_time_days,
        young,
        color="#f912f2",
        linewidth=1.3,
        label="3-month F mouse (2011)",
    )
    axis.plot(
        reference_time_days,
        cycling_temperature,
        color=COLORS["cycling"],
        linewidth=1.5,
        label="Cycling E2 model (clock anchored)",
    )
    axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title="Cycling E2 vs 2011 young female — day starts 06:00 (ZT0)",
    )
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_manuscript_cycling_comparison(
    series: dict,
    output_path: Path,
) -> None:
    """Plot the waveform-aligned cycling-E2/mousebodytemp comparison."""
    reference_time_h, nonpregnant, _ = read_mousebody_reference()
    reference_time_days = np.asarray(reference_time_h) / 24.0
    cycling_aligned, shift_h = phase_align_model(
        reference_time_h,
        nonpregnant,
        series["cycling"]["hour"],
        series["cycling"]["temperature_c"],
    )
    circular_offset_h = min(shift_h, 24.0 - shift_h)

    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    axis.plot(
        reference_time_days,
        nonpregnant,
        color="#f912f2",
        linewidth=1.3,
        label="Non-pregnant mouse",
    )
    axis.plot(
        reference_time_days,
        cycling_aligned,
        color=COLORS["cycling"],
        linewidth=1.5,
        label=(
            "Cycling E2 model "
            f"(phase aligned; {circular_offset_h:.1f} h offset)"
        ),
    )
    axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title=(
            "Cycling E2 model vs non-pregnant mouse temperature — "
            "waveform aligned"
        ),
    )
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_cycling_2011_waveform_aligned(
    series: dict,
    output_path: Path,
) -> None:
    """Plot cycling E2 and the young 2011 trace after waveform alignment."""
    reference_time_h, young, _ = read_mouse_reference()
    reference_time_days = np.asarray(reference_time_h) / 24.0
    cycling_aligned, shift_h = phase_align_model(
        reference_time_h,
        young,
        series["cycling"]["hour"],
        series["cycling"]["temperature_c"],
    )
    circular_offset_h = min(shift_h, 24.0 - shift_h)

    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    axis.plot(
        reference_time_days,
        young,
        color="#f912f2",
        linewidth=1.3,
        label="3-month F mouse (2011)",
    )
    axis.plot(
        reference_time_days,
        cycling_aligned,
        color=COLORS["cycling"],
        linewidth=1.5,
        label=(
            "Cycling E2 model "
            f"(phase aligned; {circular_offset_h:.1f} h offset)"
        ),
    )
    axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title="Cycling E2 vs 2011 young female — waveform aligned",
    )
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_constant_low_2011_waveform_aligned(
    series: dict,
    output_path: Path,
) -> None:
    """Plot constant-low E2 and the aged 2011 trace after alignment."""
    reference_time_h, _, aged = read_mouse_reference()
    reference_time_days = np.asarray(reference_time_h) / 24.0
    low_aligned, shift_h = phase_align_model(
        reference_time_h,
        aged,
        series["constant_low"]["hour"],
        series["constant_low"]["temperature_c"],
    )
    circular_offset_h = min(shift_h, 24.0 - shift_h)

    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    axis.plot(
        reference_time_days,
        aged,
        color="#f912f2",
        linewidth=1.3,
        label="24-month F mouse (2011)",
    )
    axis.plot(
        reference_time_days,
        low_aligned,
        color=COLORS["constant_low"],
        linewidth=1.5,
        label=(
            "Constant low E2 model "
            f"(phase aligned; {circular_offset_h:.1f} h offset)"
        ),
    )
    axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title="Constant low E2 vs 2011 aged female — waveform aligned",
    )
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def plot_model_amplitude_comparison(
    series: dict,
    output_path: Path,
) -> None:
    """Compare modeled temperature amplitudes across the three E2 states."""
    time_days = series["cycling"]["hour"] / 24.0
    figure, trace_axis = plt.subplots(figsize=(9.0, 4.8))
    conditions = [
        ("constant_low", "Constant low E2", "--"),
        ("cycling", "Cycling E2", "-"),
        (
            "pregnancy",
            "Constant high E2 (pregnancy hormones)",
            ":",
        ),
    ]
    for name, label, linestyle in conditions:
        values = series[name]["temperature_c"]
        trough = float(np.min(values))
        peak = float(np.max(values))
        amplitude = peak - trough
        trace_axis.plot(
            time_days,
            values,
            color=COLORS[name],
            linewidth=1.5,
            linestyle=linestyle,
            label=f"{label} (ΔT {amplitude:.2f} °C)",
        )

    trace_axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title="Modeled temperature across E2 conditions",
    )
    trace_axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    trace_axis.grid(alpha=0.22)
    trace_axis.legend(frameon=False, ncol=1)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    series = run()
    clock_anchored_series = run(output_start_clock_h=6.0)
    timeseries_path = args.output_dir / "final_model_timeseries.png"
    e2_timeseries_path = args.output_dir / "final_e2_timeseries.png"
    comparison_path = (
        args.output_dir / "final_temperature_reference_comparison.png"
    )
    cycling_2011_path = (
        args.output_dir / "cycling_e2_2011_clock_anchored_comparison.png"
    )
    manuscript_cycling_path = (
        args.output_dir / "manuscript cycling e2 comparison.png"
    )
    cycling_2011_aligned_path = (
        args.output_dir / "cycling e2 2011 waveform aligned.png"
    )
    constant_low_2011_aligned_path = (
        args.output_dir / "constant low e2 2011 waveform aligned.png"
    )
    amplitude_comparison_path = (
        args.output_dir / "model e2 condition amplitude comparison.png"
    )
    plot_timeseries(series, timeseries_path)
    plot_e2_timeseries(series, e2_timeseries_path)
    plot_reference_comparison(series, clock_anchored_series, comparison_path)
    plot_cycling_2011_comparison(clock_anchored_series, cycling_2011_path)
    plot_manuscript_cycling_comparison(series, manuscript_cycling_path)
    plot_cycling_2011_waveform_aligned(series, cycling_2011_aligned_path)
    plot_constant_low_2011_waveform_aligned(
        series,
        constant_low_2011_aligned_path,
    )
    plot_model_amplitude_comparison(series, amplitude_comparison_path)
    print(f"Wrote: {timeseries_path}")
    print(f"Wrote: {e2_timeseries_path}")
    print(f"Wrote: {comparison_path}")
    print(f"Wrote: {cycling_2011_path}")
    print(f"Wrote: {manuscript_cycling_path}")
    print(f"Wrote: {cycling_2011_aligned_path}")
    print(f"Wrote: {constant_low_2011_aligned_path}")
    print(f"Wrote: {amplitude_comparison_path}")


if __name__ == "__main__":
    main()
