"""Plot the side-by-side reference-amplitude model without changing MODEL.py."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plot_model as plots
from MODEL_closer_to_manuscript import run


HERE = Path(__file__).resolve().parent


def plot_manuscript_waveform_comparison(
    series: dict,
    condition_name: str,
    reference_values: list[float],
    reference_label: str,
    model_label: str,
    title: str,
    output_path: Path,
) -> None:
    """Write one standalone waveform-aligned manuscript panel."""
    reference_time_h, _, _ = plots.read_mousebody_reference()
    reference_time_days = np.asarray(reference_time_h) / 24.0
    model = series[condition_name]
    aligned, shift_h = plots.phase_align_model(
        reference_time_h,
        reference_values,
        model["hour"],
        model["temperature_c"],
    )

    figure, axis = plt.subplots(figsize=(9.0, 4.8))
    axis.plot(
        reference_time_days,
        reference_values,
        color="#f912f2",
        linewidth=1.3,
        label=reference_label,
    )
    axis.plot(
        reference_time_days,
        aligned,
        color=plots.COLORS[condition_name],
        linewidth=1.5,
        label=f"{model_label} (phase aligned by {shift_h:.1f} h)",
    )
    axis.set(
        xlim=(0.0, 5.0),
        xlabel="Relative time (days)",
        ylabel="Body temperature (°C)",
        title=title,
    )
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def print_metrics(series: dict, clock_anchored_series: dict) -> None:
    """Print amplitudes and comparison RMSEs for quick verification."""
    for name in ("cycling", "constant_low", "pregnancy"):
        values = series[name]["temperature_c"]
        print(
            f"{name}: mean={np.mean(values):.3f} C, "
            f"amplitude={np.ptp(values):.3f} C, "
            f"range={np.min(values):.3f}-{np.max(values):.3f} C"
        )

    mousebody_time_h, nonpregnant, pregnant_control = (
        plots.read_mousebody_reference()
    )
    age_time_h, _, aged = plots.read_mouse_reference()
    model_time_h = series["cycling"]["hour"]
    cycling_aligned, _ = plots.phase_align_model(
        mousebody_time_h,
        nonpregnant,
        model_time_h,
        series["cycling"]["temperature_c"],
    )
    pregnancy_aligned, _ = plots.phase_align_model(
        mousebody_time_h,
        pregnant_control,
        model_time_h,
        series["pregnancy"]["temperature_c"],
    )
    aged_anchored = np.interp(
        age_time_h,
        clock_anchored_series["constant_low"]["hour"],
        clock_anchored_series["constant_low"]["temperature_c"],
    )
    print(
        "comparison RMSE: "
        f"cycling={np.sqrt(np.mean((cycling_aligned - nonpregnant) ** 2)):.3f} C, "
        f"constant_low={np.sqrt(np.mean((aged_anchored - aged) ** 2)):.3f} C, "
        f"pregnancy={np.sqrt(np.mean((pregnancy_aligned - pregnant_control) ** 2)):.3f} C"
    )


def main() -> None:
    series = run()
    clock_anchored_series = run(output_start_clock_h=6.0)
    _, nonpregnant, pregnant_control = plots.read_mousebody_reference()
    outputs = {
        "timeseries": HERE / "reference_amplitude_model_timeseries.png",
        "comparison": (
            HERE / "reference_amplitude_temperature_comparison.png"
        ),
        "cycling_2011": (
            HERE / "reference_amplitude_cycling_2011_comparison.png"
        ),
        "cycling_2011_waveform_aligned": (
            HERE
            / "reference_amplitude_cycling_2011_waveform_aligned.png"
        ),
        "constant_low_2011_waveform_aligned": (
            HERE
            / "reference_amplitude_constant_low_2011_waveform_aligned.png"
        ),
        "cycling_mousebody": (
            HERE / "reference_amplitude_cycling_mousebody_comparison.png"
        ),
        "manuscript_cycling": (
            HERE / "manuscript_cycling_e2_waveform_aligned.png"
        ),
        "manuscript_pregnancy": (
            HERE / "manuscript_pregnancy_high_e2_waveform_aligned.png"
        ),
    }
    plots.plot_timeseries(series, outputs["timeseries"])
    plots.plot_reference_comparison(
        series,
        clock_anchored_series,
        outputs["comparison"],
    )
    plots.plot_cycling_2011_comparison(
        clock_anchored_series,
        outputs["cycling_2011"],
    )
    plots.plot_cycling_2011_waveform_aligned(
        series,
        outputs["cycling_2011_waveform_aligned"],
    )
    plots.plot_constant_low_2011_waveform_aligned(
        series,
        outputs["constant_low_2011_waveform_aligned"],
    )
    plots.plot_manuscript_cycling_comparison(
        series,
        outputs["cycling_mousebody"],
    )
    plot_manuscript_waveform_comparison(
        series,
        "cycling",
        nonpregnant,
        "Non-pregnant mouse",
        "Cycling E2 model",
        "Cycling E2 model vs non-pregnant mouse — waveform aligned",
        outputs["manuscript_cycling"],
    )
    plot_manuscript_waveform_comparison(
        series,
        "pregnancy",
        pregnant_control,
        "Pregnant-control mouse",
        "Pregnancy/high E2 model",
        "Pregnancy/high E2 model vs pregnant-control mouse — waveform aligned",
        outputs["manuscript_pregnancy"],
    )
    print_metrics(series, clock_anchored_series)
    for path in outputs.values():
        print(f"Wrote: {path}")


if __name__ == "__main__":
    main()
