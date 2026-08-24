"""Plot the final dynamic time series from dynamic_temp_model.py.

Run:
    python plot_temp_timeseries.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from dynamic_temp_model import calculate_dynamic_timeseries


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "final_model_timeseries.png"


def plot_timeseries(output_path: Path = DEFAULT_OUTPUT) -> None:
    """Calculate the three model traces and write their time-series plot."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Plotting requires matplotlib") from error

    cycling, menopause, pregnancy = calculate_dynamic_timeseries()
    days = [row["hour"] / 24.0 for row in cycling]

    figure, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(
        days,
        [row["temperature_c"] for row in cycling],
        color="#5146c7",
        linewidth=1.5,
        label="Cycling E2",
    )
    axis.plot(
        days,
        [row["temperature_c"] for row in menopause],
        color="#14833b",
        linewidth=1.5,
        linestyle="--",
        label=(
            f"Constant low E2: {menopause[0]['e2_pg_ml']:.3f} pg/mL, "
            "P 0.2 ng/mL"
        ),
    )
    axis.plot(
        days,
        [row["temperature_c"] for row in pregnancy],
        color="#c2410c",
        linewidth=1.5,
        linestyle=":",
        label="Pregnancy: E2 215 pg/mL, P 25 ng/mL",
    )
    axis.set_xlabel("Time (days)")
    axis.set_ylabel("Body temperature (°C)")
    axis.set_title("Final dynamic temperature model")
    axis.legend(frameon=False)
    axis.grid(alpha=0.22)
    figure.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="PNG output path (default: final_model_timeseries.png beside this file).",
    )
    args = parser.parse_args()
    plot_timeseries(args.output)
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
