"""Render the Meeting 8 temperature, estradiol, and Licinio comparisons.

All curves come directly from ``MODEL.py``. Running this file writes three PNG
figures beside it: the four-condition temperature and estradiol time series
and the three-day Licinio comparison folded onto one 24-hour axis.
"""

from pathlib import Path

import matplotlib

# Agg renders figures directly to files instead of opening pop-up windows.
# It must be selected before importing pyplot, which initializes the backend.
matplotlib.use("Agg")

# mpimg reads the Licinio source PNG; pyplot creates figures and axes.
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from MODEL import (
    AMBIENT_TEMPERATURE_C,
    E2_HEAT_HALF_PG_ML,
    E2_HEAT_HILL_COEFFICIENT,
    E2_WHOLE_BODY_THERMOGENIC_FRACTION,
    OVX_BACKGROUND_HEAT_KCAL_H,
    PITUITARY_E2_HALF_PG_ML,
    PITUITARY_E2_HILL_COEFFICIENT,
    WARMUP_H,
    run,
    simulate_endocrine,
)


# Keep all fixed inputs and outputs relative to this script, so running
# ``python PLOT.py`` works regardless of the current working directory.
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "temperature_timeseries_four_conditions.png"
DEFAULT_ESTROGEN_OUTPUT = HERE / "estradiol_timeseries_four_conditions.png"
DEFAULT_LICINIO_OUTPUT = HERE / "licinio_three_day_folded_comparison.png"
DEFAULT_LICINIO_TARGET = HERE.parent / "meeting6" / "[1] target.png"

# Each curve description is:
# (key in MODEL.run() result, legend label, hexadecimal color, line style).
SERIES = (
    ("constant_low_temperature_c", "Constant aged E2", "#16803B", "--"),
    ("constant_high_temperature_c", "Constant high E2", "#7C3AED", "-"),
    ("nonpregnant_temperature_c", "Cycling E2 / nonpregnant", "#2563EB", "-"),
    (
        "pregnant_temperature_c",
        "Constant high E2 + reduced activity",
        "#C2410C",
        "-",
    ),
)

# Four temperature conditions use only three distinct E2 inputs because both
# constant-high-E2 temperature conditions share the same hormone input.
ESTROGEN_SERIES = (
    ("constant_low_E2_pg_ml", "Constant aged E2", "#16803B", "--"),
    (
        "constant_high_E2_pg_ml",
        "Constant high E2 (both activity levels)",
        "#7C3AED",
        "-",
    ),
    ("E2_pg_ml", "Cycling E2 / nonpregnant", "#2563EB", "-"),
)

DAY_COLORS = ("#2563EB", "#DD7023", "#844CC4")
TEXT_COLOR = "#18202F"
MUTED_COLOR = "#5B6471"
GRID_COLOR = "#D8DDE4"
AXIS_COLOR = "#4A525E"


def configure_style():
    """Set the shared Matplotlib styling used by every figure."""
    # rcParams are global defaults applied to every figure created afterward.
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.edgecolor": AXIS_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "axes.linewidth": 1.0,
        "xtick.color": AXIS_COLOR,
        "ytick.color": AXIS_COLOR,
        "legend.labelcolor": TEXT_COLOR,
    })


def style_axes(ax):
    # An Axes is one rectangular plotting panel inside the complete Figure.
    # Keep every panel's grid and tick labels visually consistent.
    ax.grid(True, color=GRID_COLOR, linewidth=0.65, alpha=0.9)
    ax.tick_params(labelsize=9)
    # Draw the grid behind the data curves instead of over them.
    ax.set_axisbelow(True)


def add_heading(fig, title, subtitle):
    # Figure-level text spans the entire image rather than one Axes panel.
    # Figure coordinates run from 0 at the left/bottom to 1 at the right/top.
    fig.suptitle(
        title,
        x=0.08,
        y=0.975,
        ha="left",
        fontsize=15,
        fontweight="bold",
        color=TEXT_COLOR,
    )
    # Put the smaller subtitle just below the main title.
    fig.text(
        0.08,
        0.925,
        subtitle,
        ha="left",
        fontsize=9,
        color=MUTED_COLOR,
    )


def shade_dark_phase(ax, hours):
    """Shade ZT12-ZT24 for output that begins at ZT0."""
    # Dark phases begin at hours 12, 36, and 60 in a 72-hour trace.
    for start_h in np.arange(12.0, hours[-1] + 0.1, 24.0):
        # axvspan shades the full panel height between two x coordinates.
        ax.axvspan(
            start_h,
            # Do not let the final band extend beyond the available data.
            min(start_h + 12.0, float(hours[-1])),
            color="#F3F4F6",
            # A low z-order keeps the shading behind grids and curves.
            zorder=0,
        )


def save_figure(fig, output_path):
    # Create the destination folder if it does not already exist.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, facecolor="white")
    # Release memory because this script creates three separate figures.
    plt.close(fig)


def plot_timeseries(result):
    """Render the four temperature conditions and return their summaries."""
    # ``result`` is the dictionary returned by MODEL.run().
    hours = np.asarray(result["hour"], dtype=float)
    # These constant arrays need only their first value for the subtitle.
    low_e2 = float(result["constant_low_E2_pg_ml"][0])
    high_e2 = float(result["constant_high_E2_pg_ml"][0])
    # Join all four curves only to find common y-axis limits.
    all_values = np.concatenate([
        np.asarray(result[key], dtype=float) for key, _, _, _ in SERIES
    ])
    # Add 0.2 C of padding and round outward to the nearest 0.2 C.
    y_min = np.floor((float(np.min(all_values)) - 0.2) * 5.0) / 5.0
    y_max = np.ceil((float(np.max(all_values)) + 0.2) * 5.0) / 5.0

    # fig is the complete image; ax is its one rectangular plotting panel.
    fig, ax = plt.subplots(figsize=(10, 5.67), dpi=150)
    shade_dark_phase(ax, hours)

    summaries = {}
    # Use the key, label, color, and line style stored in SERIES for each curve.
    for key, label, color, linestyle in SERIES:
        values = np.asarray(result[key], dtype=float)
        ax.plot(
            hours,
            values,
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=2.2,
            zorder=3,
        )
        summaries[label] = {
            "mean": float(np.mean(values)),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            # np.ptp is peak-to-peak range: maximum minus minimum.
            "amplitude": float(np.ptp(values)),
        }

    ax.set(
        xlim=(0.0, float(hours[-1])),
        ylim=(y_min, y_max),
        xlabel="Time after warm-up (hours)",
        ylabel="Core temperature (°C)",
    )
    ax.set_xticks(np.arange(0.0, hours[-1] + 0.1, 12.0))
    ax.set_yticks(np.arange(np.ceil(y_min * 2) / 2, y_max + 0.01, 0.5))
    style_axes(ax)

    add_heading(
        fig,
        "Meeting 8 temperature model: four E2/activity conditions",
        f"Ta={AMBIENT_TEMPERATURE_C:g} °C; "
        f"B={OVX_BACKGROUND_HEAT_KCAL_H:.2f} kcal/h; "
        f"E2 heat={E2_WHOLE_BODY_THERMOGENIC_FRACTION:.2f} B; "
        f"aged/high E2={low_e2:.1f}/{high_e2:.2f} pg/mL; "
        f"Hill K={E2_HEAT_HALF_PG_ML:g}, "
        f"n={E2_HEAT_HILL_COEFFICIENT:g}; gray = dark phase",
    )
    ax.legend(
        # Anchor the legend just above the Axes so it does not cover data.
        loc="lower left",
        bbox_to_anchor=(0.0, 1.015),
        ncol=4,
        frameon=False,
        fontsize=8,
        handlelength=3.2,
        columnspacing=1.6,
    )
    # Reserve figure space above for the legend/title and below for labels.
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.13, top=0.79)
    save_figure(fig, DEFAULT_OUTPUT)
    return summaries


def plot_estrogen_timeseries(result):
    """Render the three distinct E2 inputs used by the four conditions."""
    hours = np.asarray(result["hour"], dtype=float)
    # Combine the three E2 arrays only to calculate one shared upper limit.
    all_values = np.concatenate([
        np.asarray(result[key], dtype=float)
        for key, _, _, _ in ESTROGEN_SERIES
    ])
    # Start at zero and round the padded maximum up to a multiple of 10 pg/mL.
    y_max = np.ceil((float(np.max(all_values)) + 2.0) / 10.0) * 10.0

    fig, ax = plt.subplots(figsize=(10, 5.67), dpi=150)
    shade_dark_phase(ax, hours)

    summaries = {}
    # This loop is the E2 counterpart of the temperature plotting loop.
    for key, label, color, linestyle in ESTROGEN_SERIES:
        values = np.asarray(result[key], dtype=float)
        ax.plot(
            hours,
            values,
            label=label,
            color=color,
            linestyle=linestyle,
            linewidth=2.2,
            zorder=3,
        )
        summaries[label] = {
            "mean": float(np.mean(values)),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
        }

    ax.set(
        xlim=(0.0, float(hours[-1])),
        ylim=(0.0, y_max),
        xlabel="Time after warm-up (hours)",
        ylabel="Estradiol (pg/mL)",
    )
    ax.set_xticks(np.arange(0.0, hours[-1] + 0.1, 12.0))
    ax.set_yticks(np.arange(0.0, y_max + 0.01, 10.0))
    style_axes(ax)

    add_heading(
        fig,
        "Meeting 8 estrogen model: E2 inputs by condition",
        "The constant-high E2 input is shared by the full- and "
        "reduced-activity conditions; gray = dark phase",
    )
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, 1.015),
        ncol=3,
        frameon=False,
        fontsize=8.5,
        handlelength=3.2,
        columnspacing=1.8,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.13, top=0.79)
    save_figure(fig, DEFAULT_ESTROGEN_OUTPUT)
    return summaries


def digitize_licinio_target():
    """Extract the black E2 trace from the Licinio reference-panel image."""
    # imread loads the entire PNG as (row, column, color-channel) pixels.
    # ``[..., 0]`` keeps channel 0 for every row and column. The source looks
    # grayscale, so this one channel is an adequate brightness map. PNG values
    # are 0-1; multiplying by 255 gives black=0 and white=255.
    gray = np.asarray(
        mpimg.imread(DEFAULT_LICINIO_TARGET)[..., 0], dtype=float
    ) * 255.0

    # Manually measured graph crop in full-image pixel coordinates:
    # x=196..754 spans the first 08:00 tick to the next 08:00 tick;
    # y=286..480 contains the black trace while excluding most labels/axes.
    x0, x1 = 196, 754
    y0, y1 = 286, 480
    # xs stores image columns; ys stores the detected trace row in each column.
    xs, ys = [], []
    for x in range(x0, x1 + 1):
        # Read one vertical column within the cropped graph region.
        column = gray[y0 : y1 + 1, x]
        # Return positions of all dark pixels; positions are relative to y0.
        indices = np.flatnonzero(column < 180)
        if indices.size:
            # Darker pixels receive more weight than gray antialiased edges.
            ink = (255.0 - column[indices]) ** 2
            xs.append(x)
            # np.average returns one weighted center row. Add y0 to convert
            # that crop-relative row back to a full-image row coordinate.
            ys.append(y0 + np.average(indices, weights=ink))

    # Build one x value for every graph column. Interpolation fills columns in
    # which no dark trace pixel was detected between neighboring known points.
    all_x = np.arange(x0, x1 + 1)
    all_y = np.interp(all_x, np.asarray(xs), np.asarray(ys))
    # Smooth one-pixel scan noise: 25% previous + 50% current + 25% next.
    # Edge padding keeps the smoothed array the same length as all_x.
    all_y = np.convolve(
        np.pad(all_y, 1, mode="edge"),
        np.array([0.25, 0.5, 0.25]),
        mode="valid",
    )
    # Convert horizontal pixels to 0-24 hours.
    time_h = (all_x - x0) / (x1 - x0) * 24.0
    # Manual y-axis calibration: row 516=0 pg/mL and row 335=40 pg/mL.
    # Image rows increase downward, hence the reversed subtraction.
    e2 = (516.0 - all_y) / (516.0 - 335.0) * 40.0
    return time_h, e2


def prepare_licinio_comparison():
    """Return the Licinio trace and three successive 24-hour model days."""
    target_t, target_e2 = digitize_licinio_target()
    # Request 73 post-warm-up hours so the first hour can be skipped while
    # retaining a complete 72-hour (three-day) comparison window.
    endocrine = simulate_endocrine(output_h=73.0)
    # Warm-up ends at lights-on (ZT0/07:00); +1 hour starts at ZT1/08:00,
    # matching the first clock label on the Licinio graph.
    start_h = WARMUP_H + 1.0
    # Boolean mask: keep only simulation times from start_h through +72 hours.
    keep = (
        (endocrine["time_h"] >= start_h)
        & (endocrine["time_h"] <= start_h + 72.0)
    )
    # Reset the retained time window so it runs from 0 through 72 hours.
    model_t = endocrine["time_h"][keep] - start_h
    model_e2 = endocrine["E2"][keep]

    days = []
    # Split 72 hours into three 24-hour arrays and reset each day to 0-24.
    for day in range(3):
        day_keep = (
            (model_t >= 24.0 * day)
            & (model_t <= 24.0 * (day + 1))
        )
        # Each list item is a (daily time array, daily E2 array) tuple.
        days.append((
            model_t[day_keep] - 24.0 * day,
            model_e2[day_keep],
        ))
    return target_t, target_e2, days


def format_licinio_axis(ax):
    # Every overlaid day uses the same 24-hour and 0-60-pg/mL scales.
    ax.set_xlim(0.0, 24.0)
    ax.set_ylim(0.0, 60.0)
    ax.set_xticks((0, 6, 12, 18, 24))
    # Model hour 0 is ZT1/08:00, so show civil-clock labels instead of 0-24.
    ax.set_xticklabels(("08:00", "14:00", "20:00", "02:00", "08:00"))
    ax.set_yticks(np.arange(0.0, 61.0, 10.0))
    style_axes(ax)


def day_summary(day_number, values):
    # Text displayed beneath the folded graph for each colored model day.
    return (
        f"Day {day_number}: mean {values.mean():.1f}\n"
        f"range {values.min():.1f}-{values.max():.1f} pg/mL"
    )


def plot_licinio_folded(comparison):
    """Fold three model days over one digitized 24-hour Licinio trace."""
    # comparison contains the digitized reference plus three (time, E2) pairs.
    target_t, target_e2, days = comparison

    # One Figure with one Axes: all four traces are overlaid in one panel.
    fig, ax = plt.subplots(figsize=(10, 6.08), dpi=150)
    # Plot the digitized published example as a dashed nearly-black line.
    ax.plot(
        target_t,
        target_e2,
        color="#23272D",
        linestyle=(0, (4, 3)),
        linewidth=1.6,
        label="Licinio",
    )
    # zip pairs day 1/2/3 with blue/orange/purple; enumerate supplies labels.
    for index, ((day_t, day_e2), color) in enumerate(
        zip(days, DAY_COLORS), start=1
    ):
        ax.plot(
            day_t,
            day_e2,
            color=color,
            linewidth=1.8,
            label=f"Model day {index}",
        )

    format_licinio_axis(ax)
    ax.set_xlabel("Time of day", fontsize=11, fontweight="bold")
    ax.set_ylabel("Estradiol (pg/mL)", fontsize=11, fontweight="bold")
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, 1.03),
        ncol=4,
        frameon=False,
        fontsize=8.5,
        handlelength=3.2,
        columnspacing=2.5,
    )
    add_heading(
        fig,
        "Licinio vs model: three simulated days folded onto 24 hours",
        "Dashed black = digitized Licinio example; colored = clock-free "
        "endocrine model days; pituitary E2 "
        f"K={PITUITARY_E2_HALF_PG_ML:g} pg/mL, "
        f"n={PITUITARY_E2_HILL_COEFFICIENT:g}",
    )
    # These are figure-relative x positions for the three summary text blocks:
    # 0=left edge, 0.5=center, and 1=right edge of the complete image.
    for center, day_number, (_, values), color in zip(
        (0.22, 0.51, 0.80),
        range(1, 4),
        days,
        DAY_COLORS,
    ):
        fig.text(
            center,
            0.045,
            day_summary(day_number, values),
            ha="center",
            color=color,
            fontsize=7.5,
            fontweight="bold",
        )

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.17, top=0.78)
    save_figure(fig, DEFAULT_LICINIO_OUTPUT)


def main():
    # Fixed workflow: run MODEL.py once, then write the three standard PNGs.
    configure_style()
    result = run()
    # Reuse the same model result for both primary time-series figures.
    temperature_summaries = plot_timeseries(result)
    estrogen_summaries = plot_estrogen_timeseries(result)
    # Prepare the digitized/model comparison once and render the overlay.
    comparison = prepare_licinio_comparison()
    plot_licinio_folded(comparison)

    for label, values in temperature_summaries.items():
        print(
            f"{label}: mean={values['mean']:.4f} C, "
            f"min={values['minimum']:.4f} C, "
            f"max={values['maximum']:.4f} C, "
            f"amplitude={values['amplitude']:.4f} C"
        )
    print(f"Wrote: {DEFAULT_OUTPUT}")
    for label, values in estrogen_summaries.items():
        print(
            f"{label}: mean={values['mean']:.4f} pg/mL, "
            f"min={values['minimum']:.4f} pg/mL, "
            f"max={values['maximum']:.4f} pg/mL"
        )
    print(f"Wrote: {DEFAULT_ESTROGEN_OUTPUT}")
    print(f"Wrote: {DEFAULT_LICINIO_OUTPUT}")


if __name__ == "__main__":
    main()
