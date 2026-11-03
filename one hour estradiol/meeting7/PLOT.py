"""Plot all four Meeting-7 core-temperature time series.

The constant-aged-E2 condition is fixed at 13.3 pg/mL, the reported group
mean for intact 24-month-old female C57BL/6J mice. Run this file to write
``temperature_timeseries_four_conditions.png`` beside it.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from MODEL import E2_HEAT_HALF_PG_ML, E2_HEAT_HILL_COEFFICIENT, run


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "temperature_timeseries_four_conditions.png"
WIDTH = 1500
HEIGHT = 850
MARGIN_LEFT = 115
MARGIN_RIGHT = 45
MARGIN_TOP = 175
MARGIN_BOTTOM = 95

SERIES = (
    ("constant_low_temperature_c", "Constant aged E2", "#16803B"),
    ("constant_high_temperature_c", "Constant high E2", "#7C3AED"),
    ("nonpregnant_temperature_c", "Cycling E2 / nonpregnant", "#2563EB"),
    ("pregnant_temperature_c", "Pregnant", "#C2410C"),
)


def font(size: int, bold: bool = False):
    """Load a common Windows font, with Pillow's built-in font as fallback."""
    names = (
        ("arialbd.ttf", "Arial Bold.ttf")
        if bold
        else ("arial.ttf", "Arial.ttf")
    )
    for name in names:
        try:
            return ImageFont.truetype(name, size=size)
        except OSError:
            pass
    return ImageFont.load_default()


def dashed_line(draw, points, fill, width=3, dash=14, gap=9):
    """Draw a dashed polyline between densely sampled points."""
    for start in range(0, len(points) - 1, dash + gap):
        stop = min(start + dash, len(points))
        if stop - start >= 2:
            draw.line(points[start:stop], fill=fill, width=width, joint="curve")


def plot_timeseries(output_path: Path = DEFAULT_OUTPUT):
    """Run the model, render the four-condition plot, and return summaries."""
    result = run()
    hours = np.asarray(result["hour"], dtype=float)
    low_e2 = float(result["constant_low_E2_pg_ml"][0])
    high_e2 = float(result["constant_high_E2_pg_ml"][0])
    all_values = np.concatenate([
        np.asarray(result[key], dtype=float) for key, _, _ in SERIES
    ])
    y_min = np.floor((float(np.min(all_values)) - 0.2) * 5.0) / 5.0
    y_max = np.ceil((float(np.max(all_values)) + 0.2) * 5.0) / 5.0
    plot_width = WIDTH - MARGIN_LEFT - MARGIN_RIGHT
    plot_height = HEIGHT - MARGIN_TOP - MARGIN_BOTTOM

    def x_pixel(hour):
        return MARGIN_LEFT + plot_width * hour / float(hours[-1])

    def y_pixel(value):
        return MARGIN_TOP + plot_height * (y_max - value) / (y_max - y_min)

    image = Image.new("RGB", (WIDTH, HEIGHT), "#FFFFFF")
    draw = ImageDraw.Draw(image)
    title_font = font(34, bold=True)
    subtitle_font = font(21)
    axis_font = font(20)
    tick_font = font(17)
    legend_font = font(19)

    draw.text(
        (MARGIN_LEFT, 34),
        "Meeting 7 temperature model: four E2/pregnancy conditions",
        fill="#111827",
        font=title_font,
    )
    draw.text(
        (MARGIN_LEFT, 86),
        f"Aged E2={low_e2:.1f}, high E2={high_e2:.2f} pg/mL; "
        f"thermal Hill n={E2_HEAT_HILL_COEFFICIENT:g}, "
        f"K={E2_HEAT_HALF_PG_ML:g} pg/mL; gray = dark phase",
        fill="#4B5563",
        font=subtitle_font,
    )

    # Dark phase is ZT12-ZT24; retained output starts at ZT0.
    for start_h in np.arange(12.0, hours[-1] + 0.1, 24.0):
        draw.rectangle(
            [
                (x_pixel(start_h), MARGIN_TOP),
                (x_pixel(min(start_h + 12.0, hours[-1])), MARGIN_TOP + plot_height),
            ],
            fill="#F3F4F6",
        )

    y_ticks = np.arange(np.ceil(y_min * 2) / 2, y_max + 0.01, 0.5)
    for value in y_ticks:
        y = y_pixel(value)
        draw.line(
            [(MARGIN_LEFT, y), (MARGIN_LEFT + plot_width, y)],
            fill="#D1D5DB",
            width=1,
        )
        label = f"{value:.1f}"
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text(
            (MARGIN_LEFT - 14 - (box[2] - box[0]), y - 10),
            label,
            fill="#374151",
            font=tick_font,
        )

    for hour in np.arange(0.0, hours[-1] + 0.1, 12.0):
        x = x_pixel(hour)
        draw.line(
            [(x, MARGIN_TOP), (x, MARGIN_TOP + plot_height)],
            fill="#E5E7EB",
            width=1,
        )
        label = f"{hour:.0f}"
        box = draw.textbbox((0, 0), label, font=tick_font)
        draw.text(
            (x - (box[2] - box[0]) / 2, MARGIN_TOP + plot_height + 15),
            label,
            fill="#374151",
            font=tick_font,
        )

    draw.rectangle(
        [
            (MARGIN_LEFT, MARGIN_TOP),
            (MARGIN_LEFT + plot_width, MARGIN_TOP + plot_height),
        ],
        outline="#6B7280",
        width=2,
    )

    # Plot every minute; this preserves shape while keeping rendering quick.
    stride = max(1, round((1.0 / 60.0) / (hours[1] - hours[0])))
    sample_indices = np.arange(0, len(hours), stride, dtype=int)
    if sample_indices[-1] != len(hours) - 1:
        sample_indices = np.append(sample_indices, len(hours) - 1)

    summaries = {}
    for series_index, (key, label, color) in enumerate(SERIES):
        values = np.asarray(result[key], dtype=float)
        points = [
            (x_pixel(float(hours[i])), y_pixel(float(values[i])))
            for i in sample_indices
        ]
        if series_index == 0:
            dashed_line(draw, points, color, width=4)
        else:
            draw.line(points, fill=color, width=4, joint="curve")
        summaries[label] = {
            "mean": float(np.mean(values)),
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            "amplitude": float(np.ptp(values)),
        }

    legend_positions = (
        (MARGIN_LEFT, 128),
        (MARGIN_LEFT + 330, 128),
        (MARGIN_LEFT + 650, 128),
        (MARGIN_LEFT + 1035, 128),
    )
    for series_index, ((_, label, color), position) in enumerate(
        zip(SERIES, legend_positions)
    ):
        x, y = position
        if series_index == 0:
            for segment in range(0, 44, 16):
                draw.line([(x + segment, y + 11), (x + segment + 10, y + 11)], fill=color, width=4)
        else:
            draw.line([(x, y + 11), (x + 44, y + 11)], fill=color, width=4)
        draw.text((x + 55, y), label, fill="#111827", font=legend_font)

    x_label = "Time after warm-up (hours)"
    x_box = draw.textbbox((0, 0), x_label, font=axis_font)
    draw.text(
        (
            MARGIN_LEFT + (plot_width - (x_box[2] - x_box[0])) / 2,
            HEIGHT - 42,
        ),
        x_label,
        fill="#111827",
        font=axis_font,
    )
    y_label = "Core temperature (°C)"
    label_image = Image.new("RGBA", (330, 50), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_image)
    label_draw.text((0, 0), y_label, fill="#111827", font=axis_font)
    label_image = label_image.rotate(90, expand=True)
    image.paste(
        label_image,
        (20, MARGIN_TOP + (plot_height - label_image.height) // 2),
        label_image,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    image.save(output_path, dpi=(180, 180))
    return low_e2, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="PNG output path",
    )
    args = parser.parse_args()
    low_e2, summaries = plot_timeseries(args.output)
    print(f"Constant aged E2: {low_e2:.6f} pg/mL")
    for label, values in summaries.items():
        print(
            f"{label}: mean={values['mean']:.4f} C, "
            f"min={values['minimum']:.4f} C, "
            f"max={values['maximum']:.4f} C, "
            f"amplitude={values['amplitude']:.4f} C"
        )
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
