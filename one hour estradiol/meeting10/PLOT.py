"""Plot the first-edition sleep and SCN-modulated HPG time series.

Running this file writes ``hpg_time_series.png`` beside it. The retained
48-hour trace begins at ZT0, which is 07:00 in the project clock. Pillow is
used so the plot has no dependency beyond the bundled scientific runtime.
"""

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import MODEL


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "hpg_time_series.png"

SCALE = 2
WIDTH = 1200
HEIGHT = 1080
LEFT = 112
RIGHT = 36
TOP = 120
BOTTOM = 72
PANEL_GAP = 18

WHITE = "#ffffff"
TEXT = "#172033"
MUTED = "#5f6878"
GRID = "#d9dee8"
SLEEP = "#e9edf5"
SCN_COLOR = "#7c3aed"
SLEEP_COLOR = "#526173"
N_COLOR = "#16803b"
GNRH_COLOR = "#c2410c"
LH_COLOR = "#2563eb"
E2_COLOR = "#b4235a"


def font(size, bold=False):
    name = "arialbd.ttf" if bold else "arial.ttf"
    path = Path("C:/Windows/Fonts") / name
    return ImageFont.truetype(str(path), size * SCALE)


def px(value):
    return int(round(value * SCALE))


def clock_label(hour_after_zt0):
    civil_hour = int(round((7.0 + hour_after_zt0) % 24.0))
    return f"{civil_hour:02d}:00"


def scale_x(hour, final_hour):
    return px(LEFT + (WIDTH - LEFT - RIGHT) * hour / final_hour)


def padded_limits(values):
    lower = float(np.min(values))
    upper = float(np.max(values))
    span = upper - lower
    padding = 0.06 * span if span > 0.0 else max(abs(lower) * 0.05, 1.0)
    return lower - padding, upper + padding


def scale_y(value, lower, upper, panel_top, panel_bottom):
    fraction = (float(value) - lower) / (upper - lower)
    return px(panel_bottom - fraction * (panel_bottom - panel_top))


def draw_polyline(draw, hour, values, color, width, bounds, final_hour):
    lower, upper, panel_top, panel_bottom = bounds
    points = [
        (
            scale_x(float(x_value), final_hour),
            scale_y(y_value, lower, upper, panel_top, panel_bottom),
        )
        for x_value, y_value in zip(hour, values)
    ]
    draw.line(points, fill=color, width=px(width), joint="curve")


def draw_panel(
    draw,
    panel_top,
    panel_bottom,
    hour,
    series,
    y_label,
    final_hour,
    limits=None,
):
    all_values = np.concatenate([np.asarray(values) for values, _, _ in series])
    lower, upper = limits if limits is not None else padded_limits(all_values)
    x_left = px(LEFT)
    x_right = px(WIDTH - RIGHT)
    y_top = px(panel_top)
    y_bottom = px(panel_bottom)

    for sleep_start in np.arange(16.0, final_hour + 24.0, 24.0):
        sleep_end = min(sleep_start + 8.0, final_hour)
        if sleep_start < final_hour:
            draw.rectangle(
                (
                    scale_x(sleep_start, final_hour),
                    y_top,
                    scale_x(sleep_end, final_hour),
                    y_bottom,
                ),
                fill=SLEEP,
            )

    for fraction in (0.0, 0.5, 1.0):
        y = px(panel_bottom - fraction * (panel_bottom - panel_top))
        draw.line((x_left, y, x_right, y), fill=GRID, width=px(1))
        value = lower + fraction * (upper - lower)
        label = f"{value:.1f}"
        box = draw.textbbox((0, 0), label, font=font(8))
        draw.text(
            (x_left - px(9) - (box[2] - box[0]), y - px(5)),
            label,
            fill=MUTED,
            font=font(8),
        )

    draw.line((x_left, y_top, x_left, y_bottom), fill=MUTED, width=px(1))
    draw.line((x_left, y_bottom, x_right, y_bottom), fill=MUTED, width=px(1))

    for values, color, line_width in series:
        draw_polyline(
            draw,
            hour,
            values,
            color,
            line_width,
            (lower, upper, panel_top, panel_bottom),
            final_hour,
        )

    draw.text(
        (px(15), px((panel_top + panel_bottom) / 2 - 8)),
        y_label,
        fill=TEXT,
        font=font(9, bold=True),
    )
    return lower, upper


def plot(output_path=OUTPUT):
    trajectory = MODEL.run(output_h=48.0)
    hour_full = trajectory["hour"]
    final_hour = float(hour_full[-1])
    keep_every = max(1, round((1.0 / 60.0) / MODEL.STEP_H))
    hour = hour_full[::keep_every]

    image = Image.new("RGB", (px(WIDTH), px(HEIGHT)), WHITE)
    draw = ImageDraw.Draw(image)

    draw.text(
        (px(LEFT), px(24)),
        "Sleep- and SCN-modulated human HPG model",
        fill=TEXT,
        font=font(17, bold=True),
    )
    draw.text(
        (px(LEFT), px(56)),
        "Gray = sleep (23:00–07:00): GnRH period ×1.57.  "
        "SCN gain peaks at 08:00/ZT1: GnRH pulse amplitude ×2.",
        fill=MUTED,
        font=font(9),
    )

    available = HEIGHT - TOP - BOTTOM - 4 * PANEL_GAP
    panel_height = available / 5.0
    panels = []
    for index in range(5):
        panel_top = TOP + index * (panel_height + PANEL_GAP)
        panels.append((panel_top, panel_top + panel_height))

    draw_panel(
        draw,
        *panels[0],
        hour,
        (
            (trajectory["SCN_drive"][::keep_every], SCN_COLOR, 2),
            (
                trajectory["sleep_period_multiplier"][::keep_every],
                SLEEP_COLOR,
                2,
            ),
        ),
        "Multiplier",
        final_hour,
        limits=(0.9, 2.1),
    )
    draw_panel(
        draw,
        *panels[1],
        hour,
        ((trajectory["N"][::keep_every], N_COLOR, 2),),
        "N state",
        final_hour,
    )
    draw_panel(
        draw,
        *panels[2],
        hour,
        ((trajectory["GnRH"][::keep_every], GNRH_COLOR, 1),),
        "GnRH",
        final_hour,
    )
    draw_panel(
        draw,
        *panels[3],
        hour,
        ((trajectory["LH_IU_L"][::keep_every], LH_COLOR, 2),),
        "LH (IU/L)",
        final_hour,
    )
    draw_panel(
        draw,
        *panels[4],
        hour,
        ((trajectory["E2_pg_ml"][::keep_every], E2_COLOR, 2),),
        "E2 (pg/mL)",
        final_hour,
    )

    legend_y = px(91)
    draw.line(
        (px(LEFT), legend_y, px(LEFT + 28), legend_y),
        fill=SCN_COLOR,
        width=px(3),
    )
    draw.text(
        (px(LEFT + 35), px(84)),
        "SCN amplitude gain",
        fill=TEXT,
        font=font(8),
    )
    draw.line(
        (px(LEFT + 190), legend_y, px(LEFT + 218), legend_y),
        fill=SLEEP_COLOR,
        width=px(3),
    )
    draw.text(
        (px(LEFT + 225), px(84)),
        "Sleep period gain",
        fill=TEXT,
        font=font(8),
    )

    bottom = panels[-1][1]
    for tick in np.arange(0.0, final_hour + 0.1, 4.0):
        x = scale_x(tick, final_hour)
        draw.line(
            (x, px(bottom), x, px(bottom + 5)),
            fill=MUTED,
            width=px(1),
        )
        label = clock_label(tick)
        box = draw.textbbox((0, 0), label, font=font(8))
        draw.text(
            (x - (box[2] - box[0]) / 2, px(bottom + 9)),
            label,
            fill=MUTED,
            font=font(8),
        )

    axis_label = "Civil time (ZT0 = 07:00)"
    box = draw.textbbox((0, 0), axis_label, font=font(9, bold=True))
    draw.text(
        ((px(WIDTH) - (box[2] - box[0])) / 2, px(HEIGHT - 28)),
        axis_label,
        fill=TEXT,
        font=font(9, bold=True),
    )

    image = image.resize((WIDTH, HEIGHT), Image.Resampling.LANCZOS)
    image.save(output_path)
    return trajectory


def main():
    trajectory = plot()
    phi = trajectory["phi"]
    pulse_indices = np.flatnonzero(np.diff(np.floor(phi)) > 0.0) + 1
    pulse_zt = trajectory["zt_h"][pulse_indices]
    sleep = (pulse_zt >= 16.0) & (pulse_zt < 24.0)
    print(f"Pulses over 48 h wake/sleep: {(~sleep).sum()}/{sleep.sum()}")
    print(
        "LH min/mean/max: "
        f"{trajectory['LH_IU_L'].min():.3f}/"
        f"{trajectory['LH_IU_L'].mean():.3f}/"
        f"{trajectory['LH_IU_L'].max():.3f} IU/L"
    )
    print(
        "E2 min/mean/max: "
        f"{trajectory['E2_pg_ml'].min():.3f}/"
        f"{trajectory['E2_pg_ml'].mean():.3f}/"
        f"{trajectory['E2_pg_ml'].max():.3f} pg/mL"
    )
    print(f"Wrote: {OUTPUT}")


if __name__ == "__main__":
    main()
