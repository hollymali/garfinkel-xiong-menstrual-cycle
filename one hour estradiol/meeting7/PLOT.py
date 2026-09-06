"""Plot the activity input and temperature output from ``MODEL.py``."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import MODEL


def shade_dark_phase(axis, duration_h):
    """Shade ZT12-ZT24 for every displayed day."""
    day_start_h = 0.0
    while day_start_h < duration_h:
        axis.axvspan(
            day_start_h + MODEL.LIGHTS_OFF_ZT,
            min(day_start_h + MODEL.ZT_PERIOD_H, duration_h),
            color="0.92",
            zorder=0,
        )
        day_start_h += MODEL.ZT_PERIOD_H


def make_figure():
    """Return a figure containing the editable activity and CBT time series."""
    result = MODEL.run()
    hour = result["hour"]
    duration_h = float(hour[-1])

    figure, (activity_axis, temperature_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": (1, 2.4)},
        constrained_layout=True,
    )

    for axis in (activity_axis, temperature_axis):
        shade_dark_phase(axis, duration_h)

    activity_axis.step(
        hour,
        result["activity_drive"],
        where="post",
        label="Non-pregnant activity",
    )
    activity_axis.step(
        hour,
        MODEL.PREGNANT_ACTIVITY_MULTIPLIER
        * result["pregnant_activity_drive"],
        where="post",
        linestyle="--",
        label="Pregnant activity",
    )
    activity_axis.set_ylabel("Relative activity")
    activity_axis.set_ylim(-0.05, 1.1)
    activity_axis.legend(frameon=False, ncol=2)

    temperature_axis.plot(
        hour,
        result["nonpregnant_temperature_c"],
        label="Non-pregnant",
    )
    temperature_axis.plot(
        hour,
        result["constant_high_temperature_c"],
        label="Constant-high E2 only",
    )
    temperature_axis.plot(
        hour,
        result["pregnant_temperature_c"],
        label="Pregnant: P4 +0.3°C; SCN–BAT effect × 0.51",
    )
    temperature_axis.set_xlabel(
        "Elapsed Zeitgeber time (h; ZT0 at 0, 24, 48, and 72)"
    )
    temperature_axis.set_ylabel("Core body temperature (°C)")
    temperature_axis.legend(frameon=False)

    figure.suptitle("Meeting 7 reduced temperature model")
    return figure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--save",
        type=Path,
        help="Save the figure to this path instead of opening it.",
    )
    args = parser.parse_args()

    figure = make_figure()
    if args.save is None:
        plt.show()
    else:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.save, dpi=180)
        print(args.save.resolve())


if __name__ == "__main__":
    main()
