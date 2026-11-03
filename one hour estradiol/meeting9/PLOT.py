"""Plot Meeting 9 estradiol and its predicted core-temperature response.

Both the endocrine trajectory and the imported Meeting 8 temperature model
come from the local ``MODEL.py``. Temperature does not feed back into the
endocrine model.

Running this file writes ``temperature_timeseries.png`` beside it.
"""

from pathlib import Path

import numpy as np

import MODEL


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "temperature_timeseries.png"


def run_temperature():
    """Return the temperature trajectory produced by the local model."""

    trajectory = MODEL.run()
    scn_bat_drive = MODEL.scn_drive(trajectory["zt_h"])
    effective_drive = (
        trajectory["E2_pg_ml"]
        + trajectory["scn_e2_scale_pg_ml"] * scn_bat_drive
    )
    return {
        "hour": trajectory["hour"],
        "E2": trajectory["E2_pg_ml"],
        "effective_drive": effective_drive,
        "scn_bat_drive": scn_bat_drive,
        "temperature": trajectory["temperature_c"],
        "e2_heat_half": MODEL.E2_HEAT_HALF_PG_ML,
        "e2_heat_hill": MODEL.E2_HEAT_HILL_COEFFICIENT,
        "scn_e2_scale": trajectory["scn_e2_scale_pg_ml"],
    }


def shade_dark_phase(axis, final_hour):
    """Shade the ZT12-ZT24 dark phase behind a 72-hour trace."""

    for start in np.arange(12.0, final_hour + 0.1, 24.0):
        axis.axvspan(start, min(start + 12.0, final_hour), color="0.95")


def plot_temperature(result, output_path=OUTPUT):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This program needs Matplotlib. Install it with:\n"
            "    python -m pip install matplotlib"
        ) from exc

    hour = result["hour"]
    e2 = result["E2"]
    effective_drive = result["effective_drive"]
    scn_bat_drive = result["scn_bat_drive"]
    temperature = result["temperature"]
    keep_every = max(1, round((1 / 60) / MODEL.STEP_H))

    figure, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for axis in axes:
        shade_dark_phase(axis, float(hour[-1]))
        axis.grid(alpha=0.25)

    axes[0].plot(
        hour[::keep_every],
        temperature[::keep_every],
        color="tab:blue",
        linewidth=1.5,
    )
    axes[0].set_ylabel("Core temperature (°C)")
    axes[0].set_title("Predicted temperature")

    axes[1].plot(
        hour[::keep_every],
        e2[::keep_every],
        color="tab:red",
        linewidth=1.2,
        label="E2",
    )
    axes[1].plot(
        hour[::keep_every],
        effective_drive[::keep_every],
        color="tab:orange",
        linewidth=1.2,
        label="E2 + k_SCN*SCN_BAT",
    )
    axes[1].set_ylabel("Thermogenic input (pg/mL)")
    axes[1].set_title("E2 and shared Hill input")
    axes[1].legend(loc="upper right")

    axes[2].plot(
        hour[::keep_every],
        scn_bat_drive[::keep_every],
        color="tab:green",
        linewidth=1.2,
    )
    axes[2].set_ylabel("SCN BAT drive")
    axes[2].set_xlabel("Time after warm-up (hours)")
    axes[2].set_xlim(0, float(hour[-1]))
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].set_title("Smooth SCN-associated BAT contribution")

    figure.suptitle("Meeting 9 shared E2-SCN thermogenesis model")
    figure.text(
        0.5,
        0.92,
        f"Shared E2 + k_SCN*SCN_BAT Hill: K={result['e2_heat_half']:g} pg/mL, "
        f"n={result['e2_heat_hill']:g}, "
        f"k_SCN={result['scn_e2_scale']:g} pg/mL; gray = dark phase",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    figure.subplots_adjust(left=0.11, right=0.98, bottom=0.08, top=0.88, hspace=0.38)
    figure.savefig(output_path, dpi=180, facecolor="white")
    plt.close(figure)


def main():
    result = run_temperature()
    plot_temperature(result)
    print(
        "E2 min/mean/max: "
        f"{result['E2'].min():.3f}/{result['E2'].mean():.3f}/"
        f"{result['E2'].max():.3f} pg/mL"
    )
    print(
        "Temperature min/mean/max: "
        f"{result['temperature'].min():.3f}/"
        f"{result['temperature'].mean():.3f}/"
        f"{result['temperature'].max():.3f} °C"
    )
    print(f"Wrote: {OUTPUT}")


if __name__ == "__main__":
    main()
