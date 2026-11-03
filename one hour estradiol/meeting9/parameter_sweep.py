"""Interactive parameter sweep for the Meeting 9 HPG model."""

import numpy as np

import MODEL


SLIDER_SPECS = (
    (
        "n_gain",
        "N gain (n_gain): maximum target level of the E2-inhibited N signal",
        r"$N^*=\mathbf{G_N}\,K_{E2\to N}^5/(K_{E2\to N}^5+E2^5)$",
        1,
        200,
    ),
    (
        "gnrh_pulse_max_jump",
        "GnRH pulse gain (gnrh_pulse_max_jump): largest possible GnRH pulse",
        r"$\Delta GnRH=\mathbf{A_{GnRH}}\,N^5/(K_A^5+N^5)$",
        1,
        200,
    ),
    (
        "lh_gnrh_gain",
        "LH gain (lh_gnrh_gain): LH production driven by GnRH",
        r"$dLH/dt=\mathbf{G_{LH}}\,GnRH-k_{LH}LH$",
        1,
        30,
    ),
    (
        "e2_gain",
        "E2 gain (e2_gain): E2 production driven by LH",
        r"$dE2/dt=\mathbf{G_{E2}}\,LH-k_{E2}E2$",
        1,
        30,
    ),
    (
        "n_e2_half_pg_ml",
        "E2 -> N midpoint (n_e2_half_pg_ml): E2 giving half-maximal N target",
        r"$N^*=G_N\,\mathbf{K_{E2\to N}}^5/(\mathbf{K_{E2\to N}}^5+E2^5)$",
        1,
        60,
    ),
)

LIVE_REDRAW_INTERVAL_MS = 250


def find_periodicity(time, values, tolerance=0.05, maximum_period=30):
    """Return the detected peaks and the smallest repeating orbit.

    Every comparable peak height and interpeak interval must agree with the
    value one candidate orbit earlier. Differences are scaled by the shared
    median height or interval, so every observation receives the same absolute
    tolerance regardless of its size.
    """

    peaks = np.flatnonzero(
        (values[1:-1] > values[:-2]) & (values[1:-1] >= values[2:])
    ) + 1
    heights = values[peaks]
    intervals = np.diff(time[peaks])
    count = min(len(heights), len(intervals))

    if count < 6:
        return peaks, 0

    heights = heights[-count:]
    intervals = intervals[-count:]
    for period in range(1, min(maximum_period, count // 2) + 1):
        height_errors = (
            np.abs(heights[period:] - heights[:-period]) / np.median(heights)
        )
        interval_errors = (
            np.abs(intervals[period:] - intervals[:-period]) / np.median(intervals)
        )
        if np.all(height_errors <= tolerance) and np.all(
            interval_errors <= tolerance
        ):
            return peaks, period

    return peaks, -1


def period_label(period):
    if period == 0:
        return "not enough peaks"
    if period == -1:
        return ">30 (high-period)"
    return f"period-{period}"


def run_model(parameters, step_h=MODEL.STEP_H):
    result = MODEL.simulate_endocrine(step_h=step_h, **parameters)
    output_points = round(MODEL.OUTPUT_H / step_h) + 1
    return (
        result["time_h"][-output_points:] - MODEL.WARMUP_H,
        result["E2"][-output_points:],
        result["LH"][-output_points:],
    )


def open_window():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button, Slider
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "This program needs Matplotlib. Install it with:\n"
            "    python -m pip install matplotlib"
        ) from exc

    figure, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    figure.subplots_adjust(
        left=0.09, right=0.97, top=0.86, bottom=0.49, hspace=0.38
    )
    figure.suptitle(
        "Meeting 9 parameter sweep: E2 inhibits N; N controls GnRH pulse timing and size",
        y=0.965,
    )
    figure.text(
        0.5,
        0.395,
        r"Fixed: $dN/dt=N^*-N$; $K_{N\to period}=K_{N\to amplitude}=1$; Hill coefficients = 5",
        ha="center",
        fontsize=9,
    )
    defaults = {
        "n_gain": MODEL.N_GAIN,
        "n_e2_half_pg_ml": MODEL.N_E2_HALF_PG_ML,
        "gnrh_pulse_max_jump": MODEL.GNRH_PULSE_MAX_JUMP,
        "lh_gnrh_gain": MODEL.LH_GNRH_GAIN,
        "e2_gain": MODEL.E2_GAIN,
    }
    sliders = {}
    for index, (field, description, equation, lower, upper) in enumerate(SLIDER_SPECS):
        row_y = 0.33 - index * 0.058
        figure.text(
            0.03,
            row_y + 0.012,
            f"{description}\n{equation}",
            va="center",
            fontsize=8.5,
        )
        slider_axis = figure.add_axes([0.70, row_y, 0.20, 0.024])
        sliders[field] = Slider(
            slider_axis,
            "",
            lower,
            upper,
            valinit=defaults[field],
            valstep=1,
            valfmt="%0.0f",
        )

    reset_button = Button(figure.add_axes([0.43, 0.025, 0.14, 0.04]), "Reset")
    redraw_state = {
        "live_pending": False,
    }

    def draw():
        parameters = {
            field: slider.val for field, slider in sliders.items()
        }
        time_h, e2, lh = run_model(parameters, MODEL.STEP_H)
        e2_peaks, e2_period = find_periodicity(time_h, e2)
        lh_peaks, lh_period = find_periodicity(time_h, lh)
        keep_every = max(1, round((3 / 60) / MODEL.STEP_H))

        for axis in axes:
            axis.clear()
            axis.grid(alpha=0.2)

        axes[0].axhspan(7, 60, color="green", alpha=0.1)
        axes[0].plot(time_h[::keep_every], e2[::keep_every], color="tab:red")
        axes[0].plot(
            time_h[e2_peaks], e2[e2_peaks], ".", color="black", markersize=4
        )
        axes[0].set_title(f"E2 orbit: {period_label(e2_period)}")
        axes[0].set_ylabel("E2 (pg/mL)")

        axes[1].axhspan(1, 10, color="green", alpha=0.1)
        axes[1].plot(time_h[::keep_every], lh[::keep_every], color="tab:blue")
        axes[1].plot(
            time_h[lh_peaks], lh[lh_peaks], ".", color="black", markersize=4
        )
        axes[1].set_title(f"LH orbit: {period_label(lh_period)}")
        axes[1].set_ylabel("LH (IU/L)")
        axes[1].set_xlabel("Time (hours)")
        axes[1].set_xlim(0, MODEL.OUTPUT_H)

        figure.canvas.draw_idle()

    def queue_live_draw(_):
        redraw_state["live_pending"] = True

    def draw_pending_live():
        if redraw_state["live_pending"]:
            redraw_state["live_pending"] = False
            draw()

    for slider in sliders.values():
        slider.on_changed(queue_live_draw)

    live_timer = figure.canvas.new_timer(interval=LIVE_REDRAW_INTERVAL_MS)
    live_timer.add_callback(draw_pending_live)
    live_timer.start()

    def finish_drag(_):
        if redraw_state["live_pending"]:
            redraw_state["live_pending"] = False
            draw()

    figure.canvas.mpl_connect("button_release_event", finish_drag)

    def reset(_):
        for slider in sliders.values():
            slider.eventson = False
            slider.reset()
            slider.eventson = True
        redraw_state["live_pending"] = False
        draw()

    reset_button.on_clicked(reset)
    draw()
    plt.show()


def main():
    open_window()


if __name__ == "__main__":
    main()
