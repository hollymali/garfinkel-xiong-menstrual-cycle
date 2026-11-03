"""Live E2/LH plots, E2 FFT power, and saved one-at-a-time bands.

Run this file to open the window. The view is always the same three days,
after the normal 96-hour warmup, at 1.5-second resolution. Failure checks cover only
those days. Green bands retain the previously tested eight-day ranges.
Gray curves show the digitized Licinio day repeated, aligned to 08:00.
The companion FFT window updates with the sliders, using seven-minute samples.
"""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import MODEL

HERE = Path(__file__).resolve().parent
REFERENCE = MODEL.REFERENCE
FIXED_KEYS = frozenset(key for key, value in REFERENCE.items() if value == 1)
VIEW_START = MODEL.WARMUP_H
FFT_STEP_H = 7 / 60
FFT_HOURS = np.arange(0., 72., FFT_STEP_H)  # Three days starting at 08:00.
FFT_WINDOWS = {
    "full": (0.5, 36.), "1 h": (0.9, 1.45), "2 h": (1.65, 2.3),
    "4-6 h": (3.5, 6.5), "8-12 h": (7., 13.),
}


def load_licinio():
    """Load the digitized Figure 2 day; its time origin is 08:00."""
    return np.loadtxt(HERE.parent / "meeting12/licinio.csv", delimiter=",",
                      skiprows=1, unpack=True)


def fft_power(values):
    """Match meeting12/compare.py: linear detrend, then |rFFT|^2 / N.

    The zero-frequency and final rFFT bins are omitted, as in that comparison.
    This is raw FFT power, without peak normalization or PSD scaling.
    """
    values = np.asarray(values, dtype=float)
    design = np.column_stack((np.arange(1, len(values) + 1) / len(values),
                              np.ones(len(values))))
    detrended = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
    coefficients = np.fft.rfft(detrended)
    return np.abs(coefficients[1:-1]) ** 2 / len(values)


def fft_spectrum(values):
    frequency = np.fft.rfftfreq(len(values), d=FFT_STEP_H)[1:-1]
    return 1 / frequency, fft_power(values)

# Saved eight-day, one-parameter-at-a-time bands; not joint bounds.
RANGES = {
    "constant_inflow": {'low': 19.106250000000003, 'high': 20.80625},
    "production_scale": {'low': 245.78125, 'high': 253.515625},
    "activation_exponent": {'low': 3.1537500000000005, 'high': 5.59},
    "inhibition_exponent": {'low': 3.55375, 'high': 4.43625},
    "n_e2_midpoint": {'low': 49.46875, 'high': 50.640625},
    "central_hill_exponent": {'low': 2.7750000000000004, 'high': 3.1350000000000002},
    "period_min": {'low': 60.876562500000006, 'high': 68.73750000000001},
    "period_max": {'low': 73.2421875, 'high': 76.40625},
    "scn_min": {'low': 38.03125, 'high': 55.046875},
    "scn_max": {'low': 243.046875, 'high': 258.359375},
    "sleep_multiplier": {'low': 1.513578125, 'high': 1.605815625},
    "scn_phase": {'low': 15.035625, 'high': 17.139375},
    "gnrh_half_life_min": {'low': 3.91, 'high': 4.09125},
    "lh_half_life_min": {'low': 28.565, 'high': 29.4621875},
    "e2_half_life_min": {'low': 22.2328125, 'high': 22.732031250000002},
}
GROUPS = {
    "E2 Constants": (
        ("constant_inflow", "Constant E2 inflow: basal ovarian E2 production",
         r"$dE_2/dt=\mathbf{b}+V R H_+(LH;K_A,a)-k_E E_2$"),
        ("production_scale", "LH-driven E2 production scale: maximal responsive production",
         r"$dE_2/dt=b+\mathbf{V}R H_+(LH;K_A,a)-k_E E_2$"),
    ),
    "Hypothalamus": (
        ("n_e2_midpoint", "E2-to-N midpoint (pg/mL): half-maximal N target",
         r"$H_-(E_2)=\mathbf{K_N}^{h}/(\mathbf{K_N}^{h}+E_2^h)$"),
        ("scn_min", "SCN minimum drive: lower circadian pulse-amplitude factor",
         r"$M_{SCN}=(S_{max}+\mathbf{S_{min}})/2+(S_{max}-\mathbf{S_{min}})\cos(2\pi(c-\theta_{clock})/24)/2$"),
        ("scn_max", "SCN maximum drive: upper circadian pulse-amplitude factor",
         r"$M_{SCN}=(\mathbf{S_{max}}+S_{min})/2+(\mathbf{S_{max}}-S_{min})\cos(2\pi(c-\theta_{clock})/24)/2$"),
    ),
    "Literature grounded": (
        ("period_min", "Fast period coefficient (min): asymptotic strong-N limit",
         r"$P_0(N)=P_{max}-(P_{max}-\mathbf{P_{min}})H_+(N;1,h)$"),
        ("period_max", "Slow period coefficient (min): weak-N limit",
         r"$P_0(N)=\mathbf{P_{max}}-(\mathbf{P_{max}}-P_{min})H_+(N;1,h)$"),
        ("sleep_multiplier", "Sleep period multiplier: slowing during 23:00-07:00",
         r"$P(N,t)=P_0(N)\,\mathbf{M_{sleep}(t)}$"),
        ("scn_phase", "SCN peak time (clock hour): 23.5 = 23:30; enter a clock hour or HH:MM",
         r"$M_{SCN}=(S_{max}+S_{min})/2+(S_{max}-S_{min})\cos(2\pi(c-\mathbf{\theta}_{clock})/24)/2$"),
    ),
    "Steepness": (
        ("central_hill_exponent", "Central Hill exponent: E2 feedback and both N pulse maps",
         r"$N^*=G_NH_-(E_2;K_N,\mathbf{h});\quad H_+(N;1,\mathbf{h})$"),
        ("activation_exponent", "LH activation exponent: steepness of stimulatory input",
         r"$H_+(LH)=LH^{\mathbf{a}}/(K_A^{\mathbf{a}}+LH^{\mathbf{a}})$"),
        ("inhibition_exponent", "LH inhibition exponent: steepness of the R target",
         r"$dR/dt=K_I^{\mathbf{m}}/(K_I^{\mathbf{m}}+LH^{\mathbf{m}})-R$"),
    ),
    "Half-lives": (
        ("gnrh_half_life_min", "GnRH half-life (min): GnRH clearance time",
         r"$dGnRH/dt=-[60\ln(2)/\mathbf{t_{1/2,G}}]GnRH$"),
        ("lh_half_life_min", "LH half-life (min): LH clearance time",
         r"$dLH/dt=G_L GnRH-[60\ln(2)/\mathbf{t_{1/2,L}}]LH$"),
        ("e2_half_life_min", "E2 half-life (min): E2 clearance time",
         r"$dE_2/dt=b+VRH_+(LH)-[60\ln(2)/\mathbf{t_{1/2,E}}]E_2$"),
    ),
    "Fixed at 1": (
        ("gnrh_pulse_jump", "GnRH pulse gain: pulse size before the N and SCN factors",
         r"$\Delta GnRH=\mathbf{A_G}H_+(N;1,h)M_{SCN}(t)$"),
        ("lh_gnrh_gain", "GnRH-to-LH gain: LH production driven by GnRH",
         r"$dLH/dt=\mathbf{G_L}GnRH-k_L LH$"),
        ("activation_midpoint", "LH activation midpoint: half-maximal stimulatory input",
         r"$H_+(LH)=LH^a/(\mathbf{K_A}^a+LH^a)$"),
        ("inhibition_midpoint", "LH inhibition midpoint: target for ovarian responsiveness R",
         r"$dR/dt=\mathbf{K_I}^{m}/(\mathbf{K_I}^{m}+LH^m)-R$"),
        ("n_gain", "N gain: maximal E2-inhibited N target",
         r"$dN/dt=\lambda_N[\mathbf{G_N}H_-(E_2;K_N,h)-N]$"),
        ("n_relaxation", "N adjustment rate (1/h): response speed of the N signal",
         r"$dN/dt=\mathbf{\lambda_N}[G_N H_-(E_2;K_N,h)-N]$"),
    ),
    "Robustness": (),
    "Failures": (),
}
PARAMETER_LABELS = {key: description.split(":")[0]
                    for rows in GROUPS.values() for key, description, _ in rows}
PARAMETER_LABELS["constant_inflow"] = "Basal E2 production"
PARAMETER_LABELS["n_gain"] = "N drive scale"
PARAMETER_LABELS["scn_phase"] = "SCN peak (clock time)"

DOMAIN = {key: (value * .25, value * 2) for key, value in REFERENCE.items()}
DOMAIN.update(
    activation_exponent=(1., 8.), inhibition_exponent=(1., 8.),
    period_min=(REFERENCE["period_min"] * .25, REFERENCE["period_max"]),
    period_max=(REFERENCE["period_min"], REFERENCE["period_max"] * 2),
    sleep_multiplier=(1., 3.14), scn_phase=(12.5, 20.5),
)

def load_pyplot():
    """Load Matplotlib from the selected Python environment."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("This program needs NumPy and Matplotlib. Install them in your VS Code Python environment with:\n"
                         "    python -m pip install numpy matplotlib") from exc
    return plt

def validate(parameters):
    if any(parameters[key] != 1 for key in FIXED_KEYS):
        raise ValueError("The six unit coefficients are fixed at 1 in this simplified explorer.")
    for key, value in parameters.items():
        low, high = DOMAIN[key]
        if not np.isfinite(value) or not low <= value <= high:
            if key == "scn_phase":
                raise ValueError(
                    f"SCN peak must be between {clock_time_label(low)} and "
                    f"{clock_time_label(high)} (across midnight).")
            raise ValueError(f"{key} must be between {low:g} and {high:g}.")
    if parameters["period_min"] > parameters["period_max"]:
        raise ValueError("Minimum pulse period exceeds the maximum.")
    if parameters["scn_min"] > parameters["scn_max"]:
        raise ValueError("Minimum SCN drive exceeds the maximum.")

def local_events(time_h, e2):
    peaks = np.flatnonzero((e2[1:-1] > e2[:-2]) & (e2[1:-1] >= e2[2:])) + 1
    troughs = np.flatnonzero((e2[1:-1] < e2[:-2]) & (e2[1:-1] <= e2[2:])) + 1
    events = []
    for peak in peaks:
        position = np.searchsorted(troughs, peak)
        if position == 0 or position == len(troughs):
            continue
        left, right = troughs[position - 1], troughs[position]
        rise, fall = float(e2[peak] - e2[left]), float(e2[peak] - e2[right])
        events.append(dict(left_h=float(time_h[left]), peak_h=float(time_h[peak]),
                           right_h=float(time_h[right]), trough_before=float(e2[left]),
                           peak=float(e2[peak]), trough_after=float(e2[right]),
                           rise=rise, fall=fall, excursion=min(rise, fall)))
    return events, peaks

def stats(v):
    return [float(np.min(v)),float(np.median(v)),float(np.max(v))] if len(v) else [0.,0.,0.]

E2_DAY_MIN, E2_DAY_MAX = 15., 30.
E2_NIGHT_MEDIAN_MIN, E2_NIGHT_MAX = 5., 15.
E2_DAY_NIGHT_RATIO = 4. / 3.
NIGHT_MIN_TOLERANCE = 1e-9
LH_DAY_MEDIAN_IPI_MIN, LH_DAY_MEDIAN_IPI_MAX = 60., 75.


def assess(trace):
    """Measure the displayed three days using the existing raw-pulse rules."""
    t = trace["time_h"]
    ee, _ = local_events(t, trace["E2"])
    le, li = local_events(t, trace["LH"])
    daily = []
    for j in range(3):
        s = VIEW_START + 1 + 24 * j  # Daytime is 08:00-20:00.
        d = [e for e in ee if s <= e['peak_h'] < s + 12 and s <= e['left_h'] and e['right_h'] < s + 12]
        n = [e for e in ee if s + 15 <= e['peak_h'] < s + 23]
        ld = [e for e in le if s <= e['peak_h'] < s + 12]
        ln = [e for e in le if s + 15 <= e['peak_h'] < s + 23]
        lp = li[(t[li] >= s) & (t[li] < s + 24)]
        cycles = []
        for a, b in zip(lp[:-1], lp[1:]):
            if not s <= t[a] < t[b] < s + 12:
                continue
            paired = [e['excursion'] for e in d if t[a] <= e['peak_h'] < t[b]]
            raw = [e for e in ee if t[a] <= e['peak_h'] < t[b]]
            if raw and not paired:
                continue
            cycles.append(max(paired, default=0.))
        daytime_peaks = t[li][(t[li] >= s) & (t[li] < s + 12)]
        row = dict(E_day=stats([e['excursion'] for e in d]),
                   E_sleep=stats([e['excursion'] for e in n]),
                   E_cycle_day=stats(cycles), L_day=stats([e['excursion'] for e in ld]),
                   L_sleep=stats([e['excursion'] for e in ln]),
                   sleep_count=int(((t[lp] >= s + 15) & (t[lp] < s + 23)).sum()),
                   LH_day_IPI=stats(np.diff(daytime_peaks) * 60))
        row['target'] = not daily_failures(row)
        daily.append(row)
    return dict(daily=daily)


def preview(parameters):
    # The extra hour supplies the same boundary context for pulse measurements.
    trace = MODEL.simulate_retained(parameters, seconds=1.5, warm=VIEW_START, output=73.)
    return dict(parameters=parameters, trace=trace, audit=assess(trace))


def number(value, bounds=()):
    value = float(value)
    text = f"{value:.6g}"
    if any(float(text) == b and value != b for b in bounds):
        text = f"{value:.12g}"
    return text

def daily_failures(day):
    """Return measured/required values for each violated leaf requirement."""
    found = []

    def add(key, label, measured, required):
        found.append(dict(key=key, criterion=label, measured=measured, required=required))

    event = float(day["E_day"][1])
    cycle = float(day["E_cycle_day"][1])
    night = float(day["E_sleep"][1])
    maximum = float(day["E_sleep"][2])
    for key, label, value in (
        ("day_event_size", "Daytime E2 median (events)", event),
        ("day_cycle_size", "Daytime E2 median (LH cycles)", cycle),
    ):
        if not E2_DAY_MIN <= value <= E2_DAY_MAX:
            add(key, label, number(value, (E2_DAY_MIN, E2_DAY_MAX)) + " pg/mL",
                f"{E2_DAY_MIN:g}-{E2_DAY_MAX:g} pg/mL")
    if night < E2_NIGHT_MEDIAN_MIN - NIGHT_MIN_TOLERANCE:
        add("night_minimum", "Nighttime E2 median minimum",
            number(night, (E2_NIGHT_MEDIAN_MIN,)) + " pg/mL",
            f">= {E2_NIGHT_MEDIAN_MIN:g} pg/mL")
    if maximum > E2_NIGHT_MAX:
        add("night_maximum", "Largest nighttime E2 excursion",
            number(maximum, (E2_NIGHT_MAX,)) + " pg/mL", f"<= {E2_NIGHT_MAX:g} pg/mL")
    for key, label, value in (
        ("event_ratio", "Day/night E2 ratio (events)", event),
        ("cycle_ratio", "Day/night E2 ratio (LH cycles)", cycle),
    ):
        if not (night > 0 and value >= E2_DAY_NIGHT_RATIO * night):
            ratio = number(value / night, (E2_DAY_NIGHT_RATIO,)) + "x" if night > 0 else "undefined"
            add(key, label, f"day {number(value)}; night {number(night)} ({ratio})",
                "day >= 4/3 night; night > 0")
    count = int(day["sleep_count"])
    if not 3 <= count <= 5:
        add("LH_count", "Nighttime LH peak count", str(count) + " peaks", "3-5 peaks")
    ld, ln = float(day["L_day"][1]), float(day["L_sleep"][1])
    if not ln > ld > 0:
        add("LH_order", "Night/day LH median excursion",
            f"night {number(ln)}; day {number(ld)} IU/L", "night > day > 0")
    ipi = day.get("LH_day_IPI", (0., 0., 0.))
    median = float(ipi[1])
    if not (median > 0 and LH_DAY_MEDIAN_IPI_MIN <= median <= LH_DAY_MEDIAN_IPI_MAX):
        measured = (number(median, (LH_DAY_MEDIAN_IPI_MIN, LH_DAY_MEDIAN_IPI_MAX)) + " min"
                    if median > 0 else "No daytime interval available")
        add("day_period", "Daytime LH median interval", measured,
            f"{LH_DAY_MEDIAN_IPI_MIN:g}-{LH_DAY_MEDIAN_IPI_MAX:g} min")
    return found

def failure_report(result):
    if result is None:
        return dict(kind="pending", records=[], failed_checks=0, total_checks=0, criteria=[])
    daily = result["audit"]["daily"]
    records = []
    for j, day in enumerate(daily):
        for reason in daily_failures(day):
            records.append(dict(reason, model_day=j + 1, seconds=1.5, warm=VIEW_START))
    return dict(kind="preview", records=records,
                failed_checks=sum(not day["target"] for day in daily), total_checks=3,
                criteria=list(dict.fromkeys(r["criterion"] for r in records)))


def percentage(value, reference):
    return 100 * (value / reference - 1)

def display_value(key, value):
    """Show the SCN phase in clock hours while retaining the solver's origin."""
    return (value + 7) % 24 if key == "scn_phase" else value

def parse_input(key, text):
    """Translate decimal clock hours or HH:MM to the stored SCN phase."""
    if key != "scn_phase":
        return float(text)
    if ":" in text:
        fields = text.strip().split(":")
        if len(fields) != 2:
            raise ValueError("Use HH:MM or decimal clock hours.")
        hour, minute = (int(field) for field in fields)
        if not (0 <= hour < 24 and 0 <= minute < 60):
            raise ValueError("Clock time must be between 00:00 and 23:59.")
        clock = hour + minute / 60
    else:
        clock = float(text)
        if not np.isfinite(clock) or not 0 <= clock < 24:
            raise ValueError("Clock hour must be at least 0 and less than 24.")
    return (clock - 7) % 24

def change_text(key, value, reference):
    if key == "scn_phase":
        return f"{60 * (value - reference):+.2f} min"
    return f"{percentage(value, reference):+.2f}%"

def band_text(key):
    band = RANGES[key]
    low, high = band["low"], band["high"]
    span = (f"{clock_time_label(low)}-{clock_time_label(high)}" if key == "scn_phase"
            else f"{low:.7g}-{high:.7g}")
    return (f"OAT pass: {span} "
            f"[{change_text(key, low, REFERENCE[key])}, {change_text(key, high, REFERENCE[key])}]")


def clock_time_label(hours, _position=None):
    """Format the displayed time axis as a 24-hour clock, starting at 07:00."""
    minutes = int(round((hours + 7) * 60)) % (24 * 60)
    return f"{minutes // 60:02d}:{minutes % 60:02d}"

class ParameterExplorer:
    def __init__(self, plt):
        from matplotlib.widgets import Button, CheckButtons
        from matplotlib.ticker import FuncFormatter, MultipleLocator
        self.parameters = dict(REFERENCE)
        self.group = next(iter(GROUPS))
        self.fine_scale = True
        self.one_at_time = True
        self.result = None
        self.row_artists = []
        self.row_axes = []
        self.row_widgets = []
        self.pending = None
        self.future = None
        self.closed = False
        self.failure_page = 0
        self.failure_table = None
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="hpg_model")
        self.figure, self.axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
        self.figure.canvas.manager.set_window_title("Simplified HPG parameter sweep - six unit coefficients fixed")
        self.figure.subplots_adjust(left=.08, right=.97, top=.845, bottom=.535, hspace=.48)
        self.plot_lines, self.plot_peaks = [], []
        for axis, label, color in zip(self.axes, ("E2 (pg/mL)", "LH (IU/L)"), ("tab:red", "tab:blue")):
            axis.grid(alpha=.2)
            for start in (16, 40, 64):
                axis.axvspan(start, start + 8, color="0.93")
            self.plot_lines.append(axis.plot([], [], color=color, linewidth=1.2,
                                             zorder=3, label="Live model")[0])
            self.plot_peaks.append(axis.plot([], [], ".", color="black", markersize=4,
                                             zorder=4)[0])
            axis.set_ylabel(label)
        day_hours, day_e2 = load_licinio()
        # Preserve the digitized sample times and include both plot boundaries.
        licinio_hours = np.unique(np.r_[0., 72., *[
            1 + day_hours + 24 * day for day in range(-1, 3)]])
        licinio_hours = licinio_hours[(licinio_hours >= 0) & (licinio_hours <= 72)]
        self.licinio_line = self.axes[0].plot(
            licinio_hours, np.interp(licinio_hours - 1, day_hours, day_e2, period=24),
            color="0.65", linewidth=2, zorder=2,
            label="Licinio (digitized day repeated)")[0]
        self.axes[0].legend(loc="upper right", fontsize=7, ncol=2)
        self.build_fft(plt, np.interp(FFT_HOURS, day_hours, day_e2, period=24))
        self.axes[1].set_xlabel("Clock time (gray = sleep 23:00-07:00)", fontsize=9)
        self.axes[1].set_xlim(0, 72)
        self.axes[1].xaxis.set_major_locator(MultipleLocator(4))
        self.axes[1].xaxis.set_major_formatter(FuncFormatter(clock_time_label))
        self.axes[0].set_ylim(0, 60)
        self.axes[0].set_yticks(np.arange(0, 61, 10))
        self.figure.suptitle("Follicular timing: integer reference; six coefficients fixed at 1", y=.966)
        self.status = self.figure.text(.5, .918, "", ha="center", fontsize=9)
        self.summary = self.figure.text(.5, .890, "", ha="center", fontsize=9, color="0.3")
        self.figure.text(.5, .477,
            "Fixed R = 1 h. Green = saved 8-day OAT bands; dotted = reference. Lo/Hi select endpoints. Live checks cover only the displayed 3 days.",
            ha="center", fontsize=8.2)
        self.group_buttons = {}
        tab_step = .62 / len(GROUPS)
        for j, name in enumerate(GROUPS):
            button = Button(self.figure.add_axes([.025 + j * tab_step, .427, tab_step - .006, .032]), name.replace("Literature grounded", "Literature\ngrounded"))
            button.label.set_fontsize(9)
            button.on_clicked(lambda _, selected=name: self.change_group(selected))
            self.group_buttons[name] = button
        self.scale_button = Button(self.figure.add_axes([.655, .427, .155, .032]), "Scale: fine around target")
        self.scale_button.label.set_fontsize(8.5)
        self.scale_button.on_clicked(self.toggle_scale)
        self.oat_check = CheckButtons(self.figure.add_axes([.831, .420, .16, .045]), ["One parameter at a time"], [True])
        self.oat_check.labels[0].set_fontsize(8.3)
        self.oat_check.on_clicked(self.toggle_oat)
        self.reset_button = Button(self.figure.add_axes([.43, .018, .14, .035]), "Reset")
        self.reset_button.on_clicked(self.reset)
        self.figure.text(.025, .035, "Fixed view: first 3 days after the 96h warmup", va="center", fontsize=8.5)
        self.timer = self.figure.canvas.new_timer(interval=50)
        self.timer.add_callback(self.poll)
        self.figure.canvas.mpl_connect("close_event", self.close)
        self.build_controls()
        self.queue()
        self.timer.start()

    def build_fft(self, plt, licinio):
        """Keep the full spectrum and four zooms readable beside the controls."""
        self.fft_figure, self.fft_axes = plt.subplot_mosaic(
            [["full", "full"], ["1 h", "2 h"], ["4-6 h", "8-12 h"]],
            figsize=(11, 8), constrained_layout=True)
        self.fft_figure.canvas.manager.set_window_title("Live E2 FFT power vs Licinio")
        self.fft_figure.suptitle(
            "E2 FFT power: live model vs Licinio\n"
            "3 days from 08:00; 7-minute samples; linear detrend; power = |FFT|² / N",
            fontsize=11)
        period, power = fft_spectrum(licinio)
        self.fft_lines = {}
        for name, (low, high) in FFT_WINDOWS.items():
            axis = self.fft_axes[name]
            keep = (period >= low) & (period <= high)
            axis.plot(period[keep], power[keep], color="0.65", linewidth=2,
                      zorder=1, label="Licinio (digitized day repeated)")
            self.fft_lines[name] = axis.plot(
                [], [], ".-", color="tab:red", linewidth=1.2, markersize=3,
                zorder=2, label="Live model")[0]
            axis.set(xlim=(low, high), ylim=(0, None), xlabel="Period (hours)",
                     ylabel="FFT power", title="Full spectrum" if name == "full"
                     else f"Zoom: {name}")
            axis.grid(alpha=.2)
        axis = self.fft_axes["full"]
        axis.set_xscale("log")
        ticks = [0.5, 1, 2, 4, 8, 12, 24, 36]
        axis.set_xticks(ticks, labels=[str(tick) for tick in ticks])
        axis.legend(fontsize=8)

    def update_fft(self, trace):
        # The retained 73rd hour supplies 07:00-08:00 on the final morning.
        values = np.interp(VIEW_START + 1 + FFT_HOURS, trace["time_h"], trace["E2"])
        period, power = fft_spectrum(values)
        for name, (low, high) in FFT_WINDOWS.items():
            keep = (period >= low) & (period <= high)
            self.fft_lines[name].set_data(period[keep], power[keep])
            axis = self.fft_axes[name]
            axis.relim()
            axis.autoscale(enable=True, axis="y")
            axis.set_ylim(bottom=0, top=None, auto=True)
        self.fft_figure.canvas.draw_idle()

    def slider_limits(self, key):
        low, high = DOMAIN[key]
        if self.fine_scale:
            band = RANGES[key]
            padding = .4 * (band["high"] - band["low"])
            low = max(low, band["low"] - padding)
            high = min(high, band["high"] + padding)
        return min(low, self.parameters[key]), max(high, self.parameters[key])

    def build_controls(self):
        from matplotlib.widgets import Button, Slider, TextBox
        for widget in self.row_widgets:
            widget.disconnect_events()
        for axis in self.row_axes:
            axis.remove()
        for artist in self.row_artists:
            artist.remove()
        self.row_artists, self.row_axes, self.row_widgets = [], [], []
        self.sliders, self.textboxes, self.boundary_buttons = {}, {}, {}
        self.delta_labels = {}
        self.robustness_table = None
        self.failure_table = None
        if self.group == "Robustness":
            self.build_robustness_table()
            for name, button in self.group_buttons.items():
                button.ax.set_facecolor("0.82" if name == self.group else "0.94")
            self.figure.canvas.draw_idle()
            return
        if self.group == "Failures":
            self.build_failure_details()
            for name, button in self.group_buttons.items():
                button.ax.set_facecolor("0.82" if name == self.group else "0.94")
            self.figure.canvas.draw_idle()
            return
        if self.group == "Literature grounded":
            self.row_artists.append(self.figure.text(.025, .405,
                "Literature-informed timing; coefficients and the exact 23:30 SCN peak remain model choices.",
                fontsize=8, color="0.3"))
        rows = GROUPS[self.group]
        row_step = min(.056, .28 / max(1, len(rows) - 1))
        for j, (key, description, equation) in enumerate(rows):
            y = .363 - j * row_step
            label = self.figure.text(.025, y + .011, f"{description}\n{equation}", va="center", fontsize=8.4)
            self.row_artists.append(label)
            if key in FIXED_KEYS:
                self.row_artists.append(self.figure.text(.77, y + .002, "1  (fixed)",
                                                        va="center", fontsize=11, color="0.3"))
                continue
            low, high = self.slider_limits(key)
            axis = self.figure.add_axes([.70, y, .185, .023])
            slider = Slider(axis, "", low, high, valinit=self.parameters[key], valfmt="%0.6g", initcolor="none")
            slider.valtext.set_visible(False)
            slider.on_changed(lambda value, field=key: self.set_value(field, float(value)))
            band = RANGES[key]
            axis.axvspan(band["low"], band["high"], ymin=.12, ymax=.88, color="tab:green", alpha=.22, zorder=2)
            axis.axvline(band["low"], color="tab:green", linewidth=1.8, zorder=4)
            axis.axvline(band["high"], color="tab:green", linewidth=1.8, zorder=4)
            axis.axvline(REFERENCE[key], color="0.25", linestyle=":", linewidth=1, zorder=4)
            self.row_artists.append(self.figure.text(.70, y + .028,
                f"Low {change_text(key, band['low'], REFERENCE[key])}",
                fontsize=7.5, color="darkgreen", va="bottom"))
            self.row_artists.append(self.figure.text(.885, y + .028,
                f"High {change_text(key, band['high'], REFERENCE[key])}",
                fontsize=7.5, color="darkgreen", ha="right", va="bottom"))
            self.sliders[key] = slider
            number_axis = self.figure.add_axes([.918, y - .001, .072, .027])
            textbox = TextBox(number_axis, "", initial=f"{display_value(key, self.parameters[key]):.11g}")
            textbox.text_disp.set_fontsize(8.5)
            textbox.on_submit(lambda value, field=key: self.submit_number(field, value))
            self.textboxes[key] = textbox
            self.row_axes += [axis, number_axis]
            self.row_widgets += [slider, textbox]
            for side, x, title in (("low", .614, "Lo"), ("high", .657, "Hi")):
                button = Button(self.figure.add_axes([x, y - .001, .035, .027]), title)
                button.label.set_fontsize(8)
                button.on_clicked(lambda _, field=key, direction=side: self.set_value(field, RANGES[field][direction]))
                self.row_axes.append(button.ax)
                self.row_widgets.append(button)
                self.boundary_buttons[key, side] = button
            self.row_artists.append(self.figure.text(.025, y - .017, band_text(key), fontsize=7.9, color="0.25"))
            delta = self.figure.text(.99, y - .017, "", ha="right", fontsize=8, color="0.25")
            self.row_artists.append(delta)
            self.delta_labels[key] = delta
        for name, button in self.group_buttons.items():
            button.ax.set_facecolor("0.82" if name == self.group else "0.94")
        self.sync_controls()
        self.figure.canvas.draw_idle()

    def build_robustness_table(self):
        """Complete reference table; the slider tabs use these same endpoints."""
        records = []
        for key, value in REFERENCE.items():
            if key in RANGES:
                band = RANGES[key]
                span = f"{band['low']:.6g} - {band['high']:.6g}"
                change = f"{percentage(band['low'], value):+.2f}% to {percentage(band['high'], value):+.2f}%"
                if key == "scn_phase":
                    span = f"{clock_time_label(band['low'])} - {clock_time_label(band['high'])}"
                    change = (f"{change_text(key, band['low'], value)} to "
                              f"{change_text(key, band['high'], value)}")
            else:
                span, change = "Fixed", "Not swept"
            default = clock_time_label(value) if key == "scn_phase" else f"{value:g}"
            records.append((key, [PARAMETER_LABELS[key], default, span, change]))
        records.append(("fixed_R_tau_h", ["R response time (h)", "1", "Fixed", "Not swept"]))
        axis = self.figure.add_axes([.025, .063, .965, .330])
        axis.set_axis_off()
        self.row_axes.append(axis)
        self.robustness_table = axis.table(
            cellText=[values for _, values in records],
            colLabels=["Parameter", "Default", "Approx. passing range", "Change from default"],
            colWidths=[.36, .10, .28, .26], cellLoc="left", colLoc="left", bbox=[0, 0, 1, 1])
        self.robustness_table.auto_set_font_size(False)
        self.robustness_table.set_fontsize(7.4)
        for (row, col), cell in self.robustness_table.get_celld().items():
            cell.set_edgecolor("0.84")
            cell.set_linewidth(.35)
            if row == 0:
                cell.set_facecolor("0.25")
                cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
            elif records[row - 1][0] not in RANGES:
                cell.set_facecolor("0.94")
            elif col >= 2:
                cell.set_facecolor("#e5f3e5")
        self.row_artists.append(self.figure.text(.025, .405,
            "8-day reference ranges: one parameter at a time; all others at default. Green = passing bands. Fixed values were not swept.",
            fontsize=8, color="0.25"))

    def current_failure_report(self):
        return failure_report(self.result)


    def build_failure_details(self):
        from matplotlib.widgets import Button
        report = self.current_failure_report()
        records = report["records"]
        self.failure_table = None
        if report["kind"] == "pending":
            self.row_artists.append(self.figure.text(.025, .385,
                "Waiting for measurements of the current settings.", fontsize=10))
            return
        header = (f"Displayed 3 days: {report['failed_checks']}/3 day checks failed; "
                  f"{len(records)} criterion violations.")
        self.row_artists.append(self.figure.text(.025, .405, header, fontsize=9,
                                                color="firebrick" if records else "darkgreen"))
        note = "Checks cover only the displayed three days; the saved bands were calculated separately."
        if records:
            note = "Failed criteria: " + "; ".join(report["criteria"])
        self.row_artists.append(self.figure.text(.025, .382, note, fontsize=7.5, color="0.3"))
        if not records:
            self.row_artists.append(self.figure.text(.025, .335,
                "All currently evaluated criteria pass.", fontsize=11, color="darkgreen"))
            return
        page_size = 8
        page_count = (len(records) + page_size - 1) // page_size
        self.failure_page = min(max(0, self.failure_page), page_count - 1)
        selected = records[self.failure_page * page_size:(self.failure_page + 1) * page_size]
        cells = [[r["criterion"], r["measured"], r["required"], str(r["model_day"]),
                  f"{r['seconds']:g}s", f"{r['warm']:g}h"] for r in selected]
        axis = self.figure.add_axes([.025, .112, .965, .252])
        axis.set_axis_off()
        self.row_axes.append(axis)
        self.failure_table = axis.table(
            cellText=cells, colLabels=["Failed criterion", "Measured", "Required", "Day*", "Step", "Warmup"],
            colWidths=[.30, .24, .25, .06, .06, .09], cellLoc="left", colLoc="left", bbox=[0, 0, 1, 1])
        self.failure_table.auto_set_font_size(False)
        self.failure_table.set_fontsize(7.9)
        for (row, col), cell in self.failure_table.get_celld().items():
            cell.set_edgecolor("0.84")
            cell.set_linewidth(.35)
            cell.set_facecolor("0.25" if row == 0 else "#fff0ed")
            if row == 0:
                cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
        self.row_artists.append(self.figure.text(.025, .100,
            "*Day after the 96h warmup. Only the first three days are checked.", fontsize=7.6, color="0.3"))
        for delta, x, title in ((-1, .025, "Prev"), (1, .11, "Next")):
            button = Button(self.figure.add_axes([x, .064, .075, .029]), title)
            button.on_clicked(lambda _, offset=delta: self.shift_failure_page(offset))
            button.set_active((self.failure_page > 0) if delta < 0 else (self.failure_page + 1 < page_count))
            self.row_axes.append(button.ax)
            self.row_widgets.append(button)
        self.row_artists.append(self.figure.text(.205, .078,
            f"Page {self.failure_page + 1}/{page_count}; all failed criteria are included.",
            va="center", fontsize=8))

    def shift_failure_page(self, delta):
        self.failure_page += delta
        self.build_controls()

    def sync_controls(self):
        for key, slider in self.sliders.items():
            value = self.parameters[key]
            if slider.val != value:
                slider.eventson = False
                slider.set_val(value)
                slider.eventson = True
            textbox = self.textboxes[key]
            text = f"{display_value(key, value):.11g}"
            if textbox.text != text:
                # TextBox.set_val forces an entire synchronous canvas redraw.
                # Update its text artist and let the single queued draw render
                # all changed controls together instead.
                textbox.text_disp.set_text(text)
                textbox.cursor_index = min(textbox.cursor_index, len(text))
                textbox.cursor.set_visible(False)
            self.delta_labels[key].set_text(f"change {change_text(key, value, REFERENCE[key])}")

    def set_value(self, key, value):
        if key in FIXED_KEYS:
            self.status.set_text(f"{key} is fixed at 1 in the simplified model.")
            self.figure.canvas.draw_idle()
            return
        candidate = dict(REFERENCE if self.one_at_time else self.parameters)
        candidate[key] = value
        try:
            validate(candidate)
        except ValueError as exc:
            self.status.set_text(str(exc))
            self.status.set_color("firebrick")
            self.sync_controls()
            self.figure.canvas.draw_idle()
            return
        self.parameters = candidate
        # A precise entry may be outside a fine slider's current view.
        if any(not slider.valmin <= candidate[field] <= slider.valmax for field, slider in self.sliders.items()):
            self.build_controls()
        else:
            self.sync_controls()
        self.queue()

    def submit_number(self, key, text):
        try:
            value = parse_input(key, text)
        except ValueError:
            message = ("Enter SCN clock time as HH:MM or a decimal hour (23:30 or 23.5)."
                       if key == "scn_phase" else f"Enter a numeric value for {key}.")
            self.status.set_text(message)
            self.status.set_color("firebrick")
            self.sync_controls()
            self.figure.canvas.draw_idle()
        else:
            self.set_value(key, value)

    def change_group(self, group):
        self.group = group
        self.build_controls()

    def toggle_scale(self, _=None):
        self.fine_scale = not self.fine_scale
        self.scale_button.label.set_text("Scale: fine around target" if self.fine_scale else "Scale: full range")
        self.build_controls()

    def toggle_oat(self, _=None):
        self.one_at_time = bool(self.oat_check.get_status()[0])
        if self.one_at_time:
            self.reset()
        else:
            self.refresh_status()

    def queue(self):
        # New slider settings replace any pending update.
        self.pending = dict(self.parameters)
        self.result = None
        self.refresh_status()
        if self.group == "Failures":
            self.build_controls()

    def poll(self):
        if self.closed:
            return
        if self.future is not None and self.future.done():
            future, self.future = self.future, None
            try:
                result = future.result()
            except Exception as exc:
                if self.running_parameters == self.parameters:
                    self.status.set_text(f"Could not update graphs: {exc}")
                    self.status.set_color("firebrick")
                    self.figure.canvas.draw_idle()
            else:
                if result["parameters"] == self.parameters:
                    self.apply_result(result)
        if self.future is None and self.pending is not None:
            self.running_parameters, self.pending = self.pending, None
            self.future = self.executor.submit(preview, self.running_parameters)

    def reset(self, _=None):
        self.parameters = dict(REFERENCE)
        self.build_controls()
        self.queue()

    def apply_result(self, result):
        self.result = result
        trace = result["trace"]
        t = trace["time_h"]
        first, last = VIEW_START, VIEW_START + 72
        mask = (t >= first) & (t <= last)
        times = t[mask]
        for index, (axis, key) in enumerate(zip(self.axes, ("E2", "LH"))):
            values = trace[key][mask]
            peaks = np.flatnonzero((values[1:-1] > values[:-2]) & (values[1:-1] >= values[2:])) + 1
            troughs = np.flatnonzero((values[1:-1] < values[:-2]) & (values[1:-1] <= values[2:])) + 1
            stride = 40  # One-minute samples at 1.5 s, plus every extremum.
            keep = np.unique(np.r_[np.arange(0, len(times), stride), len(times) - 1, peaks, troughs])
            self.plot_lines[index].set_data(times[keep] - first, values[keep])
            self.plot_peaks[index].set_data(times[peaks] - first, values[peaks])
            axis.relim()
            axis.autoscale_view(scalex=False)
        self.update_fft(trace)
        daily = result["audit"]["daily"]
        e_day = [row["E_day"][1] for row in daily]
        e_sleep = [row["E_sleep"][1] for row in daily]
        e_cap = max(row["E_sleep"][2] for row in daily)
        self.axes[0].set_title(f"Displayed E2 excursions: day median {min(e_day):.2f}-{max(e_day):.2f}; sleep median {min(e_sleep):.2f}-{max(e_sleep):.2f}; largest sleep {e_cap:.2f} pg/mL", fontsize=10)
        nights = [row["sleep_count"] for row in daily]
        day_ipi = [row["LH_day_IPI"][1] for row in daily]
        self.axes[1].set_title("LH sleep peaks: " + " / ".join(map(str, nights))
                               + f"; daytime median interval {min(day_ipi):.2f}-{max(day_ipi):.2f} min (target: 60-75)", fontsize=10)
        self.refresh_status()
        if self.group == "Failures":
            self.failure_page = 0
            self.build_controls()

    def refresh_status(self):
        report = self.current_failure_report()
        if self.result is None:
            self.status.set_text("Updating graphs... the curve still shows the preceding settings.")
            self.status.set_color("0.3")
        else:
            failed_days = sorted({row["model_day"] for row in report["records"]})
            outcome = "FAIL (days " + ", ".join(map(str, failed_days)) + ")" if failed_days else "PASS"
            self.status.set_text(f"First 3-day target: {outcome}. Graphs and checks: 1.5 s.")
            self.status.set_color("firebrick" if failed_days else "darkgreen")
        if report["records"]:
            row = report["records"][0]
            self.summary.set_text(f"Day {row['model_day']}: {row['criterion']}: {row['measured']}; required {row['required']}. See Failures for all details.")
            self.summary.set_color("firebrick")
        else:
            self.summary.set_text("Green bands apply one parameter at a time, with all others at default."
                                  if self.one_at_time else
                                  "Joint changes: reference-based green bands do not predict their combined outcome.")
            self.summary.set_color("0.3")
        self.figure.canvas.draw_idle()

    def close(self, _=None):
        if not self.closed:
            self.closed = True
            self.timer.stop()
            self.executor.shutdown(wait=False, cancel_futures=True)
            # Closing the explorer also closes its companion spectrum window.
            import matplotlib.pyplot as plt
            plt.close(self.fft_figure)


if __name__ == "__main__":
    plt = load_pyplot()
    explorer = ParameterExplorer(plt)
    plt.show()
    explorer.close()
