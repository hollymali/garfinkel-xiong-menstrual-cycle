"""
Add an ordered [N] number + the model's equations beneath each figure,
and re-save under a numbered filename.  All metrics are E2-based (no GnRH CV):
  E2_amplitude_var, E2_period_var = MAD/median of E2 peak heights / spacings.
Skips any figure whose source PNG isn't present (regenerate that .py first).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DPI = 130

entries = [
 ("bodytempestradiol.png", "[1] target.png",
  "[1]  TARGET  (experimental reference data -- no model)\n"
  "Panel A: body temperature over several days.\n"
  "Panel B: pulsatile estradiol across one follicular day."),

 ("fig_pi_cascade_amp.png", "[2] amp_only_sweep.png",
  "[2]  AMPLITUDE-ONLY SWEEP: 3-stage cascade + E2->amplitude (negative) feedback\n"
  "dGnRH/dt = -kG*GnRH   (+ a_n at each pulse, fixed interval T0)\n"
  "dFSH/dt  =  cF*GnRH - kF*FSH\n"
  "dE2/dt   =  cE*FSH  - kE*E2\n"
  "a_n = A / (1 + (E2/K)^n_amp)        (higher E2 -> smaller next pulse; sweep n_amp)\n"
  "title metric: orbit = period-# = how many E2 peaks until the pattern repeats\n"
  "(period-1 regular, period-2 alternating, >30 = high-period/irregular)"),

 ("fig_pi_cascade_freqpos.png", "[3] freq_only_sweep.png",
  "[3]  FREQUENCY-ONLY SWEEP: cascade + E2->frequency (positive) feedback (amplitude fixed)\n"
  "cascade  GnRH -> FSH -> E2  as in [2]   (pulse amplitude fixed = A)\n"
  "period(E2) = Pmin + (Pmax - Pmin)/(1 + (E2/K)^n_freq)   (higher E2 -> shorter period; sweep n_freq)"),

 ("fig_pi_cascade_combined.png", "[4] combined_freq_sweep.png",
  "[4]  COMBINED -- SWEEP FREQUENCY: both feedbacks; n_amp fixed = 8, sweep n_freq\n"
  "cascade  GnRH -> FSH -> E2  as in [2]\n"
  "a_n        = A / (1 + (E2/K_amp)^n_amp)                     (neg amplitude, fixed)\n"
  "period(E2) = Pmin + (Pmax - Pmin)/(1 + (E2/K_freq)^n_freq)  (pos frequency, swept)"),

 ("fig_ampsweep.png", "[5] combined_amp_sweep.png",
  "[5]  COMBINED -- SWEEP AMPLITUDE: both feedbacks; n_freq fixed = 8, sweep n_amp\n"
  "same model as [4], but n_freq held at 8 and n_amp swept.\n"
  "complement of [4]: E2 irregularity needs BOTH feedbacks steep (n >~ 8 each)."),

 ("fig_tol_compare.png", "[6] heatmaps.png",
  "[6]  ORBIT-PERIOD HEATMAP over (n_freq x n_amp): tolerance 3% (left) vs 1% (right)\n"
  "color = E2 orbit period = # of E2 peaks until the pattern repeats (discrete bands)\n"
  "period-1/2 = regular/alternating ; high-period (yellow) = ragged (both feedbacks steep)\n"
  "structure is robust to tolerance; 1% flips a few borderline cells to high-period. star = time series (14,16)."),

 ("fig_bifurcation.png", "[7] bifurcation.png",
  "[7]  DOUBLE BIFURCATION (E2): combined cascade, n_amp = n_freq = n swept\n"
  "top:    settled E2 peak HEIGHTS   (amplitude attractor)\n"
  "bottom: settled E2 peak-to-peak PERIODS in min   (timing attractor)\n"
  "for each: 1 dot = period-1, 2 dots = period-2, band = high-period"),

 ("fig_chaos_timeseries.png", "[8] chaos_timeseries.png",
  "[8]  n = 14 E2 TIME SERIES   (top: no circadian ; bottom: + circadian)\n"
  "Looks aperiodic, but IC-sensitivity test gives |dE2| -> 0\n"
  "=> high-period / quasiperiodic, NOT chaos."),

 ("fig_pi_cascade_circadian.png", "[9] cascade_circadian.png",
  "[9]  COMBINED CASCADE + CIRCADIAN IMPOSED LAST   (fullest panel-B result)\n"
  "cascade + both E2 feedbacks as in [4];   C(t)=0.5*[1+cos(2pi(t-14h)/24h)]  (SCN clock)\n"
  "A_circ(t) = A*(0.35 + 0.65*C(t))        (clock scales the pulse amplitude)\n"
  "a_n       = A_circ(t) / (1 + (E2/K_amp)^n_amp)"),
]

for src, dst, text in entries:
    if not os.path.exists(src):
        print(f"skip {dst} (source {src} not present)"); continue
    img = plt.imread(src)
    h, w = img.shape[:2]
    fig_w = w/DPI; img_h = h/DPI
    nlines = text.count("\n") + 1
    text_h = 0.30*nlines + 0.45
    fig = plt.figure(figsize=(fig_w, img_h + text_h), dpi=DPI)
    frac = img_h/(img_h + text_h)
    ax = fig.add_axes([0, 1-frac, 1, frac]); ax.imshow(img); ax.axis("off")
    fig.text(0.012, (1-frac)*0.5, text, ha="left", va="center",
             family="monospace", fontsize=10.5,
             bbox=dict(boxstyle="round,pad=0.5", fc="#f4f4f2", ec="#bbbbbb"))
    fig.savefig(dst, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"wrote {dst}")
