"""
Side-by-side with the goal: model E2 (delayed feedback + circadian envelope)
vs the measured estradiol (panel B of bodytempestradiol.png).
Uses MODEL.run with the physiological 70-80 min GnRH band and circ=True.
"""
import numpy as np
import matplotlib.pyplot as plt
import MODEL

DAY = MODEL.DAY

# model: delayed feedback (n=4, tau=2h) + circadian envelope, 70-80 min pulses
t, G, L, E = MODEL.run(4, 4, tau=120, t_end=8*DAY, circ=True)
day_start = 5*DAY + 8*60                       # pick a settled day, starting 08:00
m = (t >= day_start) & (t < day_start + DAY)
tod = (t[m] - day_start)/60.0 + 8.0            # hours, 08:00 -> 08:00 next day
print(f"model day: E2 {E[m].min():.1f}..{E[m].max():.1f} pg/mL")

# crop panel B (estradiol) from the goal figure
goal = plt.imread("bodytempestradiol.png")
h = goal.shape[0]
panelB = goal[int(h*0.46):, :]

fig, ax = plt.subplots(2, 1, figsize=(11, 7.5))
ax[0].imshow(panelB); ax[0].axis("off")
ax[0].set_title("GOAL: measured estradiol over one day (panel B, Bao et al.)", fontsize=11)
ax[1].plot(tod, E[m], lw=1.0, color="C4")
ax[1].set_xticks([8, 14, 20, 26, 32]); ax[1].set_xticklabels(["08:00","14:00","20:00","02:00","08:00"])
ax[1].set_ylabel("E2 (pg/mL)"); ax[1].set_xlabel("time of day")
ax[1].set_title("MODEL: GnRH->LH->E2, 70-80 min pulses, delayed feedback (n=4, tau=2h) + circadian",
                fontsize=11)
ax[1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_comparison.png", dpi=130)
print("wrote fig_comparison.png")
