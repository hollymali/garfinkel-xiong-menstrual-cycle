
"""Check the CSV points against the E2 trace in Licinio Figure 2."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

folder = Path(__file__).resolve().parent
hours, e2 = np.loadtxt(folder / "licinio.csv", delimiter=",", skiprows=1, unpack=True)
figure_image = plt.imread(folder / "licinio_e2_source.png")

# On this image, 0–24 hours spans x=287–1145 pixels.
# The vertical axis goes upward: 0–60 pg/mL spans y=826–220 pixels.
x = 287 + hours / 24 * (1145 - 287)
y = 826 - e2 / 60 * (826 - 220)

fig, ax = plt.subplots(figsize=(12, 8), layout="constrained")
ax.imshow(figure_image, cmap="gray")
ax.plot(x, y, ".", color="crimson", label="licinio.csv points")
# Keep the paper's printed axes visible. Put our calculated ticks above and to the right.
ax.set(xlim=(0, 1225), ylim=(950, 145), title="Licinio Figure 2 · axis calibration check")
ax.set_xticks(287 + np.arange(0, 25, 6) * (1145 - 287) / 24,
              labels=["08:00", "14:00", "20:00", "02:00", "08:00"])
ax.set_yticks(826 - np.arange(0, 61, 20) * (826 - 220) / 60,
              labels=["0", "20", "40", "60"])
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.tick_params(colors="tab:blue")
ax.grid(color="tab:blue", linestyle="--", linewidth=0.8, alpha=0.6)
ax.legend()
fig.savefig(folder / "licinio_digitization_check.png", dpi=150)
plt.show()


