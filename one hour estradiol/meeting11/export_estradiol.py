"""Export three days of default E2 after the model's 96-hour warmup.

time_h is hours after warmup: 0 = 07:00, 72 = 07:00 three days later.
Samples are 1.5 seconds apart; estradiol is in pg/mL.
Run this file with MODEL.py in the same folder.
"""

from pathlib import Path

import numpy as np
import MODEL

# Same reference parameters, warmup, step, and three-day span as the sweep.
trace = MODEL.simulate_endocrine(**MODEL.REFERENCE)
keep = trace["time_h"] >= MODEL.WARMUP_H
time_h = trace["time_h"][keep] - MODEL.WARMUP_H

output = Path(__file__).with_name("estradiol_3days.csv")
np.savetxt(output, np.column_stack((time_h, trace["E2"][keep])),
           delimiter=",", header="time_h,estradiol_pg_mL", comments="",
           fmt=("%.10f", "%.15g"))
print(f"Saved {output}")
