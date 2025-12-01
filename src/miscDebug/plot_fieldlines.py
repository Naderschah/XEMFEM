#!/usr/bin/env python3

import sys
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
# --------------------------------------------------------------------
# Command-line argument handling
# --------------------------------------------------------------------
if len(sys.argv) < 2:
    print("Usage: python plot_paths.py <electron_paths_debug.csv>")
    sys.exit(1)

csv_path = Path(sys.argv[1])
if not csv_path.is_file():
    print(f"Error: file does not exist: {csv_path}")
    sys.exit(1)

output_path = csv_path.with_suffix(".png")   # same directory, same base name

# --------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------
paths = defaultdict(list)

with csv_path.open() as f:
    reader = csv.reader(f)
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        seed_id   = int(row[0])
        step      = int(row[1])
        r         = float(row[2])
        z         = float(row[3])
        exit_code = int(row[4])
        paths[seed_id].append((r, z, exit_code))

# --------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------
code_colors = {
    0: ("None",          "black"),
    1: ("HitCathode",    "blue"),
    2: ("HitLiquidGas",  "green"),
    3: ("HitWall",       "orange"),
    4: ("LeftVolume",    "magenta"),
    5: ("MaxSteps",      "red"),
    6: ("DegenerateDt",  "purple"),
    7: ("HitAxis",       "yellow"),
}

plt.figure()

max_r = 0.0
max_z = 0.0

for seed_id, pts in paths.items():
    r = np.array([p[0] for p in pts])
    z = np.array([p[1] for p in pts])
    exit_code = pts[-1][2]

    # track global extents (absolute values)
    if r.size > 0:
        max_r = max(max_r, float(np.nanmax(np.abs(r))))
    if z.size > 0:
        max_z = max(max_z, float(np.nanmax(np.abs(z))))

    label, color = code_colors.get(exit_code, ("Unknown", "gray"))
    if label != "HitCathode": continue

    plt.plot(r, z, '-', linewidth=0.8, color=color, alpha=0.5)
    plt.scatter(r[-1], z[-1], s=12, marker="x", color=color)

    if np.isclose((r - np.roll(r, 1))[1:-2], np.zeros(r.shape)[1:-2]).all() and np.isclose((z - np.roll(z, 1))[1:-2], np.zeros(z.shape)[1:-2]).all():
      print("Pathline Did not move in exit code", exit_code)
      plt.scatter([r[0]], [z[0]], s=8, marker='x')

ax = plt.gca()

# Switch to log scale if extents are large
# (assuming r,z > 0 in your geometry)
if max_r > 1.0:
    ax.set_xlim(0,2)
if max_z > 10.0:
    ax.set_ylim(-1.7,0.2)

plt.xlabel("r")
plt.ylabel("z")
plt.grid(True)
plt.title("Debug Electron Paths")


plt.tight_layout()
plt.savefig(output_path, dpi=180)
print(f"Saved plot to: {output_path}")

# If interactive viewing is wanted:
# plt.show()
