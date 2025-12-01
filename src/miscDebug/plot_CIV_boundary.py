#!/usr/bin/env python3

import sys
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------
# Command-line argument
# --------------------------------------------------------------------
if len(sys.argv) < 2:
    print("Usage: python save_civ_boundary_plot.py <civ_interface_boundary.csv>")
    sys.exit(1)

csv_path = Path(sys.argv[1])
if not csv_path.is_file():
    print(f"Error: file does not exist: {csv_path}")
    sys.exit(1)

out_path = csv_path.with_suffix(".png")

# --------------------------------------------------------------------
# Load interface curve
# --------------------------------------------------------------------
r_vals = []
z_vals = []

with csv_path.open() as f:
    reader = csv.reader(f)
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        r = float(row[0])
        z = float(row[1])
        r_vals.append(r)
        z_vals.append(z)

r_vals = np.array(r_vals)
z_vals = np.array(z_vals)

if r_vals.size == 0:
    print("No points found in interface CSV")
    sys.exit(1)

# --------------------------------------------------------------------
# Define TPC boundaries (hardcoded or read from config)
# --------------------------------------------------------------------
r_min = 0.0
r_max = 0.664
z_min = -1.5025
z_max =  0.004

# --------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------
plt.figure(figsize=(5.5, 7))

# TPC boundary rectangle (optional outline)
plt.plot([r_min, r_max, r_max, r_min, r_min],
         [z_min, z_min, z_max, z_max, z_min],
         'k--', linewidth=0.8, label='TPC boundary')

# Interface curve
plt.plot(r_vals, z_vals, color='red', linewidth=2.0, label='CIV boundary')
plt.scatter(r_vals, z_vals, color='red', s=12)

# --------------------------------------------------------------------
# Axis scaling: focus on the curve vertical extent
# --------------------------------------------------------------------
r_lo = max(r_min, np.min(r_vals) - 0.05)
r_hi = min(r_max, np.max(r_vals) + 0.05)

# Only zoom vertically to the z extent of the interface curve
z_lo = np.min(z_vals) - 0.05
z_hi = np.max(z_vals) + 0.05

plt.xlim(r_lo, r_hi)
plt.ylim(z_lo, z_hi)

plt.xlabel("r")
plt.ylabel("z")
plt.title("Charge-Insensitive Boundary (CIV Interface)")
plt.grid(True)
plt.legend(loc='best', fontsize=8)

plt.tight_layout()
plt.savefig(out_path, dpi=180)
print(f"Saved CIV boundary plot to: {out_path}")
