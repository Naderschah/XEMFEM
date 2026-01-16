#!/usr/bin/env python3
# FIXME Version mismatch inside container of numpy and matplotlib
# python3 -m pip install --no-cache-dir --force-reinstall "numpy<2"
# pip install h5py
import sys
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import meshio
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

plot_TPC = True
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

def extract_bc_groups_from_med(mesh: meshio.Mesh, *, prefer_cell_groups=True):
    """
    Extract boundary-condition groups from a meshio MED mesh.

    Returns
    -------
    bc : dict
        Mapping: group_name -> dict with:
          - "points": (N, dim) float array (same as mesh.points)
          - "cells":  list of (cell_type, connectivity ndarray) for that group
          - "cell_block_indices": list of indices into mesh.cells that contributed
          - "family_ids": sorted list of MED family ids used for that group

    Notes
    -----
    meshio's MED reader stores:
      - per-cell family ids in mesh.cell_data["cell_tags"] (one array per cell block)
      - mapping family_id -> [group names] in mesh.cell_tags
      - (optionally) per-point family ids in mesh.point_data["point_tags"]
        with mapping in mesh.point_tags
    """
    if prefer_cell_groups:
        tags_map = getattr(mesh, "cell_tags", {}) or {}
        data_key = "cell_tags"
        data = mesh.cell_data.get(data_key, None)
    else:
        tags_map = getattr(mesh, "point_tags", {}) or {}
        data_key = "point_tags"
        data = mesh.point_data.get(data_key, None)

    if not tags_map:
        raise ValueError("No MED family-name mapping found on mesh.(cell_tags|point_tags).")
    if data is None:
        raise ValueError(f"No tag array found in mesh.{('cell_data' if prefer_cell_groups else 'point_data')}['{data_key}'].")

    # Build reverse map: group_name -> set(family_ids)
    group_to_famids = {}
    for fam_id, names in tags_map.items():
        for name in names:
            group_to_famids.setdefault(name, set()).add(int(fam_id))

    bc = {}
    if prefer_cell_groups:
        # mesh.cell_data["cell_tags"] is a list aligned with mesh.cells (one tag array per cell block)
        if not isinstance(data, list) or len(data) != len(mesh.cells):
            raise ValueError("mesh.cell_data['cell_tags'] must be a list aligned with mesh.cells.")

        for group_name, fam_ids in group_to_famids.items():
            fam_ids = set(fam_ids)
            group_cells = []
            contributing_blocks = []
            used_fam_ids = set()

            for bidx, (cell_block, tag_arr) in enumerate(zip(mesh.cells, data)):
                tag_arr = np.asarray(tag_arr).astype(int)
                mask = np.isin(tag_arr, list(fam_ids))
                if not np.any(mask):
                    continue
                group_cells.append((cell_block.type, np.asarray(cell_block.data)[mask]))
                contributing_blocks.append(bidx)
                used_fam_ids.update(np.unique(tag_arr[mask]).tolist())

            if group_cells:
                bc[group_name] = {
                    "points": np.asarray(mesh.points),
                    "cells": group_cells,
                    "cell_block_indices": contributing_blocks,
                    "family_ids": sorted(used_fam_ids),
                }
    else:
        # Point groups: return point indices for each group
        tag_arr = np.asarray(data).astype(int)
        for group_name, fam_ids in group_to_famids.items():
            fam_ids = set(fam_ids)
            mask = np.isin(tag_arr, list(fam_ids))
            if not np.any(mask):
                continue
            bc[group_name] = {
                "points": np.asarray(mesh.points),
                "point_indices": np.where(mask)[0],
                "family_ids": sorted(np.unique(tag_arr[mask]).tolist()),
            }

    return bc


def bc_group_to_line_segments(mesh: meshio.Mesh, group_cells):
    """
    Convert a group's 1D cells (line/line3/...) into simple 2-node segments for plotting.

    Returns
    -------
    segs : (M, 2) int array of node indices
    """
    seg_list = []
    for cell_type, conn in group_cells:
        conn = np.asarray(conn).astype(int)

        if cell_type == "line":
            # (n,2)
            seg_list.append(conn[:, [0, 1]])
        elif cell_type == "line3":
            # (n,3): plot as 2 segments: 0-1 and 1-2
            seg_list.append(conn[:, [0, 1]])
            seg_list.append(conn[:, [1, 2]])
        else:
            # ignore non-1D cells in this helper
            continue

    if not seg_list:
        return np.zeros((0, 2), dtype=int)
    return np.vstack(seg_list)


# ---- example usage ----
mesh = meshio.read("/home/felix/work/geometry/mesh/mesh.med")

# Extract BC groups from cell families (typical for boundary conditions)
bc = extract_bc_groups_from_med(mesh, prefer_cell_groups=True)
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

    plt.plot(r, z, '-', linewidth=0.8, color=color, alpha=0.5)
    if len(paths.keys())==1:
        plt.scatter(r[:-2], z[:-2], s=1, marker="o", color=color)
    plt.scatter(r[-1], z[-1], s=3, marker="x", color=color)

    if np.isclose((r - np.roll(r, 1))[1:-2], np.zeros(r.shape)[1:-2]).all() and np.isclose((z - np.roll(z, 1))[1:-2], np.zeros(z.shape)[1:-2]).all():
      print("Pathline Did not move in exit code", exit_code)
      plt.scatter([r[0]], [z[0]], s=8, marker='x')

ax = plt.gca()

if plot_TPC:
    ax.set_aspect('equal', 'box')
    pts_mesh = mesh.points
    for g in bc.values():
        for cell_type, conn in g["cells"]:
            conn = np.asarray(conn, dtype=int)

            if cell_type == "line":
                for n0, n1 in conn:
                    plt.plot(
                        [pts_mesh[n0, 0], pts_mesh[n1, 0]],
                        [pts_mesh[n0, 1], pts_mesh[n1, 1]],
                        color="black",
                        linewidth=0.6,
                        alpha=0.9,
                    )

            elif cell_type == "line3":
                for n0, n1, n2 in conn:
                    plt.plot(
                        [pts_mesh[n0, 0], pts_mesh[n1, 0]],
                        [pts_mesh[n0, 1], pts_mesh[n1, 1]],
                        color="black",
                        linewidth=0.6,
                        alpha=0.9,
                    )
                    plt.plot(
                        [pts_mesh[n1, 0], pts_mesh[n2, 0]],
                        [pts_mesh[n1, 1], pts_mesh[n2, 1]],
                        color="black",
                        linewidth=0.6,
                        alpha=0.9,
                    )
else:
    ax.set_xlim(0, ax.get_xlim()[1])
# Switch to log scale if extents are large
# (assuming r,z > 0 in your geometry)
if not plot_TPC:
    if max_r > 1.0:
        ax.set_xlim(0,2)
    if max_z > 10.0:
        ax.set_ylim(-1.7,0.2)

plt.xlabel("r")
plt.ylabel("z")
plt.grid(True)
plt.title("Debug Electron Paths")

# Build legend entries: an "x" in the right color + label text
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0],
           marker='x', linestyle='None',
           color=color, markeredgecolor=color,
           markersize=4, markeredgewidth=1.5,
           label=label)
    for _, (label, color) in sorted(code_colors.items())
]
# Place legend outside to the right
ax.legend(
    handles=handles,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),   # (x,y) in axes fraction coords
    borderaxespad=0.0,
    frameon=True,
    title="Exit code"
)

# Make room on the right so the legend isn't cut off
plt.tight_layout(rect=[0, 0, 0.82, 1])

if plot_TPC:
    plt.savefig(output_path, dpi=1800)
else:
    plt.savefig(output_path, dpi=180)
print(f"Saved plot to: {output_path}")

# If interactive viewing is wanted:
# plt.show()
