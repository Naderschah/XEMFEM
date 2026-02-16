#!/usr/bin/env python3
# python3 -m pip install --no-cache-dir --force-reinstall "numpy<2" && pip install h5py pyyaml meshio matplotlib

import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import meshio

import h5py
import yaml


# --------------------------------------------------------------------
# Path resolution
# --------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent

# Defaults relative to this script (not CWD)
DEFAULT_SIM_RESULTS = (SCRIPT_DIR / "../../sim_results").resolve()
DEFAULT_CFG = (SCRIPT_DIR / "../../geometry/config.yaml").resolve()
DEFAULT_MED = (SCRIPT_DIR / "../../geometry/mesh/mesh.med").resolve()

# Default per-run relative paths inside sim_results[/run_xxxx]
REL_CSV = Path("electron_paths_debug.csv")
REL_GRID = Path("interpolated/E_interpolated.h5")


# --------------------------------------------------------------------
# Defaults / CLI
# --------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=(
        "Plot traced electron paths; optionally overlay E-field from interpolated grid.\n"
        "If sim_results contains run_xxxx directories, iterates ONLY those and saves outputs within each run dir."
    )
)

parser.add_argument(
    "--sim-results",
    type=Path,
    default=DEFAULT_SIM_RESULTS,
    help=f"Base sim_results directory (default: {DEFAULT_SIM_RESULTS})",
)

# Optional: single-run mode (overrides run discovery)
parser.add_argument(
    "--csv",
    type=Path,
    default=None,
    help="CSV with traced paths. If set, runs only this CSV (no run_xxxx iteration).",
)
parser.add_argument(
    "-G",
    "--grid",
    type=Path,
    default=None,
    help="HDF5 grid file written by WriteGridSampleBinary. If set, used with --csv for single-run mode.",
)

parser.add_argument(
    "-c",
    "--config",
    type=Path,
    default=DEFAULT_CFG,
    help=f"Geometry config YAML containing optimize bounds (default: {DEFAULT_CFG})",
)
parser.add_argument(
    "--mesh-med",
    type=Path,
    default=DEFAULT_MED,
    help=f"MED mesh for TPC outline (default: {DEFAULT_MED})",
)

parser.add_argument(
    "-E",
    "--efield",
    action="store_true",
    help="Overlay E-field (|E| background + sparse quiver) inside bounds.",
)
grad_group = parser.add_mutually_exclusive_group()
grad_group.add_argument(
    "-dEdx",
    action="store_true",
    help="Overlay d|E|/dx (x is the first plotted coordinate, i.e. r).",
)
grad_group.add_argument(
    "-dEdy",
    action="store_true",
    help="Overlay d|E|/dy (y is the second plotted coordinate, i.e. z).",
)

parser.add_argument(
    "--no-tpc",
    action="store_true",
    help="Do not draw TPC mesh lines.",
)
parser.add_argument(
    "--dpi",
    type=int,
    default=None,
    help="Override DPI (default: 1800 with TPC, else 180).",
)
parser.add_argument(
    "--no-paths",
    action="store_true",
    help="Do not draw electron path lines/markers (useful to inspect E-field/TPC only).",
)

# Filter by exit code(s)
parser.add_argument(
    "--exit",
    type=str,
    default=None,
    help=(
        "Filter which paths to plot by exit code. "
        "Examples: --exit MaxSteps  |  --exit 5  |  --exit MaxSteps,HitWall  |  --exit 5,3"
    ),
)

# Plot bounds
parser.add_argument(
    "--xlim",
    type=float,
    nargs=2,
    metavar=("XMIN", "XMAX"),
    default=None,
    help="Override x-axis bounds (r). Example: --xlim 0 500",
)
parser.add_argument(
    "--ylim",
    type=float,
    nargs=2,
    metavar=("YMIN", "YMAX"),
    default=None,
    help="Override y-axis bounds (z). Example: --ylim -600 0",
)

args = parser.parse_args()
plot_TPC = (not args.no_tpc)


# --------------------------------------------------------------------
# sim_results/meta.txt loader (optional)
# --------------------------------------------------------------------
def load_meta(sim_results_dir: Path):
    """
    Parses sim_results/meta.txt with format like:

      SR3XENONnT

      /work/sim_results/run_0000:
        Preset = PMTsOn_RestGround

    Returns (sim_id: str|None, run_to_preset: dict[str,str]).
    """
    meta_path = sim_results_dir / "meta.txt"
    if not meta_path.is_file():
        return None, {}

    lines = meta_path.read_text().splitlines()

    # First non-empty line that isn't a section header is sim_id
    sim_id = None
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if s.endswith(":") and "run_" in s:
            continue
        sim_id = s
        break

    run_to_preset: dict[str, str] = {}
    current_run = None

    for raw in lines:
        s = raw.strip()
        if not s:
            continue

        # Section header: ".../run_0000:"
        if s.endswith(":") and "run_" in s:
            hdr = s[:-1]
            run_name = hdr.split("/")[-1]
            current_run = run_name if run_name.startswith("run_") else None
            continue

        # Content line: "Preset = X"
        if current_run is not None and "=" in s:
            left, right = s.split("=", 1)
            if left.strip() == "Preset":
                run_to_preset[current_run] = right.strip()

    return sim_id, run_to_preset


def build_plot_title(*, default_title: str, sim_id: str | None, preset: str | None):
    parts = []
    if sim_id:
        parts.append(sim_id)
    if preset:
        parts.append(preset)
    return " - ".join(parts) if parts else default_title


# Load meta once (applies to both batch and single-run)
sim_results_for_meta = args.sim_results.resolve()
meta_path = sim_results_for_meta / "meta.txt"
meta_sim_id, meta_run_to_preset = load_meta(sim_results_for_meta)

if not meta_path.is_file():
    print(f"[meta] No meta.txt found at: {meta_path} (using default titles)")
else:
    if meta_sim_id is None:
        print(f"[meta] meta.txt found but sim_id missing/parse failed: {meta_path}")
    if not meta_run_to_preset:
        print(f"[meta] meta.txt found but no run presets parsed: {meta_path}")


# --------------------------------------------------------------------
# YAML bounds (optimize box)
# --------------------------------------------------------------------
def _find_optimize_bounds(cfg_obj):
    required = {"r_min", "r_max", "z_min", "z_max"}

    def walk(x):
        if isinstance(x, dict):
            if required.issubset(x.keys()):
                return x
            for v in x.values():
                got = walk(v)
                if got is not None:
                    return got
        elif isinstance(x, list):
            for v in x:
                got = walk(v)
                if got is not None:
                    return got
        return None

    node = walk(cfg_obj)
    if node is None:
        return None
    return (
        float(node["r_min"]),
        float(node["r_max"]),
        float(node["z_min"]),
        float(node["z_max"]),
    )


def load_bounds_from_config(cfg_path: Path):
    if not cfg_path.is_file():
        raise FileNotFoundError(f"Config file does not exist: {cfg_path}")
    with cfg_path.open("r") as f:
        cfg = yaml.safe_load(f)
    b = _find_optimize_bounds(cfg)
    if b is None:
        raise KeyError("Could not find optimize bounds (r_min,r_max,z_min,z_max) anywhere in config YAML.")
    return b


def resolve_plot_bounds(args_, cfg_path: Path):
    """
    Returns (xmin, xmax, ymin, ymax).

    Priority:
      1) CLI --xlim/--ylim if provided
      2) If any overlay is enabled and CLI didn't specify that axis, use optimize bounds from config
      3) Otherwise leave as None (matplotlib autoscale)
    """
    xmin = xmax = ymin = ymax = None

    if args_.xlim is not None:
        xmin, xmax = float(args_.xlim[0]), float(args_.xlim[1])
    if args_.ylim is not None:
        ymin, ymax = float(args_.ylim[0]), float(args_.ylim[1])

    need_defaults_from_cfg = (args_.efield or args_.dEdx or args_.dEdy) and (
        (xmin is None and xmax is None) or (ymin is None and ymax is None)
    )
    if need_defaults_from_cfg:
        rmin, rmax, zmin, zmax = load_bounds_from_config(cfg_path)
        if xmin is None and xmax is None:
            xmin, xmax = rmin, rmax
        if ymin is None and ymax is None:
            ymin, ymax = zmin, zmax

    return xmin, xmax, ymin, ymax


def apply_axes_bounds(ax, bounds):
    xmin, xmax, ymin, ymax = bounds
    if xmin is not None and xmax is not None:
        ax.set_xlim(xmin, xmax)
    if ymin is not None and ymax is not None:
        ax.set_ylim(ymin, ymax)


# --------------------------------------------------------------------
# HDF5 loader for WriteGridSampleBinary output
# --------------------------------------------------------------------
def load_grid_for_xy_plot(grid_path: Path):
    if not grid_path.is_file():
        raise FileNotFoundError(f"Grid file does not exist: {grid_path}")

    with h5py.File(grid_path, "r") as h5:
        if "/E/field" not in h5:
            raise KeyError("Missing dataset /E/field")
        if "/grid" not in h5:
            raise KeyError("Missing group /grid")

        E = h5["/E/field"][...]  # (Nz, Ny, Nx, dim)
        if E.ndim != 4:
            raise ValueError(f"/E/field must be rank-4, got shape {E.shape}")

        Nz, Ny, Nx, dim = E.shape
        if dim < 2:
            raise ValueError(f"Expected at least 2 field components (Ex,Ey), got dim={dim}")

        g = h5["/grid"]
        if "x" not in g or "y" not in g or "z" not in g:
            raise KeyError("Expected /grid/x, /grid/y, /grid/z (writer now always saves x,y,z).")

        x = np.asarray(g["x"][...], dtype=float)
        y = np.asarray(g["y"][...], dtype=float)
        z = np.asarray(g["z"][...], dtype=float)

        if x.shape != (Nx,):
            raise ValueError(f"/grid/x shape {x.shape} does not match Nx={Nx}")
        if y.shape != (Ny,):
            raise ValueError(f"/grid/y shape {y.shape} does not match Ny={Ny}")
        if z.shape != (Nz,):
            raise ValueError(f"/grid/z shape {z.shape} does not match Nz={Nz}")

        k = Nz // 2
        Ex2 = E[k, :, :, 0]
        Ey2 = E[k, :, :, 1]

        X2 = np.tile(x[None, :], (Ny, 1))
        Y2 = np.tile(y[:, None], (1, Nx))

        return X2, Y2, Ex2, Ey2


def overlay_efield(ax, grid_path: Path, cfg_path: Path, bounds):
    rmin_cfg, rmax_cfg, zmin_cfg, zmax_cfg = load_bounds_from_config(cfg_path)
    xmin, xmax, ymin, ymax = bounds

    rmin = xmin if (xmin is not None and xmax is not None) else rmin_cfg
    rmax = xmax if (xmin is not None and xmax is not None) else rmax_cfg
    zmin = ymin if (ymin is not None and ymax is not None) else zmin_cfg
    zmax = ymax if (ymin is not None and ymax is not None) else zmax_cfg

    R2, Z2, Er2, Ez2 = load_grid_for_xy_plot(grid_path)

    inb = (R2 >= rmin) & (R2 <= rmax) & (Z2 >= zmin) & (Z2 <= zmax)
    if not np.any(inb):
        rmin_grid = float(np.nanmin(R2))
        rmax_grid = float(np.nanmax(R2))
        zmin_grid = float(np.nanmin(Z2))
        zmax_grid = float(np.nanmax(Z2))
        raise ValueError(
            "No grid points fall inside requested bounds.\n"
            f"  Requested r: [{rmin}, {rmax}]\n"
            f"  Requested z: [{zmin}, {zmax}]\n"
            f"  Grid r extent: [{rmin_grid}, {rmax_grid}]\n"
            f"  Grid z extent: [{zmin_grid}, {zmax_grid}]"
        )

    ii, jj = np.where(inb)
    i0, i1 = int(np.min(ii)), int(np.max(ii))
    j0, j1 = int(np.min(jj)), int(np.max(jj))

    R = R2[i0 : i1 + 1, j0 : j1 + 1]
    Z = Z2[i0 : i1 + 1, j0 : j1 + 1]
    Er = Er2[i0 : i1 + 1, j0 : j1 + 1]
    Ez = Ez2[i0 : i1 + 1, j0 : j1 + 1]

    Emag = np.sqrt(Er * Er + Ez * Ez)

    im = ax.pcolormesh(R, Z, Emag, shading="auto", alpha=0.35)
    plt.colorbar(im, ax=ax, pad=0.02, label="|E|", location="left")

    ny, nx = Emag.shape
    step = max(1, int(max(ny, nx) / 80))
    ax.quiver(
        R[::step, ::step],
        Z[::step, ::step],
        Er[::step, ::step],
        Ez[::step, ::step],
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0015,
        alpha=0.55,
    )

    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)


def overlay_efield_gradient(ax, grid_path: Path, cfg_path: Path, which: str, bounds):
    rmin_cfg, rmax_cfg, zmin_cfg, zmax_cfg = load_bounds_from_config(cfg_path)
    xmin, xmax, ymin, ymax = bounds

    rmin = xmin if (xmin is not None and xmax is not None) else rmin_cfg
    rmax = xmax if (xmin is not None and xmax is not None) else rmax_cfg
    zmin = ymin if (ymin is not None and ymax is not None) else zmin_cfg
    zmax = ymax if (ymin is not None and ymax is not None) else zmax_cfg

    R2, Z2, Er2, Ez2 = load_grid_for_xy_plot(grid_path)

    inb = (R2 >= rmin) & (R2 <= rmax) & (Z2 >= zmin) & (Z2 <= zmax)
    if not np.any(inb):
        raise ValueError("No grid points fall inside requested bounds (or config bounds if not overridden).")

    ii, jj = np.where(inb)
    i0, i1 = int(np.min(ii)), int(np.max(ii))
    j0, j1 = int(np.min(jj)), int(np.max(jj))

    R = R2[i0 : i1 + 1, j0 : j1 + 1]
    Z = Z2[i0 : i1 + 1, j0 : j1 + 1]
    Er = Er2[i0 : i1 + 1, j0 : j1 + 1]
    Ez = Ez2[i0 : i1 + 1, j0 : j1 + 1]

    Emag = np.sqrt(Er * Er + Ez * Ez)

    x = R[0, :]
    y = Z[:, 0]
    dEmag_dy, dEmag_dx = np.gradient(Emag, y, x, edge_order=2)

    if which == "x":
        field = dEmag_dx
        cbar_label = "d|E|/dx"
    elif which == "y":
        field = dEmag_dy
        cbar_label = "d|E|/dy"
    else:
        raise ValueError("which must be 'x' or 'y'")

    im = ax.pcolormesh(R, Z, field, shading="auto", alpha=0.55)
    plt.colorbar(im, ax=ax, pad=0.02, label=cbar_label, location="left")

    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)


# --------------------------------------------------------------------
# CSV loader
# --------------------------------------------------------------------
def load_paths(csv_path: Path):
    paths = defaultdict(list)
    with csv_path.open() as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            seed_id = int(row[0])
            r = float(row[2])
            z = float(row[3])
            exit_code = int(row[4])
            paths[seed_id].append((r, z, exit_code))
    return paths


# --------------------------------------------------------------------
# MED BC groups for TPC outline
# --------------------------------------------------------------------
def extract_bc_groups_from_med(mesh: meshio.Mesh, *, prefer_cell_groups=True):
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
        raise ValueError(
            f"No tag array found in mesh.{('cell_data' if prefer_cell_groups else 'point_data')}['{data_key}']."
        )

    group_to_famids = {}
    for fam_id, names in tags_map.items():
        for name in names:
            group_to_famids.setdefault(name, set()).add(int(fam_id))

    bc = {}
    if prefer_cell_groups:
        if not isinstance(data, list) or len(data) != len(mesh.cells):
            raise ValueError("mesh.cell_data['cell_tags'] must be a list aligned with mesh.cells.")

        for group_name, fam_ids in group_to_famids.items():
            fam_ids = set(fam_ids)
            group_cells = []
            used_fam_ids = set()

            for cell_block, tag_arr in zip(mesh.cells, data):
                tag_arr = np.asarray(tag_arr).astype(int)
                mask = np.isin(tag_arr, list(fam_ids))
                if not np.any(mask):
                    continue
                group_cells.append((cell_block.type, np.asarray(cell_block.data)[mask]))
                used_fam_ids.update(np.unique(tag_arr[mask]).tolist())

            if group_cells:
                bc[group_name] = {
                    "points": np.asarray(mesh.points),
                    "cells": group_cells,
                    "family_ids": sorted(used_fam_ids),
                }
    else:
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


def draw_tpc_outline(ax, med_path: Path):
    if not med_path.is_file():
        raise FileNotFoundError(f"MED mesh does not exist: {med_path}")
    mesh = meshio.read(str(med_path))
    bc = extract_bc_groups_from_med(mesh, prefer_cell_groups=True)

    pts = mesh.points
    for g in bc.values():
        for cell_type, conn in g["cells"]:
            conn = np.asarray(conn, dtype=int)

            if cell_type == "line":
                for n0, n1 in conn:
                    ax.plot(
                        [pts[n0, 0], pts[n1, 0]],
                        [pts[n0, 1], pts[n1, 1]],
                        linewidth=0.6,
                        alpha=0.9,
                        color="black",
                    )
            elif cell_type == "line3":
                for n0, n1, n2 in conn:
                    ax.plot(
                        [pts[n0, 0], pts[n1, 0]],
                        [pts[n0, 1], pts[n1, 1]],
                        linewidth=0.6,
                        alpha=0.9,
                        color="black",
                    )
                    ax.plot(
                        [pts[n1, 0], pts[n2, 0]],
                        [pts[n1, 1], pts[n2, 1]],
                        linewidth=0.6,
                        alpha=0.9,
                        color="black",
                    )


# --------------------------------------------------------------------
# Plotting helpers
# --------------------------------------------------------------------
code_colors = {
    0: ("None", "black"),
    1: ("HitCathode", "blue"),
    2: ("HitLiquidGas", "green"),
    3: ("HitWall", "orange"),
    4: ("LeftVolume", "magenta"),
    5: ("MaxSteps", "red"),
    6: ("DegenerateDt", "purple"),
    7: ("HitAxis", "yellow"),
}


def parse_exit_filter(s: str, code_map: dict):
    """
    Returns a set of allowed exit codes, or None if no filter.
    Accepts comma-separated tokens, each token can be integer or a known label key from code_colors.
    Matching is case-insensitive for labels.
    """
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None

    name_to_code = {label.lower(): code for code, (label, _) in code_map.items()}

    allowed = set()
    for tok in s.split(","):
        t = tok.strip()
        if not t:
            continue
        if t.lstrip("+-").isdigit():
            allowed.add(int(t))
            continue
        key = t.lower()
        if key not in name_to_code:
            known = ", ".join(sorted(name_to_code.keys()))
            raise SystemExit(f"Unknown exit label '{t}'. Known labels: {known}")
        allowed.add(name_to_code[key])

    return allowed if allowed else None


def default_title_for_csv(csv_path: Path, sim_results_dir: Path):
    return csv_path.parent.name if csv_path.parent != sim_results_dir else "Debug Electron Paths"


def plot_one(csv_path: Path, grid_path: Path | None, out_png: Path, *, title: str | None = None):
    if not csv_path.is_file():
        print(f"[skip] CSV missing: {csv_path}")
        return

    if (args.efield or args.dEdx or args.dEdy) and (grid_path is None or not grid_path.is_file()):
        print(f"[skip] Grid missing but overlay requested: {grid_path}")
        return

    allowed_exit_codes = parse_exit_filter(args.exit, code_colors)
    paths = load_paths(csv_path)
    bounds = resolve_plot_bounds(args, args.config)

    plt.figure()
    ax = plt.gca()

    if args.dEdx:
        overlay_efield_gradient(ax, grid_path, args.config, which="x", bounds=bounds)
    elif args.dEdy:
        overlay_efield_gradient(ax, grid_path, args.config, which="y", bounds=bounds)
    elif args.efield:
        overlay_efield(ax, grid_path, args.config, bounds=bounds)

    if not args.no_paths:
        for seed_id, pts in paths.items():
            if not pts:
                continue

            exit_code = int(pts[-1][2])
            if allowed_exit_codes is not None and exit_code not in allowed_exit_codes:
                continue

            r = np.array([p[0] for p in pts], dtype=float)
            z = np.array([p[1] for p in pts], dtype=float)

            label, color = code_colors.get(exit_code, ("Unknown", "gray"))
            ax.plot(r, z, "-", linewidth=0.8, color=color, alpha=0.5)

            if len(paths) == 1 and r.size >= 3:
                ax.scatter(r[:-2], z[:-2], s=1, marker="o", color=color)

            ax.scatter([r[-1]], [z[-1]], s=12, marker="x", color=color)

            if r.size >= 4:
                dr = (r - np.roll(r, 1))[1:-2]
                dz = (z - np.roll(z, 1))[1:-2]
                if np.isclose(dr, 0.0).all() and np.isclose(dz, 0.0).all():
                    ax.scatter([r[0]], [z[0]], s=18, marker="x", color=color)

    if plot_TPC:
        ax.set_aspect("equal", "box")
        draw_tpc_outline(ax, args.mesh_med)

    ax.set_xlabel("r")
    ax.set_ylabel("z")
    ax.grid(True)

    if title is None:
        title = default_title_for_csv(csv_path, args.sim_results)
    ax.set_title(title)

    apply_axes_bounds(ax, bounds)

    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker="x",
            linestyle="None",
            color=color,
            markeredgecolor=color,
            markersize=6,
            markeredgewidth=1.5,
            label=label,
        )
        for _, (label, color) in sorted(code_colors.items())
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        frameon=True,
        title="Exit code",
    )

    plt.tight_layout(rect=[0, 0, 0.82, 1])

    dpi = args.dpi if args.dpi is not None else (1800 if plot_TPC else 180)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi)
    plt.close()

    print(f"Saved plot to: {out_png}")


def discover_runs(sim_results_dir: Path):
    if not sim_results_dir.is_dir():
        raise FileNotFoundError(f"sim_results directory does not exist: {sim_results_dir}")
    return sorted([p for p in sim_results_dir.iterdir() if p.is_dir() and p.name.startswith("run_")])


# --------------------------------------------------------------------
# Main dispatch
# --------------------------------------------------------------------
if args.csv is not None:
    # Single-run mode (explicit)
    csv_path = args.csv
    grid_path = args.grid
    out_png = csv_path.with_suffix(".png")

    run_name = csv_path.parent.name
    preset = meta_run_to_preset.get(run_name)
    if preset is None and meta_run_to_preset:
        print(f"[meta] No preset entry for {run_name} in meta.txt (using default title)")

    title = build_plot_title(
        default_title=default_title_for_csv(csv_path, args.sim_results),
        sim_id=meta_sim_id,
        preset=preset,
    )
    plot_one(csv_path, grid_path, out_png, title=title)
    sys.exit(0)

# Batch mode: default to sim_results; iterate run_xxxx if present, else use sim_results root paths
sim_results = args.sim_results.resolve()
run_dirs = discover_runs(sim_results)

if run_dirs:
    for run_dir in run_dirs:
        csv_path = run_dir / REL_CSV
        grid_path = run_dir / REL_GRID
        out_png = csv_path.with_suffix(".png")  # saves inside the run_dir

        preset = meta_run_to_preset.get(run_dir.name)
        if preset is None and meta_run_to_preset:
            print(f"[meta] No preset entry for {run_dir.name} in meta.txt (using default title)")

        title = build_plot_title(
            default_title=run_dir.name,
            sim_id=meta_sim_id,
            preset=preset,
        )
        plot_one(csv_path, grid_path, out_png, title=title)
else:
    csv_path = sim_results / REL_CSV
    grid_path = sim_results / REL_GRID
    out_png = csv_path.with_suffix(".png")  # saves inside sim_results

    title = build_plot_title(
        default_title="Debug Electron Paths",
        sim_id=meta_sim_id,
        preset=None,
    )
    plot_one(csv_path, grid_path, out_png, title=title)
