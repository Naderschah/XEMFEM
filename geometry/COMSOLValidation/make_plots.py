"""
Compare COMSOL vs XEMFEM on the XEMFEM grid and produce plots per case:
  - comparison_voltage.png
  - comparison_E_magnitude.png
  - comparison_E_direction.png
Inputs per case directory:
  - one COMSOL*.txt export containing "NASTRAN" in filename
    (e.g. COMSOL_Grid_NASTRAN.txt or COMSOLGrid_NASTRAN.txt),
    with columns: x y V Ex Ey
      Header example:
        % x   y   V (V)   es.Ex (V/m)   es.Ey (V/m)
  - XEMFEM H5 at: interpolated/interpolated.h5 containing:
      /V/field   Dataset {Nz, Ny, Nx}  (your example: {1,1200,2000})
      /E/field   Dataset {Nz, Ny, Nx, 2}  (your example: {1,1200,2000,2})
      /grid/x    Dataset {Nx}
      /grid/y    Dataset {Ny}
      /grid/z    Dataset {Nz}

We STOP computing E from gradients. We use E directly from:
  - COMSOL text columns Ex/Ey
  - XEMFEM /E/field (component 0 -> Ex, component 1 -> Ey)

Note if running in the XEMFEM env do:
python -m pip install --no-cache-dir --force-reinstall \
  "numpy<2" \
  "matplotlib<4" \
  "h5py" \
  "meshio" \
  "pyvista"
"""

from __future__ import annotations

import argparse
import glob
import os
import xml.etree.ElementTree as ET

import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

try:
    import meshio
except ImportError:
    meshio = None

try:
    import pyvista as pv
except ImportError:
    pv = None


IMSHOW_INTERPOLATIONS: tuple[str, ...] = (
    "none",
    "nearest",
    "bilinear",
    "bicubic",
    "spline16",
    "spline36",
    "hanning",
    "hamming",
    "hermite",
    "kaiser",
    "quadric",
    "catrom",
    "gaussian",
    "bessel",
    "mitchell",
    "sinc",
    "lanczos",
    "blackman",
)


# -----------------------------
# Text loading (COMSOL)
# -----------------------------

def _find_first_data_row(txt_path: str, max_scan_lines: int = 20000) -> int:
    """
    Returns the 0-based line index where numeric data starts (first token parses as float).
    Skips COMSOL '%' comment lines and the column header line.
    """
    with open(txt_path, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if i >= max_scan_lines:
                raise RuntimeError(
                    f"Could not find data start within first {max_scan_lines} lines."
                )
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            tok = s.split(None, 1)[0]
            try:
                float(tok)
                return i
            except ValueError:
                continue
    raise RuntimeError("Empty file or no numeric data found.")


def load_comsol_to_target_grid_V_Ex_Ey(
    comsol_txt_path: str,
    x_target: np.ndarray,
    y_target: np.ndarray,
    *,
    chunksize: int = 2_000_000,
    coord_tol: float = 1e-6,
    dtype_out=np.float64,
):
    """
    Strict loader: requires COMSOL x,y to be EXACTLY on the target grid values.
    Fails if:
      - any (x,y) is not exactly in (x_target,y_target)
      - any grid cell is missing
      - any grid cell is written more than once

    Returns:
      Vt, Ext, Eyt: arrays (Ny, Nx)
    """
    x_target = np.asarray(x_target)
    y_target = np.asarray(y_target)
    if x_target.ndim != 1 or y_target.ndim != 1:
        raise ValueError("x_target and y_target must be 1D arrays.")

    Nx = x_target.size
    Ny = y_target.size
    if Nx == 0 or Ny == 0:
        raise ValueError("Empty target grid.")

    Vt = np.empty((Ny, Nx), dtype=dtype_out)
    Ext = np.empty((Ny, Nx), dtype=dtype_out)
    Eyt = np.empty((Ny, Nx), dtype=dtype_out)

    Vt.fill(np.nan)
    Ext.fill(np.nan)
    Eyt.fill(np.nan)

    seen = np.zeros((Ny, Nx), dtype=np.bool_)
    seen_x = np.zeros((Nx,), dtype=np.bool_)
    seen_y = np.zeros((Ny,), dtype=np.bool_)

    # Diagnostics buffers:
    # - coordinate deltas to nearest target grid point
    # - tolerance-binned unique coordinates to detect grid-size mismatch
    dx_samples: list[float] = []
    dy_samples: list[float] = []
    x_bins: set[int] = set()
    y_bins: set[int] = set()

    data_start = _find_first_data_row(comsol_txt_path)

    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for efficient chunked reading of very large COMSOL text files.\n"
            "Install with: pip install pandas"
        ) from e

    reader = pd.read_csv(
        comsol_txt_path,
        sep=r"\s+",
        header=None,
        names=["x", "y", "V", "Ex", "Ey"],
        comment="%",
        skiprows=data_start,
        engine="c",
        dtype=np.float64,
        chunksize=chunksize,
    )

    tol = float(coord_tol)
    for chunk in reader:
        xs = chunk["x"].to_numpy()
        ys = chunk["y"].to_numpy()
        vs = chunk["V"].to_numpy()
        exs = chunk["Ex"].to_numpy()
        eys = chunk["Ey"].to_numpy()

        # Approximate unique coordinate counts with tolerance binning.
        # This avoids false inflation from tiny text jitter below `tol`.
        x_bins.update(np.rint(xs / tol).astype(np.int64).tolist())
        y_bins.update(np.rint(ys / tol).astype(np.int64).tolist())

        # Map each point via exact dictionary lookup
        # (Loop is intentional to guarantee exact matching and provide good error messages.)
        for x, y, v, ex, ey in zip(xs, ys, vs, exs, eys):
            xv = float(x)
            yv = float(y)

            # ---- x coordinate ----
            ix_closest = int(np.argmin(np.abs(x_target - xv)))
            x_closest = x_target[ix_closest]
            dx = abs(xv - x_closest)

            if dx > tol:
                raise ValueError(
                    "COMSOL x value not on target grid (within tolerance): "
                    f"x={xv!r}. Closest grid x[{ix_closest}]={x_closest!r}, |Δx|={dx:.3e} > {tol:.1e}"
                )

            # snap to grid
            ix = ix_closest

            # ---- y coordinate ----
            iy_closest = int(np.argmin(np.abs(y_target - yv)))
            y_closest = y_target[iy_closest]
            dy = abs(yv - y_closest)

            if dy > tol:
                raise ValueError(
                    "COMSOL y value not on target grid (within tolerance): "
                    f"y={yv!r}. Closest grid y[{iy_closest}]={y_closest!r}, |Δy|={dy:.3e} > {tol:.1e}"
                )

            # snap to grid
            iy = iy_closest
            if seen[iy, ix]:
                raise ValueError(f"Duplicate COMSOL sample at grid cell (iy={iy}, ix={ix}) for (x={xv}, y={yv})")

            dx_samples.append(xv - x_closest)
            dy_samples.append(yv - y_closest)

            Vt[iy, ix] = v
            Ext[iy, ix] = ex
            Eyt[iy, ix] = ey
            seen[iy, ix] = True
            seen_x[ix] = True
            seen_y[iy] = True

    # Ensure fully populated
    if not np.all(seen):
        missing = np.argwhere(~seen)
        # show first few missing cells
        preview = missing[:10]
        raise ValueError(
            f"COMSOL grid does not fully cover target grid: missing {missing.shape[0]} cells. "
            f"First missing indices (iy,ix): {preview.tolist()}"
        )

    # Explicit grid size checks (with tolerance-binned uniqueness)
    nx_comsol = len(x_bins)
    ny_comsol = len(y_bins)
    if nx_comsol != Nx or ny_comsol != Ny:
        raise ValueError(
            "Grid size mismatch between COMSOL and target grid: "
            f"COMSOL approx unique (Nx,Ny)=({nx_comsol},{ny_comsol}), "
            f"target (Nx,Ny)=({Nx},{Ny}), tol={tol:.1e}. "
            "Likely a different grid or coordinate jitter above tolerance."
        )

    if not np.all(seen_x) or not np.all(seen_y):
        raise ValueError(
            "Grid axis coverage mismatch: not all target x/y indices were hit by COMSOL samples."
        )

    # Coordinate-difference diagnostics: classify as aligned/systematic/noisy.
    def _classify_delta(name: str, arr: np.ndarray, atol: float) -> str:
        arr = np.asarray(arr, dtype=np.float64)
        if arr.size == 0:
            return f"{name}: no samples"
        mean = float(np.mean(arr))
        std = float(np.std(arr))
        maxabs = float(np.max(np.abs(arr)))
        if maxabs <= atol:
            kind = "aligned (within tolerance)"
        elif std <= 5.0 * atol:
            kind = "systematic offset (nearly constant shift)"
        elif abs(mean) > 3.0 * std:
            kind = "systematic offset + jitter"
        else:
            kind = "noise/jitter dominated"
        return (
            f"{name}: {kind}; "
            f"mean={mean:.3e}, std={std:.3e}, max|d|={maxabs:.3e}, tol={atol:.1e}"
        )

    print("[GRID CHECK] size match: "
          f"COMSOL (Nx,Ny)=({nx_comsol},{ny_comsol}) vs target ({Nx},{Ny})")
    print("[GRID CHECK] " + _classify_delta("x delta", np.asarray(dx_samples), tol))
    print("[GRID CHECK] " + _classify_delta("y delta", np.asarray(dy_samples), tol))

    return Vt, Ext, Eyt


def load_comsol_native_regular_grid_V_Ex_Ey(
    comsol_txt_path: str,
    *,
    chunksize: int = 2_000_000,
    coord_tol: float = 1e-6,
    dtype_out=np.float64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Load COMSOL txt on its native regular grid.
    Returns:
      xg, yg, Vg, Exg, Eyg, Dg where:
        - Vg/Exg/Eyg shape is (Ny, Nx)
        - Dg is material.domain on the same grid if present, else None.
    """
    data_start = _find_first_data_row(comsol_txt_path)

    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for efficient chunked reading of very large COMSOL text files.\n"
            "Install with: pip install pandas"
        ) from e

    reader = pd.read_csv(
        comsol_txt_path,
        sep=r"\s+",
        header=None,
        names=["x", "y", "V", "Ex", "Ey", "matdom"],
        comment="%",
        skiprows=data_start,
        engine="c",
        dtype=np.float64,
        chunksize=chunksize,
    )

    x_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    v_parts: list[np.ndarray] = []
    ex_parts: list[np.ndarray] = []
    ey_parts: list[np.ndarray] = []
    dom_parts: list[np.ndarray] = []

    for chunk in reader:
        x_parts.append(chunk["x"].to_numpy())
        y_parts.append(chunk["y"].to_numpy())
        v_parts.append(chunk["V"].to_numpy())
        ex_parts.append(chunk["Ex"].to_numpy())
        ey_parts.append(chunk["Ey"].to_numpy())
        dom_parts.append(chunk["matdom"].to_numpy())

    if not x_parts:
        raise ValueError("COMSOL text file has no data rows.")

    xs = np.concatenate(x_parts)
    ys = np.concatenate(y_parts)
    vs = np.concatenate(v_parts)
    exs = np.concatenate(ex_parts)
    eys = np.concatenate(ey_parts)
    doms = np.concatenate(dom_parts)

    tol = float(coord_tol)
    x_bins = np.rint(xs / tol).astype(np.int64)
    y_bins = np.rint(ys / tol).astype(np.int64)

    _, x_inv = np.unique(x_bins, return_inverse=True)
    _, y_inv = np.unique(y_bins, return_inverse=True)

    x_counts = np.bincount(x_inv)
    y_counts = np.bincount(y_inv)
    if np.any(x_counts == 0) or np.any(y_counts == 0):
        raise ValueError("Internal error: empty bin while building COMSOL native grid.")

    x_vals = np.bincount(x_inv, weights=xs) / x_counts
    y_vals = np.bincount(y_inv, weights=ys) / y_counts

    x_order = np.argsort(x_vals)
    y_order = np.argsort(y_vals)

    x_rank = np.empty_like(x_order)
    y_rank = np.empty_like(y_order)
    x_rank[x_order] = np.arange(x_order.size)
    y_rank[y_order] = np.arange(y_order.size)

    ix = x_rank[x_inv]
    iy = y_rank[y_inv]

    xg = x_vals[x_order].astype(dtype_out, copy=False)
    yg = y_vals[y_order].astype(dtype_out, copy=False)
    Nx = int(xg.size)
    Ny = int(yg.size)

    Vg = np.full((Ny, Nx), np.nan, dtype=dtype_out)
    Exg = np.full((Ny, Nx), np.nan, dtype=dtype_out)
    Eyg = np.full((Ny, Nx), np.nan, dtype=dtype_out)
    Dg: np.ndarray | None = None
    has_domain = np.isfinite(doms).any()
    if has_domain:
        Dg = np.full((Ny, Nx), np.nan, dtype=dtype_out)
    seen = np.zeros((Ny, Nx), dtype=np.bool_)

    for p in range(xs.size):
        i = int(ix[p])
        j = int(iy[p])
        if seen[j, i]:
            raise ValueError(
                f"Duplicate COMSOL sample at native grid cell (iy={j}, ix={i}), "
                f"raw point=({xs[p]}, {ys[p]})"
            )
        Vg[j, i] = vs[p]
        Exg[j, i] = exs[p]
        Eyg[j, i] = eys[p]
        if has_domain and Dg is not None:
            Dg[j, i] = doms[p]
        seen[j, i] = True

    if not np.all(seen):
        missing = np.argwhere(~seen)
        preview = missing[:10]
        raise ValueError(
            f"COMSOL native grid not fully populated: missing {missing.shape[0]} cells. "
            f"First missing indices (iy,ix): {preview.tolist()}"
        )

    return xg, yg, Vg, Exg, Eyg, Dg


def _check_axis_strictly_increasing(name: str, axis: np.ndarray) -> None:
    d = np.diff(axis)
    if np.any(d <= 0):
        raise ValueError(f"Axis '{name}' must be strictly increasing for bilinear resampling.")


def resample_regular_grid_bilinear(
    field_src: np.ndarray,
    x_src: np.ndarray,
    y_src: np.ndarray,
    x_dst: np.ndarray,
    y_dst: np.ndarray,
    *,
    fill_value: float = np.nan,
    coord_tol: float = 1e-12,
) -> np.ndarray:
    """
    Bilinear resampling from source regular grid (x_src,y_src) to destination regular grid (x_dst,y_dst).
    field_src must have shape (Ny_src, Nx_src).
    """
    field_src = np.asarray(field_src, dtype=np.float64)
    x_src = np.asarray(x_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    x_dst = np.asarray(x_dst, dtype=np.float64)
    y_dst = np.asarray(y_dst, dtype=np.float64)

    if field_src.shape != (y_src.size, x_src.size):
        raise ValueError(
            f"field_src shape {field_src.shape} does not match (Ny,Nx)=({y_src.size},{x_src.size})"
        )
    if x_src.ndim != 1 or y_src.ndim != 1 or x_dst.ndim != 1 or y_dst.ndim != 1:
        raise ValueError("x/y source and destination axes must be 1D.")

    _check_axis_strictly_increasing("x_src", x_src)
    _check_axis_strictly_increasing("y_src", y_src)
    _check_axis_strictly_increasing("x_dst", x_dst)
    _check_axis_strictly_increasing("y_dst", y_dst)

    epsx = float(coord_tol) + 1e-12 * max(1.0, abs(x_src[-1] - x_src[0]))
    epsy = float(coord_tol) + 1e-12 * max(1.0, abs(y_src[-1] - y_src[0]))

    inside_x = (x_dst >= x_src[0] - epsx) & (x_dst <= x_src[-1] + epsx)
    inside_y = (y_dst >= y_src[0] - epsy) & (y_dst <= y_src[-1] + epsy)

    xq = np.clip(x_dst, x_src[0], x_src[-1])
    yq = np.clip(y_dst, y_src[0], y_src[-1])

    ix1 = np.searchsorted(x_src, xq, side="right")
    iy1 = np.searchsorted(y_src, yq, side="right")
    ix0 = np.clip(ix1 - 1, 0, x_src.size - 1)
    iy0 = np.clip(iy1 - 1, 0, y_src.size - 1)
    ix1 = np.clip(ix1, 0, x_src.size - 1)
    iy1 = np.clip(iy1, 0, y_src.size - 1)

    dx = x_src[ix1] - x_src[ix0]
    dy = y_src[iy1] - y_src[iy0]
    wx = np.zeros_like(xq)
    wy = np.zeros_like(yq)
    mx = dx > 0
    my = dy > 0
    wx[mx] = (xq[mx] - x_src[ix0[mx]]) / dx[mx]
    wy[my] = (yq[my] - y_src[iy0[my]]) / dy[my]

    f00 = field_src[np.ix_(iy0, ix0)]
    f10 = field_src[np.ix_(iy0, ix1)]
    f01 = field_src[np.ix_(iy1, ix0)]
    f11 = field_src[np.ix_(iy1, ix1)]

    wx2 = wx[None, :]
    wy2 = wy[:, None]
    out = (
        (1.0 - wx2) * (1.0 - wy2) * f00
        + wx2 * (1.0 - wy2) * f10
        + (1.0 - wx2) * wy2 * f01
        + wx2 * wy2 * f11
    )

    valid = inside_y[:, None] & inside_x[None, :]
    out[~valid] = fill_value
    return out


def align_field_to_target_grid(
    field_src: np.ndarray,
    x_src: np.ndarray,
    y_src: np.ndarray,
    x_dst: np.ndarray,
    y_dst: np.ndarray,
    *,
    coord_tol: float,
    label_src: str,
    label_dst: str,
) -> np.ndarray:
    """
    Return field on destination grid.
    - If source and destination axes match within tolerance, return source field.
    - Otherwise bilinearly resample source -> destination.
    """
    field_src = np.asarray(field_src, dtype=np.float64)
    x_src = np.asarray(x_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    x_dst = np.asarray(x_dst, dtype=np.float64)
    y_dst = np.asarray(y_dst, dtype=np.float64)

    if field_src.shape != (y_src.size, x_src.size):
        raise ValueError(
            f"{label_src} field shape {field_src.shape} does not match source grid "
            f"(Ny,Nx)=({y_src.size},{x_src.size})"
        )

    same_sizes = (x_src.size == x_dst.size) and (y_src.size == y_dst.size)
    if same_sizes:
        max_dx = float(np.max(np.abs(x_src - x_dst))) if x_src.size else 0.0
        max_dy = float(np.max(np.abs(y_src - y_dst))) if y_src.size else 0.0
        if max_dx <= coord_tol and max_dy <= coord_tol:
            return field_src

    out = resample_regular_grid_bilinear(
        field_src, x_src, y_src, x_dst, y_dst, coord_tol=coord_tol
    )
    print(
        f"  Aligned '{label_src}' onto '{label_dst}' grid via bilinear resampling "
        f"({x_src.size}x{y_src.size} -> {x_dst.size}x{y_dst.size})."
    )
    return out


def _nearest_axis_deltas(query: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """
    For each point in query, return delta to nearest point in reference.
    delta = query - nearest(reference).
    """
    query = np.asarray(query, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if reference.size == 0:
        return np.full_like(query, np.nan, dtype=np.float64)

    idx_r = np.searchsorted(reference, query, side="left")
    idx_l = np.clip(idx_r - 1, 0, reference.size - 1)
    idx_r = np.clip(idx_r, 0, reference.size - 1)

    dl = np.abs(query - reference[idx_l])
    dr = np.abs(query - reference[idx_r])
    use_r = dr < dl
    nearest = reference[idx_l]
    nearest[use_r] = reference[idx_r[use_r]]
    return query - nearest


def _classify_delta_series(name: str, deltas: np.ndarray, atol: float) -> str:
    deltas = np.asarray(deltas, dtype=np.float64)
    if deltas.size == 0:
        return f"{name}: no samples"
    mean = float(np.mean(deltas))
    std = float(np.std(deltas))
    maxabs = float(np.max(np.abs(deltas)))
    if maxabs <= atol:
        kind = "aligned (within tolerance)"
    elif std <= 5.0 * atol:
        kind = "systematic offset (nearly constant shift)"
    elif abs(mean) > 3.0 * std:
        kind = "systematic offset + jitter"
    else:
        kind = "noise/jitter dominated"
    return (
        f"{name}: {kind}; "
        f"mean={mean:.3e}, std={std:.3e}, max|d|={maxabs:.3e}, tol={atol:.1e}"
    )


def print_grid_alignment_report(
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    x_cmp: np.ndarray,
    y_cmp: np.ndarray,
    *,
    coord_tol: float,
    label_ref: str = "reference",
    label_cmp: str = "comparison",
) -> None:
    tol = float(coord_tol)
    print(
        f"[GRID CHECK] sizes: {label_ref}(Nx,Ny)=({x_ref.size},{y_ref.size}) "
        f"vs {label_cmp}(Nx,Ny)=({x_cmp.size},{y_cmp.size})"
    )
    dx = _nearest_axis_deltas(x_cmp, x_ref)
    dy = _nearest_axis_deltas(y_cmp, y_ref)
    print("[GRID CHECK] " + _classify_delta_series(f"x delta {label_cmp}->{label_ref}", dx, tol))
    print("[GRID CHECK] " + _classify_delta_series(f"y delta {label_cmp}->{label_ref}", dy, tol))


def assert_same_meshgrid_coordinates(
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    x_cmp: np.ndarray,
    y_cmp: np.ndarray,
    *,
    coord_tol: float,
    label_ref: str = "reference",
    label_cmp: str = "comparison",
) -> None:
    tol = float(coord_tol)
    if x_ref.size != x_cmp.size or y_ref.size != y_cmp.size:
        raise ValueError(
            f"Grid size mismatch: {label_ref}(Nx,Ny)=({x_ref.size},{y_ref.size}) "
            f"vs {label_cmp}(Nx,Ny)=({x_cmp.size},{y_cmp.size})."
        )

    dx = x_cmp - x_ref
    dy = y_cmp - y_ref
    max_dx = float(np.max(np.abs(dx))) if dx.size else 0.0
    max_dy = float(np.max(np.abs(dy))) if dy.size else 0.0

    print("[GRID CHECK] " + _classify_delta_series(f"x direct {label_cmp}-{label_ref}", dx, tol))
    print("[GRID CHECK] " + _classify_delta_series(f"y direct {label_cmp}-{label_ref}", dy, tol))

    if max_dx > tol or max_dy > tol:
        raise ValueError(
            "Meshgrid coordinates are not aligned within tolerance: "
            f"max|dx|={max_dx:.3e}, max|dy|={max_dy:.3e}, tol={tol:.1e}. "
            "Use --resample-xemfem-to-comsol to compare on COMSOL grid."
        )

# -----------------------------
# H5 loading (XEMFEM)
# -----------------------------

def load_V_midplane(h5_path: str, k: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load /V/field and /grid/{x,y,z} and return a z-midplane slice.

    Supports /V/field shapes:
      - (Nz, Ny, Nx)
      - (Ny, Nx)
      - (dim, Ny, Nx) with dim==1 (rare; handled)
      - (Nz, Ny, Nx, dim) with dim==1 (rare; handled)

    Returns:
      X2, Y2, V2 with V2 shape (Ny, Nx)
    """
    with h5py.File(h5_path, "r") as h5:
        if "/V/field" not in h5:
            raise KeyError("Missing dataset /V/field")
        if "/grid" not in h5:
            raise KeyError("Missing group /grid")

        V = h5["/V/field"][...]
        g = h5["/grid"]

        for axis in ("x", "y", "z"):
            if axis not in g:
                raise KeyError("Expected /grid/x, /grid/y, /grid/z")

        x = np.asarray(g["x"][...], dtype=float)
        y = np.asarray(g["y"][...], dtype=float)
        z = np.asarray(g["z"][...], dtype=float)

        if V.ndim == 3:
            Nz, Ny, Nx = V.shape
            if k is None:
                k = Nz // 2
            if not (0 <= k < Nz):
                raise ValueError(f"k={k} out of range for Nz={Nz}")
            V2 = V[k, :, :]
        elif V.ndim == 2:
            Ny, Nx = V.shape
            V2 = V
        elif V.ndim == 4:
            Nz, Ny, Nx, dim = V.shape
            if dim < 1:
                raise ValueError(f"/V/field has invalid dim={dim}")
            if k is None:
                k = Nz // 2
            if not (0 <= k < Nz):
                raise ValueError(f"k={k} out of range for Nz={Nz}")
            V2 = V[k, :, :, 0]
        else:
            raise ValueError(f"/V/field must be rank-2/3/4, got shape {V.shape}")

        if x.shape != (Nx,):
            raise ValueError(f"/grid/x shape {x.shape} does not match Nx={Nx}")
        if y.shape != (Ny,):
            raise ValueError(f"/grid/y shape {y.shape} does not match Ny={Ny}")

        X2 = np.tile(x[None, :], (Ny, 1))
        Y2 = np.tile(y[:, None], (1, Nx))
        return X2, Y2, np.asarray(V2, dtype=np.float64)

def E_direction_mismatch_acos_deg(
    Ex_a: np.ndarray,
    Ey_a: np.ndarray,
    Ex_b: np.ndarray,
    Ey_b: np.ndarray,
    *,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Angle between vectors A and B in degrees:
      acos( dot(A,B) / (|A||B|) )
    Returns NaN where either magnitude is too small.
    """
    Ex_a = np.asarray(Ex_a, dtype=np.float64)
    Ey_a = np.asarray(Ey_a, dtype=np.float64)
    Ex_b = np.asarray(Ex_b, dtype=np.float64)
    Ey_b = np.asarray(Ey_b, dtype=np.float64)

    dot = Ex_a * Ex_b + Ey_a * Ey_b
    mag_a = np.hypot(Ex_a, Ey_a)
    mag_b = np.hypot(Ex_b, Ey_b)
    denom = mag_a * mag_b

    out = np.full(dot.shape, np.nan, dtype=np.float64)
    valid = denom > eps

    cosang = np.empty_like(dot, dtype=np.float64)
    cosang[valid] = dot[valid] / denom[valid]
    cosang[valid] = np.clip(cosang[valid], -1.0, 1.0)

    out[valid] = np.degrees(np.arccos(cosang[valid]))
    return out

def load_E_midplane(h5_path: str, k: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Load /E/field and /grid/{x,y,z} and return a z-midplane slice.

    Expected /E/field shapes:
      - (Nz, Ny, Nx, 2)  (your case)
    Also supports:
      - (Ny, Nx, 2)
      - (2, Ny, Nx)

    Returns:
      X2, Y2, Ex2, Ey2 (each (Ny, Nx))
    """
    with h5py.File(h5_path, "r") as h5:
        if "/E/field" not in h5:
            raise KeyError("Missing dataset /E/field")
        if "/grid" not in h5:
            raise KeyError("Missing group /grid")

        E = h5["/E/field"][...]
        g = h5["/grid"]

        for axis in ("x", "y", "z"):
            if axis not in g:
                raise KeyError("Expected /grid/x, /grid/y, /grid/z")

        x = np.asarray(g["x"][...], dtype=float)
        y = np.asarray(g["y"][...], dtype=float)
        z = np.asarray(g["z"][...], dtype=float)

        if E.ndim == 4:
            Nz, Ny, Nx, dim = E.shape
            if dim < 2:
                raise ValueError(f"/E/field last dimension must be 2, got {dim}")
            if k is None:
                k = Nz // 2
            if not (0 <= k < Nz):
                raise ValueError(f"k={k} out of range for Nz={Nz}")
            E2 = E[k, :, :, :]  # (Ny, Nx, 2)
        elif E.ndim == 3:
            a, b, c = E.shape
            if a == 2 and b == y.size and c == x.size:
                # (2, Ny, Nx)
                E2 = np.moveaxis(E, 0, -1)  # (Ny, Nx, 2)
                Ny, Nx, dim = E2.shape
            elif c == 2 and a == y.size and b == x.size:
                # (Ny, Nx, 2)
                E2 = E
                Ny, Nx, dim = E2.shape
            else:
                raise ValueError(f"Unsupported /E/field shape {E.shape}; expected (...,2) components")
        else:
            raise ValueError(f"/E/field must be rank-3/4, got shape {E.shape}")

        if x.shape != (Nx,):
            raise ValueError(f"/grid/x shape {x.shape} does not match Nx={Nx}")
        if y.shape != (Ny,):
            raise ValueError(f"/grid/y shape {y.shape} does not match Ny={Ny}")

        X2 = np.tile(x[None, :], (Ny, 1))
        Y2 = np.tile(y[:, None], (1, Nx))

        Ex2 = np.asarray(E2[:, :, 0], dtype=np.float64)
        Ey2 = np.asarray(E2[:, :, 1], dtype=np.float64)
        return X2, Y2, Ex2, Ey2


# -----------------------------
# Metrics / utilities
# -----------------------------

def symmetric_error(A: np.ndarray, B: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Symmetric error: 2(A-B) / max(|A|+|B|, eps)
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    eps = float(eps)
    if not np.isfinite(eps) or eps <= 0.0:
        eps = np.finfo(np.float64).eps
    denom = np.clip(np.abs(A) + np.abs(B), eps, None)
    return np.abs(2.0 * (A - B) / denom)


def relative_error_to_reference(A: np.ndarray, B: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Signed relative error to reference:
      (A - B) / B
    Denominator is clipped in magnitude to avoid divide-by-zero.
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    eps = float(eps)
    if not np.isfinite(eps) or eps <= 0.0:
        eps = np.finfo(np.float64).eps
    denom = np.where(np.abs(B) >= eps, B, np.where(B < 0.0, -eps, eps))
    return (A - B) / denom


def E_magnitude(Ex: np.ndarray, Ey: np.ndarray) -> np.ndarray:
    return np.hypot(np.asarray(Ex, dtype=np.float64), np.asarray(Ey, dtype=np.float64))


def E_direction_deg(Ex: np.ndarray, Ey: np.ndarray) -> np.ndarray:
    return np.degrees(np.arctan2(np.asarray(Ey, dtype=np.float64), np.asarray(Ex, dtype=np.float64)))


def angle_diff_deg(a_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    """
    Smallest signed angular difference (a - b) in degrees, wrapped to [-180, 180].
    """
    a = np.asarray(a_deg, dtype=np.float64)
    b = np.asarray(b_deg, dtype=np.float64)
    return (a - b + 180.0) % 360.0 - 180.0


def _shared_linear_limits_for_maps(maps: list[np.ndarray]) -> tuple[float, float] | None:
    vals: list[np.ndarray] = []
    for m in maps:
        arr = np.asarray(m, dtype=np.float64)
        f = arr[np.isfinite(arr)]
        if f.size > 0:
            vals.append(f.ravel())
    if not vals:
        return None

    allv = np.concatenate(vals)
    vmin = float(np.min(allv))
    vmax = float(np.max(allv))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return None
    if vmax <= vmin:
        eps = max(1e-12, abs(vmin) * 1e-6 + 1e-12)
        vmax = vmin + eps
    return vmin, vmax


def _shared_log_limits_for_maps(
    maps: list[np.ndarray],
    *,
    enabled: bool,
    floor: float = 1e-12,
) -> tuple[float, float] | None:
    if not enabled:
        return None

    vals: list[np.ndarray] = []
    for m in maps:
        arr = np.asarray(m, dtype=np.float64)
        f = arr[np.isfinite(arr)]
        if f.size > 0:
            vals.append(f.ravel())
    if not vals:
        return None

    allv = np.concatenate(vals)
    if np.any(allv < 0.0):
        return None

    pos = allv[allv > 0.0]
    if pos.size == 0:
        return None

    vmin = max(float(np.min(pos)), float(floor))
    vmax = float(np.max(allv))
    if not np.isfinite(vmax) or vmax <= 0.0:
        return None
    if vmax <= vmin:
        vmax = vmin * 10.0

    return vmin, vmax


def _build_bracketed_log_norm(
    limits: tuple[float, float] | None,
    *,
    bins: int,
) -> tuple[mcolors.Normalize | None, np.ndarray | None]:
    if limits is None:
        return None, None

    vmin, vmax = limits
    nbins = int(max(2, bins))
    boundaries = np.geomspace(vmin, vmax, nbins + 1)
    norm = mcolors.BoundaryNorm(boundaries, ncolors=nbins, clip=True)
    return norm, boundaries


def _build_sym_err_decade_log_norm(
    *,
    floor: float = 1e-12,
    cmap_bins: int | None = None,
) -> tuple[mcolors.Normalize, np.ndarray, np.ndarray, list[str]]:
    """
    Fixed logarithmic bins for symmetric error:
      - hard upper limit at 1 (100%)
      - one color per decade down to `floor`
    """
    f = float(floor)
    if not np.isfinite(f) or f <= 0.0:
        f = 1e-12
    f = min(f, 1.0)

    exp_min = int(np.floor(np.log10(f)))
    exp_min = min(exp_min, 0)
    exps = np.arange(exp_min, 1, dtype=np.int64)  # ...,-2,-1,0
    boundaries = (10.0 ** exps.astype(np.float64)).astype(np.float64)
    if boundaries.size < 2:
        boundaries = np.array([0.1, 1.0], dtype=np.float64)

    n_intervals = int(boundaries.size) - 1
    if cmap_bins is None:
        norm_ncolors = n_intervals
    else:
        norm_ncolors = max(n_intervals, int(max(2, cmap_bins)))
    norm = mcolors.BoundaryNorm(boundaries, ncolors=norm_ncolors, clip=True)

    # Show ticks as percentages, top-to-bottom: 100%, 10%, 1%, ...
    ticks = boundaries[::-1]
    tick_labels = [f"{100.0 * t:g}%" for t in ticks]
    return norm, boundaries, ticks, tick_labels


def _prepare_map_for_log_norm(data: np.ndarray, norm: mcolors.Normalize | None) -> np.ndarray:
    if norm is None:
        return data
    out = np.asarray(data, dtype=np.float64).copy()
    m = np.isfinite(out) & (out <= 0.0)
    vmin = None
    if hasattr(norm, "vmin") and getattr(norm, "vmin") is not None:
        vmin = float(getattr(norm, "vmin"))
    elif hasattr(norm, "boundaries"):
        b = np.asarray(getattr(norm, "boundaries"), dtype=np.float64)
        if b.size > 0:
            vmin = float(np.min(b))
    if vmin is None:
        return out
    out[m] = vmin
    return out


def _masked_invalid(data: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_invalid(np.asarray(data, dtype=np.float64))


def _cmap_with_transparent_bad(cmap_like):
    cmap = plt.get_cmap(cmap_like).copy()
    cmap.set_bad((0.0, 0.0, 0.0, 0.0))
    return cmap


def _discretize_cmap(cmap, bins: int, *, transparent_bad: bool) -> mcolors.Colormap:
    nbins = int(max(2, bins))
    colors = cmap(np.linspace(0.0, 1.0, nbins))
    out = mcolors.ListedColormap(colors, name=f"{cmap.name}_{nbins}")
    if transparent_bad:
        out.set_bad((0.0, 0.0, 0.0, 0.0))
    return out


def _visible_max_for_panel(
    data: np.ndarray,
    *,
    norm,
    vmax_override: float | None,
) -> float | None:
    # Report the true data maximum (finite values), independent of color scaling.
    arr = np.asarray(data, dtype=np.float64)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    return float(np.max(finite))


def _annotate_panel_max(ax, vmax: float | None) -> None:
    if vmax is None or not np.isfinite(vmax):
        return
    ax.text(
        0.98,
        0.02,
        f"max data: {vmax:.4g}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.7, "edgecolor": "0.5"},
    )


def _estimate_visible_vmax_via_render(
    ax,
    data: np.ndarray,
    *,
    extent: list[float],
    interpolation: str,
    transparent_nan: bool,
    vmin: float,
    vmax: float,
    log_scale: bool,
) -> float | None:
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        return None

    arr = np.asarray(data, dtype=np.float64)
    arr_plot = _masked_invalid(arr) if transparent_nan else arr
    gray_cmap = _cmap_with_transparent_bad("gray") if transparent_nan else plt.get_cmap("gray")

    tiny = np.finfo(np.float64).tiny
    norm = None
    vvmin = float(vmin)
    vvmax = float(vmax)
    if log_scale:
        vvmin = max(vvmin, tiny)
        vvmax = max(vvmax, vvmin * (1.0 + 1e-12))
        norm = mcolors.LogNorm(vmin=vvmin, vmax=vvmax, clip=True)

    tmp = ax.imshow(
        arr_plot,
        origin="lower",
        aspect=ax.get_aspect(),
        extent=extent,
        cmap=gray_cmap,
        norm=norm,
        vmin=None if norm is not None else vvmin,
        vmax=None if norm is not None else vvmax,
        interpolation=interpolation,
        interpolation_stage="rgba",
        resample=False,
        zorder=-1.0e9,
    )

    renderer = ax.figure.canvas.get_renderer()
    try:
        made = tmp.make_image(renderer, unsampled=False)
    finally:
        tmp.remove()

    rgba = made[0] if isinstance(made, tuple) else made
    rgba = np.asarray(rgba)
    if rgba.ndim != 3 or rgba.shape[2] < 4:
        return None

    alpha = rgba[..., 3]
    visible = alpha > 0
    if not np.any(visible):
        return None

    # Gray colormap: red channel tracks normalized intensity [0,1].
    tmax = float(np.max(rgba[..., 0][visible].astype(np.float64) / 255.0))
    tmax = float(np.clip(tmax, 0.0, 1.0))

    if log_scale:
        lvmin = np.log(vvmin)
        lvmax = np.log(vvmax)
        return float(np.exp(lvmin + tmax * (lvmax - lvmin)))

    return float(vvmin + tmax * (vvmax - vvmin))


def _tighten_panel_cbar_to_visible(
    ax,
    im,
    cb,
    data: np.ndarray,
    *,
    extent: list[float],
    interpolation: str,
    transparent_nan: bool,
    log_cbar_bins: int,
) -> None:
    if cb is None:
        return

    norm = im.norm
    log_scale = False
    bounds = None
    bins = int(max(2, log_cbar_bins))

    if isinstance(norm, mcolors.BoundaryNorm):
        b = np.asarray(getattr(norm, "boundaries", None), dtype=np.float64)
        if b.size < 2:
            return
        vmin_cur = float(np.min(b))
        vmax_cur = float(np.max(b))
        log_scale = True
        bounds = b
        bins = int(max(2, b.size - 1))
    elif isinstance(norm, mcolors.LogNorm):
        if norm.vmin is None or norm.vmax is None:
            return
        vmin_cur = float(norm.vmin)
        vmax_cur = float(norm.vmax)
        log_scale = True
    else:
        clim = im.get_clim()
        if clim is None:
            return
        vmin_cur = float(clim[0])
        vmax_cur = float(clim[1])

    vis_vmax = _estimate_visible_vmax_via_render(
        ax,
        data,
        extent=extent,
        interpolation=interpolation,
        transparent_nan=transparent_nan,
        vmin=vmin_cur,
        vmax=vmax_cur,
        log_scale=log_scale,
    )
    if vis_vmax is None or not np.isfinite(vis_vmax):
        return

    # Only tighten when there is a meaningful gap.
    rel_gap = (vmax_cur - vis_vmax) / max(abs(vmax_cur), np.finfo(np.float64).eps)
    if rel_gap <= 1e-3:
        return

    if log_scale:
        tiny = np.finfo(np.float64).tiny
        vmin_new = max(vmin_cur, tiny)
        vmax_new = max(vis_vmax, vmin_new * (1.0 + 1e-9))
        if vmax_new <= vmin_new:
            return

        if bounds is not None:
            new_norm, new_bounds = _build_bracketed_log_norm((vmin_new, vmax_new), bins=bins)
            if new_norm is None or new_bounds is None:
                return
            im.set_norm(new_norm)
            cb.boundaries = new_bounds
            cb.update_normal(im)
        else:
            im.set_norm(mcolors.LogNorm(vmin=vmin_new, vmax=vmax_new, clip=True))
            cb.update_normal(im)
        return

    eps = max(1e-12, abs(vmin_cur) * 1e-12)
    vmax_new = max(vis_vmax, vmin_cur + eps)
    if vmax_new <= vmin_cur:
        return
    im.set_clim(vmin_cur, vmax_new)
    cb.update_normal(im)


# -----------------------------
# Plotting
# -----------------------------

def plot_comparison_2x2(
    x: np.ndarray,
    y: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    sym_err: np.ndarray,
    hist_vals: np.ndarray,
    *,
    title_A: str,
    title_B: str,
    title_err: str,
    suptitle: str,
    out_path: str,
    title_abs: str = r"Absolute difference: $|A-B|$",
    dpi: int = 200,
    fig_width: float = 12.0,
    fig_height: float = 9.0,
    imshow_interpolation: str = "nearest",
    axes_aspect: str = "equal",
    top_cmap: str = "viridis",
    top_vmin: float | None = None,
    top_vmax: float | None = None,
    show_abs_diff: bool = True,
    log_lower_rows: bool = False,
    log_floor: float = 1e-12,
    log_cbar_bins: int = 16,
    lower_rows_shared_cbar: bool = False,
    annotate_panel_max: bool = True,
    transparent_nan: bool = True,
    info_text: str | None = None,
    extra_rows: list[dict[str, np.ndarray | str]] | None = None,
):
    """
    Base 2x2 comparison plot:
      [0,0] A
      [0,1] B
      [1,0] symmetric error map
      [1,1] absolute difference map |A-B| (or hidden when show_abs_diff=False)

    Optional extra rows:
      each item in extra_rows appends one row.
      Supported item formats:
        1) {"left_map","right_map","title_left","title_right"}
        2) {"sym_err","abs_diff","title_err","title_abs"}  (legacy alias)
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    extra_rows = [] if extra_rows is None else list(extra_rows)

    a_finite = A[np.isfinite(A)]
    b_finite = B[np.isfinite(B)]
    if a_finite.size == 0 or b_finite.size == 0:
        raise ValueError("plot_comparison_2x2: A and B must contain at least one finite value.")

    a_min = float(np.min(a_finite))
    b_min = float(np.min(b_finite))
    a_max = float(np.max(a_finite))
    b_max = float(np.max(b_finite))

    if (top_vmin is not None) or (top_vmax is not None):
        if top_vmin is None or top_vmax is None:
            raise ValueError("plot_comparison_2x2: both top_vmin and top_vmax must be set together.")
        vmin_top = float(top_vmin)
        vmax_top = float(top_vmax)
    else:
        # First-row colorbars: use the maximum of the dataset with the smaller maximum.
        # i.e. shared vmax = min(max(A), max(B)).
        vmin_top = min(a_min, b_min)
        vmax_top = min(a_max, b_max)
        if not np.isfinite(vmin_top) or not np.isfinite(vmax_top):
            raise ValueError("plot_comparison_2x2: invalid finite range for top-row color scale.")
        if vmin_top >= vmax_top:
            # Fallback for degenerate data ranges.
            vmin_top = min(a_min, b_min)
            vmax_top = max(a_max, b_max)

    nrows = 2 + len(extra_rows)
    fig_height_eff = fig_height + 2.2 * len(extra_rows)
    fig, axs = plt.subplots(nrows, 2, figsize=(fig_width, fig_height_eff), constrained_layout=True)
    extent = [x[0], x[-1], y[0], y[-1]]
    panel_artists: list[tuple[object, object, object, np.ndarray]] = []

    top_cmap_obj = _cmap_with_transparent_bad(top_cmap) if transparent_nan else plt.get_cmap(top_cmap)
    lower_cmap_base = _cmap_with_transparent_bad("viridis") if transparent_nan else plt.get_cmap("viridis")
    lower_cmap_obj = (
        _discretize_cmap(lower_cmap_base, log_cbar_bins, transparent_bad=transparent_nan)
        if log_lower_rows
        else lower_cmap_base
    )
    sym_decade_norm = sym_decade_bounds = sym_decade_ticks = sym_decade_labels = None
    if log_lower_rows:
        sym_cmap_bins = int(max(2, getattr(lower_cmap_obj, "N", log_cbar_bins)))
        sym_decade_norm, sym_decade_bounds, sym_decade_ticks, sym_decade_labels = _build_sym_err_decade_log_norm(
            floor=log_floor,
            cmap_bins=sym_cmap_bins,
        )

    A_plot = _masked_invalid(A) if transparent_nan else A
    B_plot = _masked_invalid(B) if transparent_nan else B
    im0 = axs[0, 0].imshow(
        A_plot,
        origin="lower",
        aspect=axes_aspect,
        extent=extent,
        vmin=vmin_top,
        vmax=vmax_top,
        cmap=top_cmap_obj,
        interpolation=imshow_interpolation,
        interpolation_stage="rgba",
        resample=False,
    )
    axs[0, 0].set_title(title_A)
    axs[0, 0].set_xlabel("x")
    axs[0, 0].set_ylabel("y")
    cb0 = fig.colorbar(im0, ax=axs[0, 0])
    panel_artists.append((axs[0, 0], im0, cb0, A))
    if annotate_panel_max:
        _annotate_panel_max(axs[0, 0], _visible_max_for_panel(A, norm=None, vmax_override=vmax_top))

    im1 = axs[0, 1].imshow(
        B_plot,
        origin="lower",
        aspect=axes_aspect,
        extent=extent,
        vmin=vmin_top,
        vmax=vmax_top,
        cmap=top_cmap_obj,
        interpolation=imshow_interpolation,
        interpolation_stage="rgba",
        resample=False,
    )
    axs[0, 1].set_title(title_B)
    axs[0, 1].set_xlabel("x")
    axs[0, 1].set_ylabel("y")
    cb1 = fig.colorbar(im1, ax=axs[0, 1])
    panel_artists.append((axs[0, 1], im1, cb1, B))
    if annotate_panel_max:
        _annotate_panel_max(axs[0, 1], _visible_max_for_panel(B, norm=None, vmax_override=vmax_top))

    # Pre-parse extra rows to compute lower-row color scales.
    # Default behavior is per-panel cbar range based on each map's visible finite data.
    # If lower_rows_shared_cbar=True, rows without `separate_cbar` share per-column scales.
    parsed_rows: list[tuple[np.ndarray, np.ndarray, str, str, bool]] = []
    for i, row in enumerate(extra_rows):
        if ("left_map" in row) and ("right_map" in row):
            left_map = np.asarray(row["left_map"], dtype=np.float64)
            right_map = np.asarray(row["right_map"], dtype=np.float64)
            title_left = str(row.get("title_left", "Left"))
            title_right = str(row.get("title_right", "Right"))
            separate_cbar = bool(row.get("separate_cbar", False))
        else:
            left_map = np.asarray(row["sym_err"], dtype=np.float64)
            right_map = np.asarray(row["abs_diff"], dtype=np.float64)
            title_left = str(row["title_err"])
            title_right = str(row.get("title_abs", r"Absolute difference: $|A-B|$"))
            separate_cbar = bool(row.get("separate_cbar", False))

        if left_map.shape != A.shape or right_map.shape != A.shape:
            raise ValueError(
                f"plot_comparison_2x2: extra row {i} shape mismatch: "
                f"left={left_map.shape}, right={right_map.shape}, expected {A.shape}"
            )
        parsed_rows.append((left_map, right_map, title_left, title_right, separate_cbar))

    base_sym = np.asarray(sym_err, dtype=np.float64)
    base_abs = np.abs(A - B)

    def _panel_scale(data_map: np.ndarray) -> tuple[mcolors.Normalize | None, np.ndarray | None, float | None, float | None]:
        arr = np.asarray(data_map, dtype=np.float64)
        log_limits = _shared_log_limits_for_maps([arr], enabled=log_lower_rows, floor=log_floor)
        norm, bounds = _build_bracketed_log_norm(log_limits, bins=log_cbar_bins)
        if norm is not None:
            return norm, bounds, None, None
        lin = _shared_linear_limits_for_maps([arr])
        if lin is None:
            return None, None, None, None
        return None, None, lin[0], lin[1]

    col1_log_norm = col1_log_bounds = None
    col2_log_norm = col2_log_bounds = None
    col1_lin = col2_lin = None
    if lower_rows_shared_cbar:
        col1_maps = [base_sym] + [r[0] for r in parsed_rows if not r[4]]
        col2_maps: list[np.ndarray] = []
        if show_abs_diff:
            col2_maps.append(base_abs)
        col2_maps.extend([r[1] for r in parsed_rows if not r[4]])

        col2_log_limits = _shared_log_limits_for_maps(col2_maps, enabled=log_lower_rows, floor=log_floor)
        if log_lower_rows:
            col1_log_norm, col1_log_bounds = sym_decade_norm, sym_decade_bounds
        else:
            col1_log_limits = _shared_log_limits_for_maps(col1_maps, enabled=log_lower_rows, floor=log_floor)
            col1_log_norm, col1_log_bounds = _build_bracketed_log_norm(col1_log_limits, bins=log_cbar_bins)
        col2_log_norm, col2_log_bounds = _build_bracketed_log_norm(col2_log_limits, bins=log_cbar_bins)
        col1_lin = _shared_linear_limits_for_maps(col1_maps)
        col2_lin = _shared_linear_limits_for_maps(col2_maps)

    if lower_rows_shared_cbar:
        sym_norm, sym_bounds = col1_log_norm, col1_log_bounds
        sym_vmin = None if col1_log_norm is not None or col1_lin is None else col1_lin[0]
        sym_vmax = None if col1_log_norm is not None or col1_lin is None else col1_lin[1]
    else:
        if log_lower_rows:
            sym_norm, sym_bounds, sym_vmin, sym_vmax = sym_decade_norm, sym_decade_bounds, None, None
        else:
            sym_norm, sym_bounds, sym_vmin, sym_vmax = _panel_scale(base_sym)

    sym_plot_raw = _prepare_map_for_log_norm(base_sym, sym_norm)
    sym_plot = _masked_invalid(sym_plot_raw) if transparent_nan else sym_plot_raw

    im2 = axs[1, 0].imshow(
        sym_plot,
        origin="lower",
        aspect=axes_aspect,
        extent=extent,
        norm=sym_norm,
        vmin=sym_vmin,
        vmax=sym_vmax,
        cmap=lower_cmap_obj,
        interpolation=imshow_interpolation,
        interpolation_stage="rgba",
        resample=False,
    )
    axs[1, 0].set_title(title_err)
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("y")
    cb2_kwargs = {"boundaries": sym_bounds} if sym_bounds is not None else {}
    if log_lower_rows and sym_decade_ticks is not None:
        cb2_kwargs["ticks"] = sym_decade_ticks
    cb2 = fig.colorbar(im2, ax=axs[1, 0], **cb2_kwargs)
    panel_artists.append((axs[1, 0], im2, cb2, base_sym))
    if log_lower_rows and sym_decade_labels is not None:
        cb2.set_ticklabels(sym_decade_labels)
    if annotate_panel_max:
        _annotate_panel_max(axs[1, 0], _visible_max_for_panel(base_sym, norm=sym_norm, vmax_override=sym_vmax))

    if show_abs_diff:
        if lower_rows_shared_cbar:
            abs_norm, abs_bounds = col2_log_norm, col2_log_bounds
            abs_vmin = None if col2_log_norm is not None or col2_lin is None else col2_lin[0]
            abs_vmax = None if col2_log_norm is not None or col2_lin is None else col2_lin[1]
        else:
            abs_norm, abs_bounds, abs_vmin, abs_vmax = _panel_scale(base_abs)

        abs_diff_raw = _prepare_map_for_log_norm(base_abs, abs_norm)
        abs_diff = _masked_invalid(abs_diff_raw) if transparent_nan else abs_diff_raw
        im3 = axs[1, 1].imshow(
            abs_diff,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            norm=abs_norm,
            vmin=abs_vmin,
            vmax=abs_vmax,
            cmap=lower_cmap_obj,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        axs[1, 1].set_title(title_abs)
        axs[1, 1].set_xlabel("x")
        axs[1, 1].set_ylabel("y")
        cb3_kwargs = {"boundaries": abs_bounds} if abs_bounds is not None else {}
        cb3 = fig.colorbar(im3, ax=axs[1, 1], **cb3_kwargs)
        panel_artists.append((axs[1, 1], im3, cb3, base_abs))
        if annotate_panel_max:
            _annotate_panel_max(axs[1, 1], _visible_max_for_panel(base_abs, norm=abs_norm, vmax_override=abs_vmax))
    else:
        axs[1, 1].axis("off")

    for i, (left_map, right_map, title_left, title_right, separate_cbar) in enumerate(parsed_rows):
        r = 2 + i
        if lower_rows_shared_cbar and not separate_cbar:
            left_norm, left_bounds = col1_log_norm, col1_log_bounds
            right_norm, right_bounds = col2_log_norm, col2_log_bounds
            left_vmin = None if col1_log_norm is not None or col1_lin is None else col1_lin[0]
            left_vmax = None if col1_log_norm is not None or col1_lin is None else col1_lin[1]
            right_vmin = None if col2_log_norm is not None or col2_lin is None else col2_lin[0]
            right_vmax = None if col2_log_norm is not None or col2_lin is None else col2_lin[1]
        else:
            if log_lower_rows:
                left_norm, left_bounds, left_vmin, left_vmax = sym_decade_norm, sym_decade_bounds, None, None
            else:
                left_norm, left_bounds, left_vmin, left_vmax = _panel_scale(left_map)
            right_norm, right_bounds, right_vmin, right_vmax = _panel_scale(right_map)

        left_plot_raw = _prepare_map_for_log_norm(left_map, left_norm)
        left_plot = _masked_invalid(left_plot_raw) if transparent_nan else left_plot_raw

        imL = axs[r, 0].imshow(
            left_plot,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            norm=left_norm,
            vmin=left_vmin,
            vmax=left_vmax,
            cmap=lower_cmap_obj,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        axs[r, 0].set_title(title_left)
        axs[r, 0].set_xlabel("x")
        axs[r, 0].set_ylabel("y")
        cbL_kwargs = {"boundaries": left_bounds} if left_bounds is not None else {}
        if log_lower_rows and sym_decade_ticks is not None:
            cbL_kwargs["ticks"] = sym_decade_ticks
        cbL = fig.colorbar(imL, ax=axs[r, 0], **cbL_kwargs)
        panel_artists.append((axs[r, 0], imL, cbL, left_map))
        if log_lower_rows and sym_decade_labels is not None:
            cbL.set_ticklabels(sym_decade_labels)
        if annotate_panel_max:
            _annotate_panel_max(axs[r, 0], _visible_max_for_panel(left_map, norm=left_norm, vmax_override=left_vmax))

        right_plot_raw = _prepare_map_for_log_norm(right_map, right_norm)
        right_plot = _masked_invalid(right_plot_raw) if transparent_nan else right_plot_raw

        imR = axs[r, 1].imshow(
            right_plot,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            norm=right_norm,
            vmin=right_vmin,
            vmax=right_vmax,
            cmap=lower_cmap_obj,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        axs[r, 1].set_title(title_right)
        axs[r, 1].set_xlabel("x")
        axs[r, 1].set_ylabel("y")
        cbR_kwargs = {"boundaries": right_bounds} if right_bounds is not None else {}
        cbR = fig.colorbar(imR, ax=axs[r, 1], **cbR_kwargs)
        panel_artists.append((axs[r, 1], imR, cbR, right_map))
        if annotate_panel_max:
            _annotate_panel_max(axs[r, 1], _visible_max_for_panel(right_map, norm=right_norm, vmax_override=right_vmax))

    # Tighten panel colorbars to the actual visible maxima after rasterization.
    fig.canvas.draw()
    for axp, imp, cbp, data_p in panel_artists:
        _tighten_panel_cbar_to_visible(
            axp,
            imp,
            cbp,
            data_p,
            extent=extent,
            interpolation=imshow_interpolation,
            transparent_nan=transparent_nan,
            log_cbar_bins=log_cbar_bins,
        )

    if info_text:
        fig.text(
            0.99,
            0.99,
            info_text,
            ha="right",
            va="top",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.5"},
        )

    fig.suptitle(suptitle)
    fig.savefig(out_path, dpi=int(dpi))
    plt.close(fig)


def plot_direction_1x3(
    x: np.ndarray,
    y: np.ndarray,
    dir_xem: np.ndarray,
    dir_com: np.ndarray,
    angle_err_deg: np.ndarray,
    *,
    suptitle: str,
    out_path: str,
    dpi: int = 200,
    fig_width: float = 16.0,
    fig_height: float = 5.5,
    imshow_interpolation: str = "nearest",
    axes_aspect: str = "equal",
    annotate_panel_max: bool = True,
    transparent_nan: bool = True,
    info_text: str | None = None,
) -> None:
    dir_xem = np.asarray(dir_xem, dtype=np.float64)
    dir_com = np.asarray(dir_com, dtype=np.float64)
    angle_err_deg = np.asarray(angle_err_deg, dtype=np.float64)
    if dir_xem.shape != dir_com.shape or dir_xem.shape != angle_err_deg.shape:
        raise ValueError(
            "plot_direction_1x3: all arrays must have same shape, got "
            f"{dir_xem.shape}, {dir_com.shape}, {angle_err_deg.shape}"
        )

    extent = [x[0], x[-1], y[0], y[-1]]
    top_cmap_obj = _cmap_with_transparent_bad("twilight_shifted") if transparent_nan else plt.get_cmap("twilight_shifted")
    err_cmap_obj = _cmap_with_transparent_bad("viridis") if transparent_nan else plt.get_cmap("viridis")

    fig, axs = plt.subplots(1, 3, figsize=(fig_width, fig_height), constrained_layout=True)
    panels = [
        (dir_xem, "XEMFEM (angle deg)", top_cmap_obj, -180.0, 180.0),
        (dir_com, "COMSOL (angle deg)", top_cmap_obj, -180.0, 180.0),
        (angle_err_deg, "Angle Err (deg)", err_cmap_obj, None, None),
    ]
    panel_artists: list[tuple[object, object, object, np.ndarray, int]] = []

    for ax, (arr, title, cmap_obj, vmin, vmax) in zip(axs, panels):
        arr_plot = _masked_invalid(arr) if transparent_nan else arr
        if vmin is None or vmax is None:
            lin = _shared_linear_limits_for_maps([arr])
            if lin is None:
                vmin_use = vmax_use = None
            else:
                vmin_use, vmax_use = lin
        else:
            vmin_use, vmax_use = vmin, vmax

        im = ax.imshow(
            arr_plot,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            cmap=cmap_obj,
            vmin=vmin_use,
            vmax=vmax_use,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cb = fig.colorbar(im, ax=ax)
        panel_artists.append((ax, im, cb, arr, 16))
        if annotate_panel_max:
            _annotate_panel_max(ax, _visible_max_for_panel(arr, norm=None, vmax_override=vmax_use))

    fig.canvas.draw()
    for axp, imp, cbp, data_p, bins in panel_artists:
        _tighten_panel_cbar_to_visible(
            axp,
            imp,
            cbp,
            data_p,
            extent=extent,
            interpolation=imshow_interpolation,
            transparent_nan=transparent_nan,
            log_cbar_bins=bins,
        )

    if info_text:
        fig.text(
            0.99,
            0.99,
            info_text,
            ha="right",
            va="top",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.5"},
        )

    fig.suptitle(suptitle)
    fig.savefig(out_path, dpi=int(dpi))
    plt.close(fig)


def save_plot_linear_and_log(
    *,
    out_path: str,
    **plot_kwargs,
) -> str:
    """
    Save two variants of the same comparison plot:
      - linear lower rows at out_path
      - log-scaled lower rows at out_path with '_log' suffix
    Returns the log-variant output path.
    """
    plot_comparison_2x2(out_path=out_path, log_lower_rows=False, **plot_kwargs)

    root, ext = os.path.splitext(out_path)
    if not ext:
        ext = ".png"
    out_log = f"{root}_log{ext}"

    plot_comparison_2x2(out_path=out_log, log_lower_rows=True, **plot_kwargs)
    return out_log


def plot_direction_comparison_2x3(
    x: np.ndarray,
    y: np.ndarray,
    dir_xem: np.ndarray,
    dir_com: np.ndarray,
    dir_com_alt: np.ndarray,
    err_xem_vs_com: np.ndarray,
    err_alt_vs_com: np.ndarray,
    err_xem_vs_alt: np.ndarray,
    *,
    suptitle: str,
    out_path: str,
    dpi: int = 200,
    fig_width: float = 16.0,
    fig_height: float = 9.0,
    imshow_interpolation: str = "nearest",
    axes_aspect: str = "equal",
    log_error_row: bool = False,
    log_floor: float = 1e-12,
    log_cbar_bins: int = 16,
    annotate_panel_max: bool = True,
    transparent_nan: bool = True,
    info_text: str | None = None,
) -> None:
    dir_xem = np.asarray(dir_xem, dtype=np.float64)
    dir_com = np.asarray(dir_com, dtype=np.float64)
    dir_com_alt = np.asarray(dir_com_alt, dtype=np.float64)
    err_xem_vs_com = np.asarray(err_xem_vs_com, dtype=np.float64)
    err_alt_vs_com = np.asarray(err_alt_vs_com, dtype=np.float64)
    err_xem_vs_alt = np.asarray(err_xem_vs_alt, dtype=np.float64)

    shape_ref = dir_xem.shape
    for name, arr in (
        ("dir_com", dir_com),
        ("dir_com_alt", dir_com_alt),
        ("err_xem_vs_com", err_xem_vs_com),
        ("err_alt_vs_com", err_alt_vs_com),
        ("err_xem_vs_alt", err_xem_vs_alt),
    ):
        if arr.shape != shape_ref:
            raise ValueError(
                f"plot_direction_comparison_2x3: shape mismatch for {name}: {arr.shape}, expected {shape_ref}"
            )

    top_cmap_obj = _cmap_with_transparent_bad("twilight_shifted") if transparent_nan else plt.get_cmap("twilight_shifted")
    err_cmap_base = _cmap_with_transparent_bad("viridis") if transparent_nan else plt.get_cmap("viridis")
    err_cmap_obj = (
        _discretize_cmap(err_cmap_base, log_cbar_bins, transparent_bad=transparent_nan)
        if log_error_row
        else err_cmap_base
    )

    fig, axs = plt.subplots(2, 3, figsize=(fig_width, fig_height), constrained_layout=True)
    extent = [x[0], x[-1], y[0], y[-1]]
    panel_artists: list[tuple[object, object, object, np.ndarray]] = []

    top_maps = [
        (dir_xem, "XEMFEM (angle deg)"),
        (dir_com, "COMSOL (angle deg)"),
        (dir_com_alt, "COMSOL (alt, XEMFEM mesh) (angle deg)"),
    ]
    for c, (arr, title) in enumerate(top_maps):
        arr_plot = _masked_invalid(arr) if transparent_nan else arr
        im = axs[0, c].imshow(
            arr_plot,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            cmap=top_cmap_obj,
            vmin=-180.0,
            vmax=180.0,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        axs[0, c].set_title(title)
        axs[0, c].set_xlabel("x")
        axs[0, c].set_ylabel("y")
        cb = fig.colorbar(im, ax=axs[0, c])
        panel_artists.append((axs[0, c], im, cb, arr))
        if annotate_panel_max:
            _annotate_panel_max(axs[0, c], _visible_max_for_panel(arr, norm=None, vmax_override=180.0))

    def _error_scale(arr: np.ndarray) -> tuple[mcolors.Normalize | None, np.ndarray | None, float | None, float | None]:
        log_limits = _shared_log_limits_for_maps([arr], enabled=log_error_row, floor=log_floor)
        norm, bounds = _build_bracketed_log_norm(log_limits, bins=log_cbar_bins)
        if norm is not None:
            return norm, bounds, None, None
        lin = _shared_linear_limits_for_maps([arr])
        if lin is None:
            return None, None, None, None
        return None, None, lin[0], lin[1]

    err_maps = [
        (err_xem_vs_com, "Angle Err: XEMFEM vs COMSOL (deg)"),
        (err_alt_vs_com, "Angle Err: COMSOL (alt) vs COMSOL (deg)"),
        (err_xem_vs_alt, "Angle Err: XEMFEM vs COMSOL (alt) (deg)"),
    ]
    for c, (arr, title) in enumerate(err_maps):
        norm, bounds, vmin, vmax = _error_scale(arr)
        arr_raw = _prepare_map_for_log_norm(arr, norm)
        arr_plot = _masked_invalid(arr_raw) if transparent_nan else arr_raw
        im = axs[1, c].imshow(
            arr_plot,
            origin="lower",
            aspect=axes_aspect,
            extent=extent,
            cmap=err_cmap_obj,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            interpolation=imshow_interpolation,
            interpolation_stage="rgba",
            resample=False,
        )
        axs[1, c].set_title(title)
        axs[1, c].set_xlabel("x")
        axs[1, c].set_ylabel("y")
        cb_kwargs = {"boundaries": bounds} if bounds is not None else {}
        cb = fig.colorbar(im, ax=axs[1, c], **cb_kwargs)
        panel_artists.append((axs[1, c], im, cb, arr))
        if annotate_panel_max:
            _annotate_panel_max(axs[1, c], _visible_max_for_panel(arr, norm=norm, vmax_override=vmax))

    fig.canvas.draw()
    for axp, imp, cbp, data_p in panel_artists:
        _tighten_panel_cbar_to_visible(
            axp,
            imp,
            cbp,
            data_p,
            extent=extent,
            interpolation=imshow_interpolation,
            transparent_nan=transparent_nan,
            log_cbar_bins=log_cbar_bins,
        )

    if info_text:
        fig.text(
            0.99,
            0.99,
            info_text,
            ha="right",
            va="top",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8, "edgecolor": "0.5"},
        )

    fig.suptitle(suptitle)
    fig.savefig(out_path, dpi=int(dpi))
    plt.close(fig)


def save_direction_plot_linear_and_log(
    *,
    out_path: str,
    **plot_kwargs,
) -> str:
    plot_direction_comparison_2x3(out_path=out_path, log_error_row=False, **plot_kwargs)

    root, ext = os.path.splitext(out_path)
    if not ext:
        ext = ".png"
    out_log = f"{root}_log{ext}"
    plot_direction_comparison_2x3(out_path=out_log, log_error_row=True, **plot_kwargs)
    return out_log


# -----------------------------
# VTU/PVD(U) comparison
# -----------------------------

def _xml_local_name(tag: str) -> str:
    if "}" in tag:
        return tag.split("}", 1)[1]
    return tag


def _parse_pvd_latest_dataset_file(pvd_path: str) -> str:
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    entries: list[tuple[float, int, str]] = []

    for i, elem in enumerate(root.iter()):
        if _xml_local_name(elem.tag) != "DataSet":
            continue
        rel = elem.attrib.get("file") or elem.attrib.get("File")
        if not rel:
            continue
        ts_raw = elem.attrib.get("timestep") or elem.attrib.get("time") or ""
        try:
            ts = float(ts_raw)
        except Exception:
            ts = float("-inf")
        entries.append((ts, i, rel))

    if not entries:
        raise ValueError(f"No DataSet entries found in {pvd_path}")

    _, _, rel_file = max(entries, key=lambda t: (t[0], t[1]))
    if os.path.isabs(rel_file):
        return rel_file
    return os.path.normpath(os.path.join(os.path.dirname(pvd_path), rel_file))


def _parse_parallel_piece_files(par_path: str) -> list[str]:
    tree = ET.parse(par_path)
    root = tree.getroot()
    pieces: list[str] = []

    for elem in root.iter():
        if _xml_local_name(elem.tag) != "Piece":
            continue
        rel = elem.attrib.get("Source") or elem.attrib.get("source")
        if not rel:
            continue
        if os.path.isabs(rel):
            pieces.append(rel)
        else:
            pieces.append(os.path.normpath(os.path.join(os.path.dirname(par_path), rel)))

    if not pieces:
        raise ValueError(f"No Piece sources found in {par_path}")
    return pieces


def _merge_piece_meshes(meshes: list["meshio.Mesh"]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if not meshes:
        raise ValueError("No VTU piece meshes to merge.")

    points_parts: list[np.ndarray] = []
    common_keys = set(meshes[0].point_data.keys())
    for m in meshes[1:]:
        common_keys &= set(m.point_data.keys())

    if not common_keys:
        raise ValueError("No common point-data fields across VTU pieces.")

    merged_data: dict[str, list[np.ndarray]] = {k: [] for k in sorted(common_keys)}

    for m in meshes:
        pts = np.asarray(m.points, dtype=np.float64)
        if pts.ndim != 2 or pts.shape[1] < 2:
            raise ValueError(f"Invalid points array shape in piece: {pts.shape}")
        points_parts.append(pts)

        npt = pts.shape[0]
        for k in merged_data.keys():
            arr = np.asarray(m.point_data[k])
            if arr.shape[0] != npt:
                raise ValueError(
                    f"Point-data '{k}' row count {arr.shape[0]} does not match points {npt}"
                )
            merged_data[k].append(arr)

    points = np.concatenate(points_parts, axis=0)
    point_data = {k: np.concatenate(v, axis=0) for k, v in merged_data.items()}
    return points, point_data


def _meshio_mesh_to_arrays(m: "meshio.Mesh") -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points = np.asarray(m.points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError(f"Invalid points array shape: {points.shape}")

    point_data = {k: np.asarray(v) for k, v in m.point_data.items()}
    return points, point_data


def _flatten_pyvista_datasets(obj, out: list) -> None:
    if obj is None:
        return

    if pv is not None and isinstance(obj, pv.MultiBlock):
        for i in range(obj.n_blocks):
            _flatten_pyvista_datasets(obj[i], out)
        return

    n_points = int(getattr(obj, "n_points", 0))
    if n_points > 0 and getattr(obj, "points", None) is not None:
        out.append(obj)


def _load_with_pyvista(path: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if pv is None:
        raise ImportError("pyvista is not available")

    ds = pv.read(path)
    blocks: list = []
    _flatten_pyvista_datasets(ds, blocks)
    if not blocks:
        raise ValueError(f"pyvista read produced no datasets for {path}")

    points_parts: list[np.ndarray] = []
    common_keys = set(blocks[0].point_data.keys())
    for b in blocks[1:]:
        common_keys &= set(b.point_data.keys())
    if not common_keys:
        raise ValueError("No common point-data fields across pyvista dataset blocks.")

    merged_data: dict[str, list[np.ndarray]] = {k: [] for k in sorted(common_keys)}
    for b in blocks:
        pts = np.asarray(b.points, dtype=np.float64)
        if pts.ndim != 2 or pts.shape[1] < 2:
            raise ValueError(f"Invalid pyvista points shape: {pts.shape}")
        points_parts.append(pts)
        npt = pts.shape[0]

        for k in merged_data.keys():
            arr = np.asarray(b.point_data[k])
            if arr.shape[0] != npt:
                raise ValueError(
                    f"pyvista point-data '{k}' row count {arr.shape[0]} does not match points {npt}"
                )
            merged_data[k].append(arr)

    points = np.concatenate(points_parts, axis=0)
    point_data = {k: np.concatenate(v, axis=0) for k, v in merged_data.items()}
    return points, point_data


def _load_with_meshio(path: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if meshio is None:
        raise ImportError("meshio is not available")

    m = meshio.read(path)
    return _meshio_mesh_to_arrays(m)


def _load_with_meshio_pieces(par_path: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if meshio is None:
        raise ImportError("meshio is not available")

    piece_files = _parse_parallel_piece_files(par_path)
    meshes = [meshio.read(pf) for pf in piece_files]
    return _merge_piece_meshes(meshes)


def _load_unstructured_dataset(path: str, _seen: set[str] | None = None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if _seen is None:
        _seen = set()

    apath = os.path.abspath(path)
    if apath in _seen:
        raise RuntimeError(f"Recursive VTU/PVD(U) reference detected at {apath}")
    _seen.add(apath)

    if not os.path.isfile(apath):
        raise FileNotFoundError(apath)

    ext = os.path.splitext(apath)[1].lower()
    if ext == ".pvd":
        next_path = _parse_pvd_latest_dataset_file(apath)
        return _load_unstructured_dataset(next_path, _seen=_seen)

    if ext in (".pvdu", ".pvtu"):
        if pv is not None:
            try:
                return _load_with_pyvista(apath)
            except Exception as pv_err:
                meshio_err: Exception | None = None
                try:
                    return _load_with_meshio_pieces(apath)
                except Exception as e:
                    meshio_err = e
                raise RuntimeError(
                    "Failed to load parallel VTK dataset with both backends. "
                    f"pyvista error: {pv_err}; meshio error: {meshio_err}"
                ) from pv_err

        meshio_err: Exception | None = None
        try:
            return _load_with_meshio_pieces(apath)
        except Exception as e:
            meshio_err = e

        raise RuntimeError(
            "Failed to load parallel VTK dataset with meshio (pyvista unavailable). "
            f"meshio error: {meshio_err}. "
            "Install pyvista (pip install pyvista) or use a meshio version supporting VTK version 2.2."
        ) from meshio_err

    if ext == ".vtu":
        if pv is not None:
            try:
                return _load_with_pyvista(apath)
            except Exception as pv_err:
                meshio_err: Exception | None = None
                try:
                    return _load_with_meshio(apath)
                except Exception as e:
                    meshio_err = e
                raise RuntimeError(
                    "Failed to load VTU dataset with both backends. "
                    f"pyvista error: {pv_err}; meshio error: {meshio_err}"
                ) from pv_err

        meshio_err: Exception | None = None
        try:
            return _load_with_meshio(apath)
        except Exception as e:
            meshio_err = e

        raise RuntimeError(
            "Failed to load VTU with meshio (pyvista unavailable). "
            f"meshio error: {meshio_err}. "
            "Install pyvista (pip install pyvista) or use a meshio version supporting VTK version 2.2."
        ) from meshio_err

    raise ValueError(f"Unsupported dataset extension '{ext}' for {apath}")


def _find_xemfem_simulation_dataset(case_dir: str) -> str | None:
    sim_dir = os.path.join(case_dir, "Simulation")
    if not os.path.isdir(sim_dir):
        return None

    preferred = ["Simulation.pvdu", "Simulation.pvtu", "Simulation.pvd", "Simulation.vtu"]
    for name in preferred:
        p = os.path.join(sim_dir, name)
        if os.path.isfile(p):
            return p

    for pattern in ("*.pvdu", "*.pvtu", "*.pvd", "*.vtu"):
        cands = sorted(glob.glob(os.path.join(sim_dir, pattern)))
        if cands:
            return cands[0]
    return None


def _canonical_field_name(name: str) -> str:
    return "".join(ch for ch in name.lower() if ch.isalnum())


def _pick_key_by_alias(
    data: dict[str, np.ndarray],
    aliases: tuple[str, ...],
) -> str | None:
    alias_index = {a: i for i, a in enumerate(aliases)}
    scored: list[tuple[int, int, str]] = []
    for key in sorted(data.keys()):
        canon = _canonical_field_name(key)
        if canon in alias_index:
            scored.append((alias_index[canon], len(key), key))
            continue
        for i, a in enumerate(aliases):
            if canon.startswith(a) or canon.endswith(a):
                scored.append((100 + i, len(key), key))
                break
    if not scored:
        return None
    scored.sort()
    return scored[0][2]


def _pick_point_field_keys(
    xemfem_data: dict[str, np.ndarray],
    comsol_data: dict[str, np.ndarray],
) -> tuple[str, str, str]:
    common = sorted(set(xemfem_data.keys()) & set(comsol_data.keys()))
    for key in ("V", "Emag", "E"):
        if key in common:
            label = "V" if key == "V" else ("|E|" if key == "Emag" else "E")
            return key, key, label

    # Semantic fallback for mismatched naming between tools, e.g.
    # XEMFEM: "V", COMSOL: "Electric_potential".
    groups: list[tuple[str, tuple[str, ...]]] = [
        (
            "V",
            (
                "v",
                "voltage",
                "potential",
                "electricpotential",
                "phi",
            ),
        ),
        (
            "|E|",
            (
                "emag",
                "electricfieldmagnitude",
                "norme",
                "abse",
                "fieldmagnitude",
            ),
        ),
        (
            "E",
            (
                "e",
                "electricfield",
                "efield",
                "electricfieldvector",
            ),
        ),
    ]
    for label, aliases in groups:
        kx = _pick_key_by_alias(xemfem_data, aliases)
        kc = _pick_key_by_alias(comsol_data, aliases)
        if kx is not None and kc is not None:
            return kx, kc, label

    if len(xemfem_data) == 1 and len(comsol_data) == 1:
        kx = next(iter(xemfem_data.keys()))
        kc = next(iter(comsol_data.keys()))
        return kx, kc, "field"

    raise ValueError(
        "Could not match point-data fields between XEMFEM and COMSOL. "
        f"XEMFEM keys={sorted(xemfem_data.keys())}, COMSOL keys={sorted(comsol_data.keys())}"
    )


def _scalarize_field_values(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        if arr.shape[1] == 1:
            return arr[:, 0]
        return np.linalg.norm(arr, axis=1)
    raise ValueError(f"Unsupported point-data shape for scalarization: {arr.shape}")


def _match_points_quantized(
    points_ref: np.ndarray,
    points_cmp: np.ndarray,
    *,
    coord_tol: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    pref = np.asarray(points_ref, dtype=np.float64)
    pcmp = np.asarray(points_cmp, dtype=np.float64)

    if pref.ndim != 2 or pcmp.ndim != 2:
        raise ValueError("Point arrays must be rank-2.")
    if pref.shape[0] == 0 or pcmp.shape[0] == 0:
        raise ValueError("Point arrays must be non-empty.")

    dim = min(pref.shape[1], pcmp.shape[1], 3)
    if dim < 1:
        raise ValueError("Point arrays must have at least one coordinate dimension.")

    # Interpret coord_tol as a relative (fractional) coordinate tolerance.
    # Example: coord_tol=1e-4 means 0.01% of coordinate scale.
    tol_rel = float(coord_tol)
    if tol_rel <= 0.0:
        tol_rel = 1e-9

    pref_d = pref[:, :dim]
    pcmp_d = pcmp[:, :dim]

    mn = np.minimum(np.min(pref_d, axis=0), np.min(pcmp_d, axis=0))
    mx = np.maximum(np.max(pref_d, axis=0), np.max(pcmp_d, axis=0))
    span = mx - mn

    absmax = np.maximum(np.max(np.abs(pref_d), axis=0), np.max(np.abs(pcmp_d), axis=0))
    scale = np.maximum(span, absmax)
    scale = np.where(np.isfinite(scale) & (scale > 0.0), scale, 1.0)

    tol_abs = tol_rel * scale
    tol_abs = np.where(np.isfinite(tol_abs) & (tol_abs > 0.0), tol_abs, 1e-12)

    # Quantize on a scale-normalized lattice so matching is based on percentage difference.
    qref = np.rint((pref_d - mn) / tol_abs).astype(np.int64)
    qcmp = np.rint((pcmp_d - mn) / tol_abs).astype(np.int64)

    key_to_ref: dict[tuple[int, ...], int] = {}
    for i in range(qref.shape[0]):
        key = tuple(qref[i, :].tolist())
        if key not in key_to_ref:
            key_to_ref[key] = i

    idx_ref: list[int] = []
    idx_cmp: list[int] = []
    for j in range(qcmp.shape[0]):
        key = tuple(qcmp[j, :].tolist())
        i = key_to_ref.get(key)
        if i is None:
            continue
        idx_ref.append(i)
        idx_cmp.append(j)

    iref = np.asarray(idx_ref, dtype=np.int64)
    icmp = np.asarray(idx_cmp, dtype=np.int64)
    match_frac = float(icmp.size) / float(pcmp.shape[0])
    return iref, icmp, match_frac


def plot_vtu_mean_deviation_scatter(
    comsol_vals: np.ndarray,
    xemfem_vals: np.ndarray,
    *,
    field_name: str,
    mean_abs_dev: float,
    mean_pct_dev: float,
    out_path: str,
    dpi: int,
) -> None:
    a = np.asarray(comsol_vals, dtype=np.float64)
    b = np.asarray(xemfem_vals, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"Scatter arrays must have same shape, got {a.shape} vs {b.shape}")

    finite = np.isfinite(a) & np.isfinite(b)
    a = a[finite]
    b = b[finite]
    if a.size == 0:
        raise ValueError("No finite values available for VTU comparison plot.")

    max_points = 250_000
    if a.size > max_points:
        step = int(np.ceil(a.size / max_points))
        a = a[::step]
        b = b[::step]

    lo = float(min(np.min(a), np.min(b)))
    hi = float(max(np.max(a), np.max(b)))
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise ValueError("Invalid finite range for VTU comparison plot.")
    if lo == hi:
        lo -= 1.0
        hi += 1.0

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 6.5), constrained_layout=True)
    ax.scatter(a, b, s=1.0, alpha=0.2, linewidths=0.0, rasterized=True)
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(f"COMSOL {field_name}")
    ax.set_ylabel(f"XEMFEM {field_name}")
    ax.set_title("COMSOL vs XEMFEM on simulation mesh")

    ax.text(
        0.98,
        0.98,
        f"mean |Δ| = {mean_abs_dev:.6g}\nmean %Δ = {mean_pct_dev:.6g}%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.5"},
    )

    fig.savefig(out_path, dpi=int(dpi))
    plt.close(fig)


def run_vtu_point_comparison(
    case_dir: str, *, coord_tol: float, dpi: int
) -> tuple[str, float, float, str] | None:
    comsol_vtu = os.path.join(case_dir, "COMSOL_XEMFEM.vtu")
    if not os.path.isfile(comsol_vtu):
        return None

    sim_dataset = _find_xemfem_simulation_dataset(case_dir)
    if sim_dataset is None:
        return None

    points_x, pdata_x = _load_unstructured_dataset(sim_dataset)
    points_c, pdata_c = _load_unstructured_dataset(comsol_vtu)

    key_x, key_c, field_label = _pick_point_field_keys(pdata_x, pdata_c)
    x_vals_all = _scalarize_field_values(pdata_x[key_x])
    c_vals_all = _scalarize_field_values(pdata_c[key_c])

    if x_vals_all.shape[0] != points_x.shape[0]:
        raise ValueError(
            f"XEMFEM field '{key_x}' length {x_vals_all.shape[0]} "
            f"does not match points {points_x.shape[0]}"
        )
    if c_vals_all.shape[0] != points_c.shape[0]:
        raise ValueError(
            f"COMSOL field '{key_c}' length {c_vals_all.shape[0]} "
            f"does not match points {points_c.shape[0]}"
        )

    iref, icmp, match_frac = _match_points_quantized(points_x, points_c, coord_tol=coord_tol)
    if iref.size == 0:
        raise ValueError("No overlapping points found between XEMFEM and COMSOL datasets.")
    match_frac_threshold = 0.95
    match_frac_ok = match_frac >= match_frac_threshold

    c_vals = c_vals_all[icmp]
    x_vals = x_vals_all[iref]
    valid = np.isfinite(c_vals) & np.isfinite(x_vals)
    c_vals = c_vals[valid]
    x_vals = x_vals[valid]
    if c_vals.size == 0:
        raise ValueError("No finite matched values to compare in VTU datasets.")

    abs_dev = np.abs(x_vals - c_vals)
    mean_abs_dev = float(np.mean(abs_dev))
    denom = np.maximum(np.abs(c_vals), np.finfo(np.float64).eps)
    mean_pct_dev = float(np.mean(abs_dev / denom) * 100.0)
    out_path = os.path.join(case_dir, "comparison_vtu_mean_deviation.png")
    plot_vtu_mean_deviation_scatter(
        c_vals,
        x_vals,
        field_name=field_label,
        mean_abs_dev=mean_abs_dev,
        mean_pct_dev=mean_pct_dev,
        out_path=out_path,
        dpi=dpi,
    )
    if not match_frac_ok:
        print(
            "  VTU compare warning:"
            f" only matched {match_frac * 100.0:.2f}% of COMSOL points "
            f"(threshold={100.0 * match_frac_threshold:.2f}%, relative tol={coord_tol:.3g}). "
            "Statistics are reported on the matched subset."
        )
    print(
        "  VTU compare:"
        f" field='{key_x}' vs '{key_c}' ({field_label}), matched={c_vals.size}/{points_c.shape[0]} "
        f"({100.0 * match_frac:.2f}%), mean|Δ|={mean_abs_dev:.6g}, mean%Δ={mean_pct_dev:.6g}%"
    )
    return out_path, mean_abs_dev, mean_pct_dev, field_label


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("parent_dir", help="Directory containing one subdirectory per validation case")
    ap.add_argument("--k", type=int, default=None, help="z-slice index (default: Nz//2)")
    ap.add_argument("--eps", type=float, default=1e-12, help="Epsilon for symmetric error denominator")
    ap.add_argument(
        "--coord-tol",
        type=float,
        default=1e-6,
        help="Absolute tolerance for grid coordinate comparisons and COMSOL grid binning.",
    )
    ap.add_argument(
        "--resample-xemfem-to-comsol",
        action="store_true",
        help="Resample XEMFEM fields onto the COMSOL native grid before comparison.",
    )
    ap.add_argument(
        "--plot-dpi",
        type=int,
        default=200,
        help="PNG output DPI for comparison plots.",
    )
    ap.add_argument(
        "--plot-width",
        type=float,
        default=12.0,
        help="Figure width in inches.",
    )
    ap.add_argument(
        "--plot-height",
        type=float,
        default=9.0,
        help="Figure height in inches.",
    )
    ap.add_argument(
        "--imshow-interpolation",
        type=str,
        default="nearest",
        choices=IMSHOW_INTERPOLATIONS,
        help="imshow interpolation mode (e.g. nearest, none, bilinear, bicubic).",
    )
    ap.add_argument(
        "--axes-aspect",
        type=str,
        default="equal",
        choices=("equal", "auto"),
        help="Data aspect for imshow axes. Use 'equal' to preserve x/y scale.",
    )
    args = ap.parse_args()

    parent_dir = os.path.abspath(args.parent_dir)
    if not os.path.isdir(parent_dir):
        raise NotADirectoryError(parent_dir)

    case_dirs = sorted(d for d in glob.glob(os.path.join(parent_dir, "*")) if os.path.isdir(d))
    if not case_dirs:
        raise RuntimeError(f"No subdirectories found in {parent_dir}")

    for case_dir in case_dirs:
        print(f"Processing {case_dir}")

        # ---- locate files ----
        txt_candidates = sorted(glob.glob(os.path.join(case_dir, "COMSOL*.txt")))
        if len(txt_candidates) == 0:
            print("  Skipping: no COMSOL*.txt found.")
            continue

        nastran_candidates = []
        for p in txt_candidates:
            b = os.path.basename(p).upper()
            if b.startswith("COMSOL") and b.endswith(".TXT") and ("NASTRAN" in b):
                nastran_candidates.append(p)
        if len(nastran_candidates) != 1:
            print(
                "  Skipping: expected exactly one COMSOL-on-XEMFEM-mesh txt (contains 'NASTRAN'), "
                f"found {len(nastran_candidates)}: {[os.path.basename(p) for p in nastran_candidates]}"
            )
            continue
        comsol_txt = nastran_candidates[0]

        h5_path = os.path.join(case_dir, "interpolated", "interpolated.h5")
        if not os.path.isfile(h5_path):
            print(f"  Skipping: missing {h5_path}")
            continue

        # ---- load XEMFEM (V and E directly) ----
        try:
            X2, Y2, V_xem = load_V_midplane(h5_path, k=args.k)
            _, _, Ex_xem, Ey_xem = load_E_midplane(h5_path, k=args.k)
        except Exception as e:
            print(f"  Skipping: XEMFEM load error: {e}")
            continue

        x_xem = X2[0, :]
        y_xem = Y2[:, 0]

        if Ex_xem.shape != V_xem.shape or Ey_xem.shape != V_xem.shape:
            print(f"  Skipping: XEMFEM shape mismatch: V {V_xem.shape}, Ex {Ex_xem.shape}, Ey {Ey_xem.shape}")
            continue

        # ---- load COMSOL-on-XEMFEM-mesh export (V, Ex, Ey) ----
        try:
            x_com, y_com, V_com_native, Ex_com_native, Ey_com_native, _ = load_comsol_native_regular_grid_V_Ex_Ey(
                comsol_txt,
                coord_tol=args.coord_tol,
            )
        except Exception as e:
            print(f"  Skipping: COMSOL load error: {e}")
            continue

        print_grid_alignment_report(
            x_xem, y_xem, x_com, y_com,
            coord_tol=args.coord_tol,
            label_ref="XEMFEM",
            label_cmp="COMSOL",
        )
        # Always compare on XEMFEM grid. Align COMSOL-on-XEMFEM export to this grid.
        x_plot = x_xem
        y_plot = y_xem
        V_x, Ex_x, Ey_x = V_xem, Ex_xem, Ey_xem
        try:
            V_ref = align_field_to_target_grid(
                V_com_native,
                x_com,
                y_com,
                x_plot,
                y_plot,
                coord_tol=args.coord_tol,
                label_src="COMSOL V",
                label_dst="XEMFEM",
            )
            Ex_ref = align_field_to_target_grid(
                Ex_com_native,
                x_com,
                y_com,
                x_plot,
                y_plot,
                coord_tol=args.coord_tol,
                label_src="COMSOL Ex",
                label_dst="XEMFEM",
            )
            Ey_ref = align_field_to_target_grid(
                Ey_com_native,
                x_com,
                y_com,
                x_plot,
                y_plot,
                coord_tol=args.coord_tol,
                label_src="COMSOL Ey",
                label_dst="XEMFEM",
            )
        except Exception as e:
            print(f"  Skipping: COMSOL->XEMFEM alignment error: {e}")
            continue

        if V_ref.shape != V_x.shape:
            print(f"  Skipping: shape mismatch V after alignment: {V_x.shape} vs {V_ref.shape}")
            continue
        if Ex_ref.shape != Ex_x.shape or Ey_ref.shape != Ey_x.shape:
            print("  Skipping: shape mismatch E after alignment")
            continue

        # ---- Voltage comparison ----
        sym_V = 100.0 * symmetric_error(V_x, V_ref, eps=args.eps)
        abs_V = np.abs(V_x - V_ref)
        vvals = sym_V[np.isfinite(sym_V)]
        if np.any(np.isfinite(sym_V)):
            iy_sym, ix_sym = np.unravel_index(np.nanargmax(sym_V), sym_V.shape)
            print(
                "  Voltage max sym err:"
                f" {sym_V[iy_sym, ix_sym]:.6g}% at"
                f" x={x_plot[ix_sym]:.6g}, y={y_plot[iy_sym]:.6g}"
            )
        if np.any(np.isfinite(abs_V)):
            iy_abs, ix_abs = np.unravel_index(np.nanargmax(abs_V), abs_V.shape)
            print(
                "  Voltage max abs diff:"
                f" {abs_V[iy_abs, ix_abs]:.6g} V at"
                f" x={x_plot[ix_abs]:.6g}, y={y_plot[iy_abs]:.6g}"
            )
        extra_rows_V: list[dict[str, np.ndarray | str]] = []
        voltage_out = os.path.join(case_dir, "comparison_voltage.png")
        plot_comparison_2x2(
            out_path=voltage_out,
            x=x_plot,
            y=y_plot,
            A=V_x,
            B=V_ref,
            sym_err=sym_V,
            hist_vals=vvals,
            title_A="XEMFEM (V)",
            title_B="COMSOL (V)",
            title_err="Sym Err (%): XEMFEM vs COMSOL",
            title_abs="Abs Diff: XEMFEM vs COMSOL",
            suptitle="Voltage",
            dpi=args.plot_dpi,
            fig_width=args.plot_width,
            fig_height=args.plot_height,
            imshow_interpolation=args.imshow_interpolation,
            axes_aspect=args.axes_aspect,
            extra_rows=extra_rows_V,
            log_lower_rows=False,
        )

        # ---- Electric field magnitude comparison (from DIRECT Ex/Ey) ----
        Emag_xem = E_magnitude(Ex_x, Ey_x)
        Emag_com = E_magnitude(Ex_ref, Ey_ref)
        sym_Emag = 100.0 * symmetric_error(Emag_xem, Emag_com, eps=args.eps)
        evals = sym_Emag[np.isfinite(sym_Emag)]
        extra_rows_E: list[dict[str, np.ndarray | str]] = []
        emag_out = os.path.join(case_dir, "comparison_E_magnitude.png")
        plot_comparison_2x2(
            out_path=emag_out,
            x=x_plot,
            y=y_plot,
            A=Emag_xem,
            B=Emag_com,
            sym_err=sym_Emag,
            hist_vals=evals,
            title_A="XEMFEM (|E|)",
            title_B="COMSOL (|E|)",
            title_err="Sym Err (%): XEMFEM vs COMSOL",
            title_abs="Abs Diff: XEMFEM vs COMSOL",
            suptitle="Electric field magnitude",
            dpi=args.plot_dpi,
            fig_width=args.plot_width,
            fig_height=args.plot_height,
            imshow_interpolation=args.imshow_interpolation,
            axes_aspect=args.axes_aspect,
            extra_rows=extra_rows_E,
            log_lower_rows=False,
        )

        # ---- Electric field direction comparison (from DIRECT Ex/Ey) ----
        dir_xem = E_direction_deg(Ex_x, Ey_x)
        dir_com = E_direction_deg(Ex_ref, Ey_ref)
        theta_err_deg = E_direction_mismatch_acos_deg(Ex_x, Ey_x, Ex_ref, Ey_ref, eps=args.eps)

        edir_out = os.path.join(case_dir, "comparison_E_direction.png")
        plot_direction_1x3(
            x=x_plot,
            y=y_plot,
            dir_xem=dir_xem,
            dir_com=dir_com,
            angle_err_deg=theta_err_deg,
            out_path=edir_out,
            suptitle="Electric field direction",
            dpi=args.plot_dpi,
            fig_width=max(args.plot_width, 16.0),
            fig_height=max(5.5, 0.6 * args.plot_height),
            imshow_interpolation=args.imshow_interpolation,
            axes_aspect=args.axes_aspect,
        )

        print("  Saved:")
        print(f"    {voltage_out}")
        print(f"    {emag_out}")
        print(f"    {edir_out}")


if __name__ == "__main__":
    main()
