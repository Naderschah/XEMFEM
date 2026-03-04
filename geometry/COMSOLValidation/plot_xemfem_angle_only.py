#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np


def load_E_midplane(h5_path: str, k: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as h5:
        if "/E/field" not in h5:
            raise KeyError("Missing dataset /E/field")

        E = h5["/E/field"][...]

        if E.ndim == 4:
            Nz, Ny, Nx, dim = E.shape
            if dim < 2:
                raise ValueError(f"/E/field last dimension must be >=2, got {dim}")
            kk = Nz // 2 if k is None else int(k)
            if kk < 0 or kk >= Nz:
                raise ValueError(f"k={kk} out of range [0,{Nz-1}]")
            Ex = np.asarray(E[kk, :, :, 0], dtype=np.float64)
            Ey = np.asarray(E[kk, :, :, 1], dtype=np.float64)
            return Ex, Ey

        if E.ndim == 3:
            a, b, c = E.shape
            if a == 2:
                Ex = np.asarray(E[0, :, :], dtype=np.float64)
                Ey = np.asarray(E[1, :, :], dtype=np.float64)
                return Ex, Ey
            if c == 2:
                Ex = np.asarray(E[:, :, 0], dtype=np.float64)
                Ey = np.asarray(E[:, :, 1], dtype=np.float64)
                return Ex, Ey

        raise ValueError(f"Unsupported /E/field shape: {E.shape}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Temporary plot: XEMFEM angle(deg) only, no axes/padding/cbar, pixel-accurate."
    )
    ap.add_argument("h5", help="Path to interpolated/interpolated.h5")
    ap.add_argument(
        "-o",
        "--out",
        default=None,
        help="Output PNG path (default: <h5_dir>/xemfem_angle_only.png)",
    )
    ap.add_argument("--k", type=int, default=None, help="z-slice index for 3D datasets (default: Nz//2)")
    ap.add_argument(
        "--require-8k",
        action="store_true",
        help="Fail unless data is exactly 7680x4320 (Nx x Ny).",
    )
    args = ap.parse_args()

    h5_path = os.path.abspath(args.h5)
    if not os.path.isfile(h5_path):
        raise FileNotFoundError(h5_path)

    Ex, Ey = load_E_midplane(h5_path, k=args.k)
    angle_deg = np.degrees(np.arctan2(Ey, Ex))

    Ny, Nx = angle_deg.shape
    if args.require_8k and (Nx != 7680 or Ny != 4320):
        raise ValueError(f"Expected 7680x4320 (Nx x Ny), got {Nx}x{Ny}")

    out_path = args.out
    if out_path is None:
        out_path = os.path.join(os.path.dirname(h5_path), "xemfem_angle_only.png")
    out_path = os.path.abspath(out_path)

    # Pixel-accurate write: output image size exactly equals data shape (Ny x Nx).
    # NaN/invalid values are rendered as solid black.
    cmap = plt.get_cmap("twilight_shifted").copy()
    cmap.set_bad((0.0, 0.0, 0.0, 1.0))
    angle_plot = np.ma.masked_invalid(angle_deg)

    plt.imsave(
        out_path,
        angle_plot,
        cmap=cmap,
        vmin=-180.0,
        vmax=180.0,
        origin="lower",
    )

    print(f"[OK] wrote {out_path} ({Nx}x{Ny} px)")


if __name__ == "__main__":
    main()
