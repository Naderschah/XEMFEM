# plot_per_angle.py
# Pure Python (no FreeCAD/Part) plotting per angle.
# Produces:
#   out_dir/all/all_<angle>.png
#   out_dir/by_bc/<BC>/<angle>.png
#   out_dir/by_component/<Component>/<angle>.png

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Pt2 = Tuple[float, float]


def iter_json_files(root: str) -> Iterable[str]:
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith(".json"):
                yield os.path.join(dirpath, fn)


def load_all_json(root_pts_out: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for p in iter_json_files(root_pts_out):
        try:
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
            data["_file"] = p
            out.append(data)
        except Exception:
            continue
    return out


def safe_name(s: str) -> str:
    s = s or "NA"
    return "".join(ch if ch.isalnum() or ch in "._-+" else "_" for ch in s).strip("._-") or "NA"


# ----------------------------
# Instruction decoding -> polyline
# ----------------------------

def _seg_start(instr: Any) -> Pt2:
    if isinstance(instr, list) and instr and isinstance(instr[0], str):
        return float(instr[1]), float(instr[2])
    return float(instr[0]), float(instr[1])


def _rot2(phi: float, x: float, y: float) -> Pt2:
    c = math.cos(phi)
    s = math.sin(phi)
    return (c * x - s * y, s * x + c * y)


def _inv_rot2(phi: float, x: float, y: float) -> Pt2:
    c = math.cos(phi)
    s = math.sin(phi)
    return (c * x + s * y, -s * x + c * y)


def _sample_circle_arc(p0: Pt2, p1: Pt2, c: Pt2, cw: bool, n: int = 64) -> List[Pt2]:
    x0, y0 = p0
    x1, y1 = p1
    cx, cy = c
    r = math.hypot(x0 - cx, y0 - cy)
    if r == 0:
        return [p0, p1]

    a0 = math.atan2(y0 - cy, x0 - cx)
    a1 = math.atan2(y1 - cy, x1 - cx)

    if cw:
        while a1 > a0:
            a1 -= 2 * math.pi
        ts = [a0 + (a1 - a0) * k / (n - 1) for k in range(n)]
    else:
        while a1 < a0:
            a1 += 2 * math.pi
        ts = [a0 + (a1 - a0) * k / (n - 1) for k in range(n)]

    return [(cx + r * math.cos(t), cy + r * math.sin(t)) for t in ts]


def _ellipse_param(p: Pt2, cx: float, cy: float, a: float, b: float, phi: float) -> float:
    x, y = p[0] - cx, p[1] - cy
    xl, yl = _inv_rot2(phi, x, y)
    aa = a if a != 0 else 1e-15
    bb = b if b != 0 else 1e-15
    return math.atan2(yl / bb, xl / aa)


def _sample_ellipse_arc(p0: Pt2, p1: Pt2, cx: float, cy: float, a: float, b: float, phi_deg: float, cw: bool, n: int = 96) -> List[Pt2]:
    phi = math.radians(phi_deg)
    t0 = _ellipse_param(p0, cx, cy, a, b, phi)
    t1 = _ellipse_param(p1, cx, cy, a, b, phi)

    if cw:
        while t1 > t0:
            t1 -= 2 * math.pi
        ts = [t0 + (t1 - t0) * k / (n - 1) for k in range(n)]
    else:
        while t1 < t0:
            t1 += 2 * math.pi
        ts = [t0 + (t1 - t0) * k / (n - 1) for k in range(n)]

    pts: List[Pt2] = []
    for t in ts:
        xl = a * math.cos(t)
        yl = b * math.sin(t)
        xr, yr = _rot2(phi, xl, yl)
        pts.append((cx + xr, cy + yr))
    return pts


def pts_list_to_polyline(pts: List[Any], samples_per_curve: int = 80) -> List[Pt2]:
    """
    Pure-python rendering.
    - bspline: control polygon only
    - hyperbola: chord only
    """
    if not pts:
        return []

    n = len(pts)
    out: List[Pt2] = []

    def stitch(seg: List[Pt2]) -> None:
        nonlocal out
        if not seg:
            return
        if out and math.hypot(out[-1][0] - seg[0][0], out[-1][1] - seg[0][1]) < 1e-12:
            out.extend(seg[1:])
        else:
            out.extend(seg)

    for i in range(n):
        instr = pts[i]
        instr_next = pts[(i + 1) % n]
        p0 = _seg_start(instr)
        p1 = _seg_start(instr_next)

        typ = "line"
        if isinstance(instr, list) and instr and isinstance(instr[0], str):
            typ = instr[0]

        if typ == "line":
            stitch([p0, p1])

        elif typ == "arc":
            cx, cy, cw = float(instr[3]), float(instr[4]), bool(instr[5])
            stitch(_sample_circle_arc(p0, p1, (cx, cy), cw, n=samples_per_curve))

        elif typ == "ellipse":
            cx, cy = float(instr[3]), float(instr[4])
            a, b = float(instr[5]), float(instr[6])
            phi_deg = float(instr[7])
            cw = bool(instr[8])
            stitch(_sample_ellipse_arc(p0, p1, cx, cy, a, b, phi_deg, cw, n=samples_per_curve))

        elif typ == "bspline":
            payload = instr[3] if len(instr) > 3 and isinstance(instr[3], dict) else {}
            poles = payload.get("poles", [])
            pole_pts: List[Pt2] = [(float(p[0]), float(p[1])) for p in poles if isinstance(p, list) and len(p) >= 2]
            stitch(pole_pts if len(pole_pts) >= 2 else [p0, p1])

        elif typ == "hyperbola":
            stitch([p0, p1])

        else:
            stitch([p0, p1])

    return out


# ----------------------------
# Grouping per angle
# ----------------------------

def collect_polylines_for_angle(data: Dict[str, Any], angle: float, samples_per_curve: int) -> List[List[Pt2]]:
    """
    If JSON contains multiple angles, include it if `angle` is in angles_deg.
    """
    angles = [float(a) for a in (data.get("angles_deg") or [])]
    if angle not in angles:
        return []

    polylines: List[List[Pt2]] = []

    if "instances" in data:
        for inst in (data.get("instances") or []):
            if inst.get("error") is not None:
                continue
            for pts in (inst.get("pts_lists") or []):
                poly = pts_list_to_polyline(pts, samples_per_curve=samples_per_curve)
                if len(poly) >= 2:
                    polylines.append(poly)
    else:
        for pts in (data.get("outers") or []):
            poly = pts_list_to_polyline(pts, samples_per_curve=samples_per_curve)
            if len(poly) >= 2:
                polylines.append(poly)

    return polylines


def unique_angles_present(all_data: List[Dict[str, Any]]) -> List[float]:
    s = set()
    for d in all_data:
        for a in (d.get("angles_deg") or []):
            try:
                s.add(float(a))
            except Exception:
                pass
    return sorted(s)


def plot_polylines(polylines: List[List[Pt2]], out_png: str, title: Optional[str] = None) -> None:
    plt.figure()
    if title:
        plt.title(title)

    any_drawn = False
    for poly in polylines:
        if len(poly) < 2:
            continue
        any_drawn = True
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        plt.plot(xs, ys, linewidth=0.8)

    if not any_drawn:
        plt.text(0.5, 0.5, "no curves", ha="center", va="center", transform=plt.gca().transAxes)

    plt.axis("equal")
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_per_angle_all_and_subsets(root_pts_out: str, out_dir: str, samples_per_curve: int = 80) -> None:
    all_data = load_all_json(root_pts_out)
    angles = unique_angles_present(all_data)

    out_all = os.path.join(out_dir, "all")
    out_bc = os.path.join(out_dir, "by_bc")
    out_comp = os.path.join(out_dir, "by_component")
    os.makedirs(out_all, exist_ok=True)
    os.makedirs(out_bc, exist_ok=True)
    os.makedirs(out_comp, exist_ok=True)

    # Pre-index by bc/component for faster grouping
    by_bc: Dict[str, List[Dict[str, Any]]] = {}
    by_comp: Dict[str, List[Dict[str, Any]]] = {}
    for d in all_data:
        bc = str(d.get("bc", "NA") or "NA")
        comp = str(d.get("component", "NA") or "NA")
        by_bc.setdefault(bc, []).append(d)
        by_comp.setdefault(comp, []).append(d)

    for ang in angles:
        ang_tag = f"{int(round(ang)):03d}"

        # all
        polys_all: List[List[Pt2]] = []
        for d in all_data:
            polys_all.extend(collect_polylines_for_angle(d, ang, samples_per_curve))
        plot_polylines(polys_all, os.path.join(out_all, f"all_{ang_tag}.png"), title=f"all @ {ang:g}°")

        # by bc
        for bc, ds in by_bc.items():
            polys: List[List[Pt2]] = []
            for d in ds:
                polys.extend(collect_polylines_for_angle(d, ang, samples_per_curve))
            if polys:
                plot_polylines(polys, os.path.join(out_bc, safe_name(bc), f"{ang_tag}.png"), title=f"BC={bc} @ {ang:g}°")

        # by component
        for comp, ds in by_comp.items():
            polys = []
            for d in ds:
                polys.extend(collect_polylines_for_angle(d, ang, samples_per_curve))
            if polys:
                plot_polylines(polys, os.path.join(out_comp, safe_name(comp), f"{ang_tag}.png"), title=f"{comp} @ {ang:g}°")


if __name__ == "__main__":
    ROOT_PTS_OUT = os.environ.get("PTS_OUT_DIR", "./out/pts_out")
    OUT_DIR = os.environ.get("PLOTS_OUT_DIR", "./out/plots_per_angle")
    plot_per_angle_all_and_subsets(ROOT_PTS_OUT, OUT_DIR, samples_per_curve=80)
