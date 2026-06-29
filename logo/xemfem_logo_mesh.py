#!/usr/bin/env python3
"""
Generate a 2-D Gmsh surface mesh for an XENONnT-style / XEMFEM logo.

The geometry is parametric rather than image-traced:
  - the upper logo is a grid of hexagonal cells, each split into 6 triangular surfaces;
  - the text is rasterized onto an equilateral triangular lattice;
  - two physical surface groups are created: BLUE and CYAN.

Usage:
  pip install gmsh
  python xemfem_logo_mesh.py --out xemfem_logo.msh --geo xemfem_logo.geo_unrolled --show
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Sequence

import gmsh

Point2 = tuple[float, float]
Segment = tuple[Point2, Point2]


BLUE = "BLUE"
CYAN = "CYAN"

# Hex row lengths for the upper logo. This is deliberately editable: it is the
# only part that encodes the current drawing-specific layout.
HEX_ROWS = [5, 9, 11, 11, 13, 13, 13, 13, 13, 11, 11, 9, 5]

# Approximate cyan cells in the upper logo, indexed as (row, column), where each
# row is independently centered. Tweak this set to match the final artwork.
CYAN_HEXES = {
    (3, 0), (3, 6), (3, 7), (3, 8), (3, 9), (3, 10),
    (4, 2), (4, 6),
    (5, 3), (5, 5),
    (6, 4), (6, 7), (6, 8), (6, 9), (6, 10),
    (7, 3), (7, 5),
    (8, 2), (8, 6),
    (9, 0), (9, 6), (9, 7), (9, 8), (9, 9), (9, 10),
}

# Simple stroke font in normalized glyph coordinates. y=0 is baseline, y=1 is cap height.
# The glyphs are rasterized onto a triangular lattice below.
GLYPH_STROKES: dict[str, list[Segment]] = {
    "X": [((0.00, 0.00), (1.00, 1.00)), ((0.00, 1.00), (1.00, 0.00))],
    "E": [((0.00, 0.00), (0.00, 1.00)), ((0.00, 1.00), (1.00, 1.00)),
          ((0.00, 0.50), (0.88, 0.50)), ((0.00, 0.00), (1.00, 0.00))],
    "M": [((0.00, 0.00), (0.00, 1.00)), ((1.00, 0.00), (1.00, 1.00)),
          ((0.00, 1.00), (0.50, 0.50)), ((0.50, 0.50), (1.00, 1.00))],
    "F": [((0.00, 0.00), (0.00, 1.00)), ((0.00, 1.00), (1.00, 1.00)),
          ((0.00, 0.50), (0.88, 0.50))],
}
GLYPH_WIDTHS = {"X": 1.00, "E": 0.95, "M": 1.15, "F": 0.95}


def round_key(x: float, ndigits: int = 10) -> float:
    """Stable keying for repeated geometric coordinates."""
    return round(float(x), ndigits)


@dataclass
class GmshSurfaceLogo:
    mesh_size: float

    def __post_init__(self) -> None:
        self.geo = gmsh.model.geo
        self.point_cache: dict[tuple[float, float], int] = {}
        # key -> (line_tag, start_point_tag, end_point_tag)
        self.line_cache: dict[frozenset[int], tuple[int, int, int]] = {}
        self.surfaces_by_group: dict[str, list[int]] = defaultdict(list)

    def point(self, xy: Point2) -> int:
        x, y = xy
        key = (round_key(x), round_key(y))
        if key not in self.point_cache:
            self.point_cache[key] = self.geo.addPoint(x, y, 0.0, self.mesh_size)
        return self.point_cache[key]

    def oriented_line(self, start: int, end: int) -> int:
        key = frozenset((start, end))
        if key not in self.line_cache:
            tag = self.geo.addLine(start, end)
            self.line_cache[key] = (tag, start, end)
            return tag
        tag, cached_start, cached_end = self.line_cache[key]
        if start == cached_start and end == cached_end:
            return tag
        return -tag

    def polygon_surface(self, vertices: Sequence[Point2]) -> int:
        if len(vertices) < 3:
            raise ValueError("A plane surface needs at least 3 vertices")
        pts = [self.point(p) for p in vertices]
        lines = [self.oriented_line(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
        loop = self.geo.addCurveLoop(lines)
        return self.geo.addPlaneSurface([loop])

    def add_triangle(self, a: Point2, b: Point2, c: Point2, group: str) -> int:
        s = self.polygon_surface([a, b, c])
        self.surfaces_by_group[group].append(s)
        return s

    def add_regular_hex(self, cx: float, cy: float, radius: float, group: str) -> list[int]:
        """Add a flat-top hexagon split into 6 triangular surfaces."""
        center = (cx, cy)
        vertices = [
            (cx + radius * math.cos(k * math.pi / 3.0),
             cy + radius * math.sin(k * math.pi / 3.0))
            for k in range(6)
        ]
        tags = []
        for k in range(6):
            # Each wedge is a true 2-D CAD surface. Shared spokes and edges are reused.
            tags.append(self.add_triangle(center, vertices[k], vertices[(k + 1) % 6], group))
        return tags


def point_segment_distance(p: Point2, a: Point2, b: Point2) -> float:
    px, py = p
    ax, ay = a
    bx, by = b
    vx, vy = bx - ax, by - ay
    wx, wy = px - ax, py - ay
    vv = vx * vx + vy * vy
    if vv == 0.0:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, (wx * vx + wy * vy) / vv))
    qx, qy = ax + t * vx, ay + t * vy
    return math.hypot(px - qx, py - qy)


def triangular_lattice_triangles(xmin: float, xmax: float, ymin: float, ymax: float, side: float) -> Iterable[tuple[Point2, Point2, Point2]]:
    """Yield equilateral triangles covering the requested bounding box."""
    h = math.sqrt(3.0) * side / 2.0
    x0 = xmin - 2.0 * side
    y0 = ymin - 2.0 * h
    n_j = int(math.ceil((ymax - y0) / h)) + 4

    def p(i: int, j: int) -> Point2:
        return (x0 + i * side + 0.5 * j * side, y0 + j * h)

    for j in range(n_j):
        # Pick a broad i range. The j-dependent skew is handled by the coordinate map p(i,j).
        i_min = int(math.floor((xmin - x0 - 0.5 * (j + 1) * side) / side)) - 3
        i_max = int(math.ceil((xmax - x0 - 0.5 * j * side) / side)) + 3
        for i in range(i_min, i_max + 1):
            # Two triangles per rhombus cell in the triangular lattice.
            t1 = (p(i, j), p(i + 1, j), p(i, j + 1))
            t2 = (p(i + 1, j), p(i + 1, j + 1), p(i, j + 1))
            yield t1
            yield t2


def transform_strokes(strokes: Sequence[Segment], x0: float, y0: float, width: float, height: float) -> list[Segment]:
    out = []
    for (ax, ay), (bx, by) in strokes:
        out.append(((x0 + ax * width, y0 + ay * height),
                    (x0 + bx * width, y0 + by * height)))
    return out


def add_triangular_stroke_glyph(
    logo: GmshSurfaceLogo,
    strokes: Sequence[Segment],
    x0: float,
    y0: float,
    width: float,
    height: float,
    triangle_side: float,
    stroke_width: float,
    group: str,
) -> None:
    scaled_strokes = transform_strokes(strokes, x0, y0, width, height)
    xmin, xmax = x0, x0 + width
    ymin, ymax = y0, y0 + height
    pad = 0.5 * stroke_width + triangle_side

    for tri in triangular_lattice_triangles(xmin - pad, xmax + pad, ymin - pad, ymax + pad, triangle_side):
        cx = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
        cy = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
        if not (xmin - pad <= cx <= xmax + pad and ymin - pad <= cy <= ymax + pad):
            continue
        d = min(point_segment_distance((cx, cy), a, b) for a, b in scaled_strokes)
        if d <= 0.5 * stroke_width:
            # Use a consistent orientation for positive area triangles.
            area2 = ((tri[1][0] - tri[0][0]) * (tri[2][1] - tri[0][1])
                     - (tri[1][1] - tri[0][1]) * (tri[2][0] - tri[0][0]))
            if area2 < 0:
                logo.add_triangle(tri[0], tri[2], tri[1], group)
            else:
                logo.add_triangle(tri[0], tri[1], tri[2], group)


def text_width(text: str, height: float, spacing: float) -> float:
    width = 0.0
    for k, ch in enumerate(text):
        width += GLYPH_WIDTHS[ch] * height
        if k != len(text) - 1:
            width += spacing
    return width


def add_text(
    logo: GmshSurfaceLogo,
    text: str = "XEMFEM",
    center_x: float = 0.0,
    baseline_y: float = -5.0,
    height: float = 1.2,
    triangle_side: float = 0.095,
    stroke_width: float = 0.18,
    spacing: float = 0.23,
    group: str = BLUE,
) -> None:
    total_w = text_width(text, height, spacing)
    x = center_x - 0.5 * total_w
    for ch in text:
        if ch not in GLYPH_STROKES:
            raise ValueError(f"Unsupported glyph {ch!r}; add it to GLYPH_STROKES")
        width = GLYPH_WIDTHS[ch] * height
        add_triangular_stroke_glyph(
            logo=logo,
            strokes=GLYPH_STROKES[ch],
            x0=x,
            y0=baseline_y,
            width=width,
            height=height,
            triangle_side=triangle_side,
            stroke_width=stroke_width,
            group=group,
        )
        x += width + spacing


def add_detector(
    logo: GmshSurfaceLogo,
    center_x: float = 0.0,
    center_y: float = 1.45,
    hex_radius: float = 0.23,
    pitch_x: float = 0.66,
    pitch_y: float = 0.58,
    cyan_radius_scale: float = 1.14,
) -> None:
    top_y = center_y + 0.5 * (len(HEX_ROWS) - 1) * pitch_y
    for row, ncols in enumerate(HEX_ROWS):
        y = top_y - row * pitch_y
        left_x = center_x - 0.5 * (ncols - 1) * pitch_x
        for col in range(ncols):
            group = CYAN if (row, col) in CYAN_HEXES else BLUE
            radius = hex_radius * cyan_radius_scale if group == CYAN else hex_radius
            logo.add_regular_hex(left_x + col * pitch_x, y, radius, group)


def add_physical_groups_and_colors(logo: GmshSurfaceLogo) -> None:
    gmsh.model.geo.synchronize()

    color_map = {
        BLUE: (46, 91, 170, 255),
        CYAN: (91, 207, 247, 255),
    }
    physical_tags = {BLUE: 1, CYAN: 2}

    for group, surfaces in logo.surfaces_by_group.items():
        tag = physical_tags[group]
        gmsh.model.addPhysicalGroup(2, surfaces, tag=tag)
        gmsh.model.setPhysicalName(2, tag, group)
        gmsh.model.setColor([(2, s) for s in surfaces], *color_map[group])


def build_model(args: argparse.Namespace) -> None:
    gmsh.initialize()
    gmsh.model.add("xemfem_logo")

    try:
        gmsh.option.setNumber("Mesh.Algorithm", args.algorithm)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.option.setNumber("Mesh.SaveAll", 0)

        logo = GmshSurfaceLogo(mesh_size=args.mesh_size)
        add_detector(
            logo,
            center_x=0.0,
            center_y=args.detector_center_y,
            hex_radius=args.hex_radius,
            pitch_x=args.hex_pitch_x,
            pitch_y=args.hex_pitch_y,
            cyan_radius_scale=args.cyan_radius_scale,
        )
        add_text(
            logo,
            text="MFEM",
            center_x=0.0,
            baseline_y=args.text_baseline_y,
            height=args.text_height,
            triangle_side=args.text_triangle_side,
            stroke_width=args.text_stroke_width,
            spacing=args.text_spacing,
            group=BLUE,
        )
        add_physical_groups_and_colors(logo)

        gmsh.model.mesh.generate(2)
        gmsh.write(args.out)
        if args.geo:
            gmsh.write(args.geo)
        if args.show:
            gmsh.fltk.run()
    finally:
        gmsh.finalize()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate a 2-D Gmsh mesh for the XEMFEM logo.")
    p.add_argument("--out", default="xemfem_logo.msh", help="Output Gmsh mesh file")
    p.add_argument("--geo", default="xemfem_logo.geo_unrolled", help="Optional unrolled geometry export; pass '' to skip")
    p.add_argument("--show", action="store_true", help="Open the Gmsh GUI after meshing")
    p.add_argument("--mesh-size", type=float, default=0.055, help="Target element size for generated surfaces")
    p.add_argument("--algorithm", type=int, default=6, help="Gmsh 2-D meshing algorithm, e.g. 5=Delaunay, 6=Frontal-Delaunay")

    # Upper logo parameters.
    p.add_argument("--detector-center-y", type=float, default=1.45)
    p.add_argument("--hex-radius", type=float, default=0.23)
    p.add_argument("--hex-pitch-x", type=float, default=0.56)
    p.add_argument("--hex-pitch-y", type=float, default=0.53)
    p.add_argument(
    "--cyan-radius-scale",
    type=float,
    default=1.14,
    help="Scale factor applied only to cyan detector hexes",
    )

    # Text parameters.
    p.add_argument("--text-baseline-y", type=float, default=-3.55)
    p.add_argument("--text-height", type=float, default=1.10)
    p.add_argument("--text-triangle-side", type=float, default=0.085)
    p.add_argument("--text-stroke-width", type=float, default=0.18)
    p.add_argument("--text-spacing", type=float, default=0.522)
    return p.parse_args()


if __name__ == "__main__":
    build_model(parse_args())
