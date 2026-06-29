#!/usr/bin/env python3
"""
Render a high-quality front view of the XEMFEM Gmsh logo mesh.

Default behavior draws the CAD/surface boundaries from the .msh file, so the
output looks like clean logo line art. Use --line-source mesh to draw every
finite-element edge instead.

Examples
--------
python render_xemfem_front.py xemfem_logo.msh --out xemfem_logo_front.svg
python render_xemfem_front.py xemfem_logo.msh --out xemfem_logo_front.png --dpi 900
python render_xemfem_front.py xemfem_logo.msh --out xemfem_logo_mesh.png --line-source mesh --fill-alpha 0.035
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


def _rounded_point_key(p: tuple[float, float], ndigits: int = 10) -> tuple[float, float]:
    return (round(float(p[0]), ndigits), round(float(p[1]), ndigits))


def _segment_key(a: tuple[float, float], b: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
    ka = _rounded_point_key(a)
    kb = _rounded_point_key(b)
    return tuple(sorted((ka, kb)))  # type: ignore[return-value]


def _default_color(name: str, fallback_index: int = 0) -> str:
    colors = {
        "BLUE": "#2E5BAA",
        "CYAN": "#5BCFF7",
    }
    if name.upper() in colors:
        return colors[name.upper()]

    fallback = [
        "#2E5BAA", "#5BCFF7", "#5A5A5A", "#9367BC",
        "#D8913C", "#409E6A", "#AA4A44", "#345A8A",
    ]
    return fallback[fallback_index % len(fallback)]


def _load_physical_surface_groups(gmsh_module):
    """Return [(physical_tag, physical_name, [surface_entity_tags]), ...]."""
    gmsh = gmsh_module
    physical_groups = gmsh.model.getPhysicalGroups(2)
    groups = []

    for _, physical_tag in physical_groups:
        name = gmsh.model.getPhysicalName(2, physical_tag) or f"physical_{physical_tag}"
        entities = list(gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag))
        groups.append((int(physical_tag), name, [int(e) for e in entities]))

    # Fallback for meshes that do not contain physical groups.
    if not groups:
        entities = [int(tag) for dim, tag in gmsh.model.getEntities(2) if dim == 2]
        groups.append((0, "BLUE", entities))

    return groups


def _node_coordinate_map(gmsh_module) -> dict[int, tuple[float, float]]:
    gmsh = gmsh_module
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords_by_tag: dict[int, tuple[float, float]] = {}
    for i, tag in enumerate(node_tags):
        coords_by_tag[int(tag)] = (float(coords[3 * i]), float(coords[3 * i + 1]))
    return coords_by_tag


def _surface_polygons_from_mesh(gmsh_module, surface_entities: Iterable[int]) -> list[list[tuple[float, float]]]:
    """Return first-order triangle/quad polygons from 2-D elements on the surfaces."""
    gmsh = gmsh_module
    coords_by_tag = _node_coordinate_map(gmsh)
    polygons: list[list[tuple[float, float]]] = []

    for surface_tag in surface_entities:
        try:
            element_types, _, element_node_tags = gmsh.model.mesh.getElements(2, int(surface_tag))
        except Exception:
            continue

        for element_type, nodes_flat in zip(element_types, element_node_tags):
            name, dim, _order, num_nodes, _local_coords, num_primary_nodes = gmsh.model.mesh.getElementProperties(int(element_type))
            if int(dim) != 2:
                continue

            # Use primary corner nodes for high-order elements, so curved/high-order
            # meshes still render as clean first-order polygons.
            n_total = int(num_nodes)
            n_primary = int(num_primary_nodes)
            if n_primary not in (3, 4):
                # Ignore non triangle/quad surface elements for this front-view renderer.
                continue

            node_list = [int(n) for n in nodes_flat]
            for offset in range(0, len(node_list), n_total):
                corner_tags = node_list[offset: offset + n_primary]
                try:
                    polygons.append([coords_by_tag[t] for t in corner_tags])
                except KeyError:
                    pass

    return polygons


def _segments_from_mesh_polygons(polygons: Iterable[list[tuple[float, float]]]) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    segments_by_key: dict[tuple[tuple[float, float], tuple[float, float]], tuple[tuple[float, float], tuple[float, float]]] = {}
    for polygon in polygons:
        if len(polygon) < 2:
            continue
        for i in range(len(polygon)):
            a = polygon[i]
            b = polygon[(i + 1) % len(polygon)]
            segments_by_key[_segment_key(a, b)] = (a, b)
    return list(segments_by_key.values())


def _point_xyz(gmsh_module, point_tag: int) -> tuple[float, float]:
    gmsh = gmsh_module
    xyz = gmsh.model.getValue(0, int(point_tag), [])
    return (float(xyz[0]), float(xyz[1]))


def _segments_from_cad_boundaries(gmsh_module, surface_entities: Iterable[int]) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Return unique straight boundary segments of CAD surface entities."""
    gmsh = gmsh_module
    segments_by_key: dict[tuple[tuple[float, float], tuple[float, float]], tuple[tuple[float, float], tuple[float, float]]] = {}
    curve_tags: set[int] = set()

    for surface_tag in surface_entities:
        try:
            boundary = gmsh.model.getBoundary([(2, int(surface_tag))], oriented=False, recursive=False)
        except Exception:
            continue
        for dim, tag in boundary:
            if int(dim) == 1:
                curve_tags.add(abs(int(tag)))

    for curve_tag in curve_tags:
        try:
            endpoints = gmsh.model.getBoundary([(1, curve_tag)], oriented=False, recursive=False)
        except Exception:
            continue

        point_tags = [abs(int(tag)) for dim, tag in endpoints if int(dim) == 0]
        if len(point_tags) < 2:
            continue

        a = _point_xyz(gmsh, point_tags[0])
        b = _point_xyz(gmsh, point_tags[1])
        segments_by_key[_segment_key(a, b)] = (a, b)

    return list(segments_by_key.values())


def _bounds_from_geometry(groups) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for data in groups.values():
        for poly in data.get("polygons", []):
            for x, y in poly:
                xs.append(float(x))
                ys.append(float(y))
        for a, b in data.get("segments", []):
            xs.extend([float(a[0]), float(b[0])])
            ys.extend([float(a[1]), float(b[1])])

    if not xs or not ys:
        raise RuntimeError("No drawable geometry was found in the mesh.")

    return min(xs), max(xs), min(ys), max(ys)


def render(args: argparse.Namespace) -> None:
    import gmsh
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    gmsh.initialize()
    try:
        gmsh.open(str(args.mesh))
        physical_groups = _load_physical_surface_groups(gmsh)

        drawable = {}
        for i, (_physical_tag, name, surface_entities) in enumerate(physical_groups):
            polygons = _surface_polygons_from_mesh(gmsh, surface_entities)
            if args.line_source == "cad":
                segments = _segments_from_cad_boundaries(gmsh, surface_entities)
                if not segments:
                    segments = _segments_from_mesh_polygons(polygons)
            else:
                segments = _segments_from_mesh_polygons(polygons)

            drawable[name] = {
                "color": args.color.get(name.upper(), _default_color(name, i)),
                "polygons": polygons,
                "segments": segments,
            }
    finally:
        gmsh.finalize()

    xmin, xmax, ymin, ymax = _bounds_from_geometry(drawable)
    width = xmax - xmin
    height = ymax - ymin
    pad = args.padding * max(width, height)
    xmin -= pad
    xmax += pad
    ymin -= pad
    ymax += pad
    width = xmax - xmin
    height = ymax - ymin

    fig_width = float(args.width_in)
    fig_height = fig_width * height / width

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=args.dpi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")

    if not args.transparent:
        fig.patch.set_facecolor(args.background)
        ax.set_facecolor(args.background)

    # Optional low-alpha face fill, useful for showing physical groups without
    # losing the line-art character of the logo.
    if args.fill_alpha > 0.0:
        for name, data in drawable.items():
            if data["polygons"]:
                pc = PolyCollection(
                    data["polygons"],
                    facecolors=data["color"],
                    edgecolors="none",
                    alpha=args.fill_alpha,
                    antialiased=True,
                )
                ax.add_collection(pc)

    # Draw darker/lighter under-stroke first if requested. This can make PNGs
    # look clean on non-white backgrounds.
    if args.underlay_width > 0.0:
        for _name, data in drawable.items():
            if data["segments"]:
                under = LineCollection(
                    data["segments"],
                    colors=args.underlay_color,
                    linewidths=args.line_width + args.underlay_width,
                    antialiased=True,
                )
                ax.add_collection(under)

    for _name, data in drawable.items():
        if data["segments"]:
            lc = LineCollection(
                data["segments"],
                colors=data["color"],
                linewidths=args.line_width,
                antialiased=True,
            )
            ax.add_collection(lc)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out,
        dpi=args.dpi,
        bbox_inches="tight",
        pad_inches=0.0,
        transparent=args.transparent,
        facecolor=("none" if args.transparent else args.background),
    )
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render a clean front view from a Gmsh .msh logo file.")
    p.add_argument("mesh", type=Path, help="Input .msh file")
    p.add_argument("--out", default="xemfem_logo_front.svg", help="Output image: .svg, .pdf, .png, etc.")
    p.add_argument("--line-source", choices=("cad", "mesh"), default="cad", help="cad = clean logo surface boundaries; mesh = every FE element edge")
    p.add_argument("--dpi", type=int, default=900, help="DPI for raster outputs; ignored by true vector backends like SVG/PDF")
    p.add_argument("--width-in", type=float, default=9.0, help="Figure width in inches before rasterization")
    p.add_argument("--padding", type=float, default=0.035, help="Relative padding around the logo")
    p.add_argument("--line-width", type=float, default=1.25, help="Stroke width in points")
    p.add_argument("--fill-alpha", type=float, default=0.0, help="Optional transparent fill for surfaces")
    p.add_argument("--transparent", action="store_true", help="Transparent background")
    p.add_argument("--background", default="#FFFFFF", help="Background color when not transparent")
    p.add_argument("--blue", default="#2E5BAA", help="Color for physical group BLUE")
    p.add_argument("--cyan", default="#5BCFF7", help="Color for physical group CYAN")
    p.add_argument("--underlay-width", type=float, default=0.0, help="Extra width of optional under-stroke in points")
    p.add_argument("--underlay-color", default="#FFFFFF", help="Color of optional under-stroke")

    args = p.parse_args()
    args.color = {
        "BLUE": args.blue,
        "CYAN": args.cyan,
    }
    return args


if __name__ == "__main__":
    render(parse_args())
