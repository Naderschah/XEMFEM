#!/usr/bin/env python3
"""
Convert COMSOLValidation meshes (e.g., .msh/.med/.mesh) to NASTRAN BDF using case config files.

Usage examples:
  python3 createCOMSOLmesh.py
  python3 createCOMSOLmesh.py -c parallel_plate_capacitor/config.yaml
  python3 createCOMSOLmesh.py --case parallel_plate_capacitor
  python3 createCOMSOLmesh.py --output-name mesh_comsol.bdf --no-overwrite
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import traceback
from typing import Iterable

import meshio
import numpy as np
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent  # .../MFEMElectrostatics


@dataclass(frozen=True)
class ConversionResult:
    case_name: str
    config_path: Path
    msh_path: Path
    out_path: Path
    wrote: bool


def _load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Config is not a mapping: {path}")
    return data


def _expand_config_path_variants(raw_path: Path, case_dir: Path) -> list[Path]:
    candidates: list[Path] = []
    candidates.append(raw_path)
    if not raw_path.is_absolute():
        candidates.append((case_dir / raw_path).resolve())
        candidates.append((REPO_ROOT / raw_path).resolve())

    raw_str = str(raw_path)
    work_prefix = "/work/"
    if raw_str.startswith(work_prefix):
        rel = raw_str[len(work_prefix):]
        candidates.append((REPO_ROOT / rel).resolve())

    return candidates


def _append_candidate(candidates: list[Path], seen: set[Path], path: Path) -> None:
    p = path.resolve()
    if p not in seen:
        seen.add(p)
        candidates.append(p)


def _resolve_mesh_path(
    raw_mesh_path: str,
    case_dir: Path,
    raw_save_path: str | None,
    prefer_amr_serial: bool,
) -> Path:
    """
    Resolve a mesh path from config.
    Supports:
      - absolute local paths
      - /work/... style paths (mapped to repo root)
      - relative paths (to case dir or repo root)
    """
    raw = Path(raw_mesh_path).expanduser()
    candidates: list[Path] = []
    seen: set[Path] = set()

    mesh_bases = _expand_config_path_variants(raw, case_dir)
    save_bases: list[Path] = []
    if raw_save_path:
        save_path = Path(raw_save_path).expanduser()
        save_bases = _expand_config_path_variants(save_path, case_dir)

    if prefer_amr_serial:
        # Prefer AMR serial mesh first when AMR is enabled in config.
        for save_base in save_bases:
            _append_candidate(candidates, seen, save_base / "amr_mesh" / "amr_mesh_serial.mesh")
            _append_candidate(candidates, seen, save_base / "amr_mesh_serial.mesh")
        for base in mesh_bases:
            _append_candidate(candidates, seen, base / "amr_mesh" / "amr_mesh_serial.mesh")
            _append_candidate(candidates, seen, base / "amr_mesh_serial.mesh")

    for base in mesh_bases:
        _append_candidate(candidates, seen, base)
        if not prefer_amr_serial:
            _append_candidate(candidates, seen, base / "amr_mesh_serial.mesh")
            _append_candidate(candidates, seen, base / "amr_mesh" / "amr_mesh_serial.mesh")

    if not prefer_amr_serial:
        for save_base in save_bases:
            _append_candidate(candidates, seen, save_base / "amr_mesh_serial.mesh")
            _append_candidate(candidates, seen, save_base / "amr_mesh" / "amr_mesh_serial.mesh")

    # Common fallbacks for this validation layout
    if prefer_amr_serial:
        _append_candidate(candidates, seen, case_dir / "mesh" / "amr_mesh" / "amr_mesh_serial.mesh")
        _append_candidate(candidates, seen, case_dir / "mesh" / "amr_mesh_serial.mesh")
    _append_candidate(candidates, seen, case_dir / "mesh" / "mesh.med")
    _append_candidate(candidates, seen, case_dir / "mesh" / "mesh.msh")
    if not prefer_amr_serial:
        _append_candidate(candidates, seen, case_dir / "mesh" / "amr_mesh_serial.mesh")
        _append_candidate(candidates, seen, case_dir / "mesh" / "amr_mesh" / "amr_mesh_serial.mesh")

    for c in candidates:
        if c.exists() and c.is_file():
            return c

    tried = "\n  - ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        f"Could not resolve mesh.path '{raw_mesh_path}' for case '{case_dir.name}'. Tried:\n  - {tried}"
    )


def _find_case_configs(base_dir: Path, selected_cases: set[str] | None) -> list[Path]:
    configs: list[Path] = []
    for d in sorted(base_dir.iterdir()):
        if not d.is_dir():
            continue
        if selected_cases is not None and d.name not in selected_cases:
            continue
        cfg = d / "config.yaml"
        if cfg.exists():
            configs.append(cfg)
    return configs


_MFEM_GEOM_TO_CELL: dict[int, tuple[str, int]] = {
    0: ("vertex", 1),
    1: ("line", 2),
    2: ("triangle", 3),
    3: ("quad", 4),
    4: ("tetra", 4),
    5: ("hexahedron", 8),
    6: ("wedge", 6),
    7: ("pyramid", 5),
}


def _clean_nonempty_lines(path: Path) -> list[str]:
    lines: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if line:
                lines.append(line)
    return lines


def _parse_mfem_cell_section(
    lines: list[str],
    idx: int,
    count: int,
    section_name: str,
    path: Path,
) -> tuple[list[tuple[int, str, list[int]]], int]:
    records: list[tuple[int, str, list[int]]] = []
    for local_i in range(count):
        if idx >= len(lines):
            raise ValueError(
                f"Unexpected EOF while reading '{section_name}' entries in MFEM mesh '{path}' "
                f"(entry {local_i + 1} of {count})."
            )

        fields = lines[idx].split()
        idx += 1
        if len(fields) < 3:
            raise ValueError(
                f"Invalid '{section_name}' row in MFEM mesh '{path}': '{lines[idx - 1]}'. "
                "Expected: '<attribute> <geometry> <vertex ids...>'."
            )

        attr = int(fields[0])
        geom = int(fields[1])
        if geom not in _MFEM_GEOM_TO_CELL:
            known = ", ".join(str(k) for k in sorted(_MFEM_GEOM_TO_CELL))
            raise ValueError(
                f"Unsupported MFEM geometry id {geom} in section '{section_name}' for '{path}'. "
                f"Known ids: {known}."
            )

        cell_type, nnodes = _MFEM_GEOM_TO_CELL[geom]
        conn = [int(v) for v in fields[2:]]
        if len(conn) < nnodes:
            raise ValueError(
                f"Section '{section_name}' in MFEM mesh '{path}' has geometry id {geom} "
                f"({cell_type}) but only {len(conn)} vertex ids, expected {nnodes}."
            )
        records.append((attr, cell_type, conn[:nnodes]))

    return records, idx


def _mfem_records_to_mesh(points: np.ndarray, records: list[tuple[int, str, list[int]]], path: Path) -> meshio.Mesh:
    by_type_conn: dict[str, list[list[int]]] = {}
    by_type_ref: dict[str, list[int]] = {}

    for attr, cell_type, conn in records:
        if cell_type == "vertex":
            # Nastran BDF export in this tool is element-based; standalone point entities
            # (typical only as 1D boundaries) are intentionally skipped.
            continue
        if cell_type not in by_type_conn:
            by_type_conn[cell_type] = []
            by_type_ref[cell_type] = []
        by_type_conn[cell_type].append(conn)
        by_type_ref[cell_type].append(attr)

    cells: list[tuple[str, np.ndarray]] = []
    refs: list[np.ndarray] = []
    npts = int(points.shape[0])
    for cell_type, conn_rows in by_type_conn.items():
        conn = np.asarray(conn_rows, dtype=np.int64)
        if conn.size > 0:
            cmin = int(conn.min())
            cmax = int(conn.max())
            if cmin < 0 or cmax >= npts:
                raise ValueError(
                    f"MFEM mesh '{path}' has invalid vertex ids for '{cell_type}': "
                    f"min={cmin}, max={cmax}, vertices={npts}."
                )
        cells.append((cell_type, conn))
        refs.append(np.asarray(by_type_ref[cell_type], dtype=np.int32))

    if not cells:
        raise ValueError(f"MFEM mesh '{path}' contains no supported element blocks.")

    return meshio.Mesh(points=points, cells=cells, cell_data={"nastran:ref": refs})


def _read_mfem_serial_mesh(path: Path) -> meshio.Mesh:
    lines = _clean_nonempty_lines(path)
    if not lines:
        raise ValueError(f"Empty mesh file: {path}")

    header = lines[0]
    if header.startswith("MFEM NC mesh"):
        raise ValueError(
            f"MFEM NC mesh format is not supported by createCOMSOLmesh.py yet: {path}"
        )
    if not header.startswith("MFEM mesh"):
        raise ValueError(
            f"Not an MFEM serial mesh header in '{path}': '{header}'."
        )

    idx = 1
    points: np.ndarray | None = None
    records: list[tuple[int, str, list[int]]] = []

    while idx < len(lines):
        section = lines[idx].lower()
        idx += 1

        if section == "dimension":
            if idx >= len(lines):
                raise ValueError(f"Missing value after 'dimension' in {path}")
            _ = int(lines[idx])
            idx += 1
            continue

        if section == "elements":
            if idx >= len(lines):
                raise ValueError(f"Missing count after 'elements' in {path}")
            ne = int(lines[idx])
            idx += 1
            parsed, idx = _parse_mfem_cell_section(lines, idx, ne, "elements", path)
            records.extend(parsed)
            continue

        if section == "boundary":
            if idx >= len(lines):
                raise ValueError(f"Missing count after 'boundary' in {path}")
            nb = int(lines[idx])
            idx += 1
            parsed, idx = _parse_mfem_cell_section(lines, idx, nb, "boundary", path)
            records.extend(parsed)
            continue

        if section == "vertices":
            if idx >= len(lines):
                raise ValueError(f"Missing count after 'vertices' in {path}")
            nv = int(lines[idx])
            idx += 1
            if idx >= len(lines):
                raise ValueError(f"Missing vertex dimension after 'vertices' in {path}")
            if lines[idx].lower() == "nodes":
                raise ValueError(
                    f"MFEM mesh '{path}' uses a 'nodes' section for coordinates, which is not supported yet."
                )
            vdim = int(lines[idx])
            idx += 1
            pts: list[list[float]] = []
            for vid in range(nv):
                if idx >= len(lines):
                    raise ValueError(
                        f"Unexpected EOF while reading vertices in '{path}' "
                        f"(vertex {vid + 1} of {nv})."
                    )
                vals = lines[idx].split()
                idx += 1
                if len(vals) < vdim:
                    raise ValueError(
                        f"Vertex row has {len(vals)} values, expected {vdim}, in '{path}'."
                    )
                pts.append([float(vals[j]) for j in range(vdim)])
            points = np.asarray(pts, dtype=np.float64)
            continue

        if section == "nodes":
            raise ValueError(
                f"MFEM mesh '{path}' includes a 'nodes' section, which is not supported yet."
            )

        if section in {"attribute_sets", "bdr_attribute_sets", "vertex_parents", "coarse_elements"}:
            if idx >= len(lines):
                raise ValueError(f"Missing count after '{section}' in {path}")
            nrows = int(lines[idx])
            idx += 1 + nrows
            continue

        # Unknown/optional sections are ignored by default.

    if points is None:
        raise ValueError(f"MFEM mesh '{path}' does not contain a 'vertices' section.")

    return _mfem_records_to_mesh(points, records, path)


def _first_nonempty_noncomment_line(path: Path) -> str:
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if line:
                return line
    return ""


def _looks_like_mfem_serial_mesh(path: Path) -> bool:
    head = _first_nonempty_noncomment_line(path)
    return head.startswith("MFEM mesh") or head.startswith("MFEM NC mesh")


def _load_input_mesh(path: Path) -> meshio.Mesh:
    if _looks_like_mfem_serial_mesh(path):
        return _read_mfem_serial_mesh(path)

    try:
        return meshio.read(str(path))
    except Exception as meshio_err:
        if path.suffix.lower() == ".mesh":
            try:
                return _read_mfem_serial_mesh(path)
            except Exception as mfem_err:
                raise RuntimeError(
                    f"Failed to load '{path}' as meshio format and as MFEM serial mesh.\n"
                    f"meshio error: {meshio_err}\n"
                    f"MFEM parser error: {mfem_err}"
                ) from mfem_err
        raise


def _write_nastran_native(
    msh_path: Path,
    out_path: Path,
    point_format: str,
    cell_format: str,
    line_element: str,
) -> None:
    mesh = _load_input_mesh(msh_path)

    # meshio Nastran writer uses cell_data["nastran:ref"] as the element property field (PID).
    # Map gmsh physical IDs into nastran:ref so COMSOL can partition imported domains/boundaries.
    if "nastran:ref" not in mesh.cell_data and "gmsh:physical" in mesh.cell_data:
        mesh.cell_data["nastran:ref"] = [
            np.asarray(block, dtype=np.int32) for block in mesh.cell_data["gmsh:physical"]
        ]

    # meshio's fixed-small output is hard to read and its fixed-large/free formatting can
    # assert on some coordinates in certain meshio versions. Use a robust free-format writer
    # for deterministic, comma-separated output.
    _write_nastran_free(mesh, out_path, line_element=line_element)


def _fmt_float_free(x: float) -> str:
    # Keep numbers compact but explicit. Avoid "-0".
    if abs(x) < 1.0e-15:
        return "0.0"
    return f"{x:.15g}"


def _cell_refs(mesh: meshio.Mesh) -> list[np.ndarray]:
    for key in ("nastran:ref", "gmsh:physical", "medit:ref", "cell_tags"):
        if key in mesh.cell_data:
            refs = [np.asarray(block, dtype=np.int64) for block in mesh.cell_data[key]]
            if len(refs) == len(mesh.cells):
                return refs
    return [np.ones(len(block.data), dtype=np.int64) for block in mesh.cells]


def _card_and_node_count(cell_type: str, line_element: str) -> tuple[str, int]:
    line_card = "CROD" if line_element == "crod" else "CBAR"
    mapping = {
        "line": (line_card, 2),
        "triangle": ("CTRIA3", 3),
        "quad": ("CQUAD4", 4),
        "tetra": ("CTETRA", 4),
        "hexahedron": ("CHEXA", 8),
        "wedge": ("CPENTA", 6),
        "pyramid": ("CPYRA", 5),
    }
    if cell_type not in mapping:
        supported = ", ".join(sorted(mapping.keys()))
        raise ValueError(f"Unsupported mesh cell type '{cell_type}' for BDF export. Supported: {supported}")
    return mapping[cell_type]


def _write_nastran_free(mesh: meshio.Mesh, out_path: Path, line_element: str) -> None:
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError(f"Invalid points array shape for BDF export: {points.shape}")
    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(points.shape[0], dtype=np.float64)])
    elif points.shape[1] > 3:
        points = points[:, :3]

    refs_by_block = _cell_refs(mesh)

    with out_path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("$ Nastran file written by createCOMSOLmesh.py (free format)\n")
        fh.write("BEGIN BULK\n")

        for nid, xyz in enumerate(points, start=1):
            x, y, z = (float(v) for v in xyz)
            fh.write(f"GRID,{nid},,{_fmt_float_free(x)},{_fmt_float_free(y)},{_fmt_float_free(z)}\n")

        eid = 1
        for block_idx, block in enumerate(mesh.cells):
            card, nnodes = _card_and_node_count(block.type, line_element=line_element)
            conn = np.asarray(block.data, dtype=np.int64)
            if conn.ndim != 2 or conn.shape[1] < nnodes:
                raise ValueError(
                    f"Cell block '{block.type}' has connectivity shape {conn.shape}; expected at least {nnodes} nodes."
                )
            pids = refs_by_block[block_idx]
            if len(pids) != len(conn):
                raise ValueError(
                    f"Cell data length mismatch for '{block.type}': {len(pids)} refs vs {len(conn)} elements."
                )

            for local_idx, nodes in enumerate(conn):
                pid = int(pids[local_idx])
                node_ids = (nodes[:nnodes] + 1).tolist()
                node_fields = ",".join(str(int(n)) for n in node_ids)
                fh.write(f"{card},{eid},{pid},{node_fields}\n")
                eid += 1

        fh.write("ENDDATA\n")


def _convert_case_config(
    config_path: Path,
    output_name: str,
    overwrite: bool,
    point_format: str,
    cell_format: str,
    line_element: str,
) -> ConversionResult:
    cfg = _load_yaml(config_path)
    case_dir = config_path.parent
    case_name = case_dir.name

    mesh_cfg = cfg.get("mesh", {})
    if not isinstance(mesh_cfg, dict):
        raise ValueError(f"Missing/invalid 'mesh' section in {config_path}")
    raw_mesh_path = mesh_cfg.get("path")
    if not isinstance(raw_mesh_path, str) or not raw_mesh_path:
        raise ValueError(f"Missing mesh.path in {config_path}")

    raw_save_path = cfg.get("save_path")
    if raw_save_path is not None and not isinstance(raw_save_path, str):
        raise ValueError(f"Invalid save_path in {config_path}; expected string if provided.")

    amr_cfg = mesh_cfg.get("AMR", {})
    if not isinstance(amr_cfg, dict):
        raise ValueError(f"Invalid mesh.AMR section in {config_path}; expected mapping if provided.")
    prefer_amr_serial = bool(amr_cfg.get("enable", False))

    msh_path = _resolve_mesh_path(
        raw_mesh_path,
        case_dir,
        raw_save_path,
        prefer_amr_serial=prefer_amr_serial,
    )
    out_path = msh_path.parent / output_name
    if out_path.exists() and not overwrite:
        return ConversionResult(case_name, config_path, msh_path, out_path, wrote=False)

    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    _write_nastran_native(
        msh_path=msh_path,
        out_path=tmp_path,
        point_format=point_format,
        cell_format=cell_format,
        line_element=line_element,
    )
    tmp_path.replace(out_path)
    return ConversionResult(case_name, config_path, msh_path, out_path, wrote=True)


def _iter_cases(
    configs: Iterable[Path],
    output_name: str,
    overwrite: bool,
    show_traceback: bool,
    point_format: str,
    cell_format: str,
    line_element: str,
) -> tuple[int, int]:
    n_ok = 0
    n_err = 0
    for cfg in configs:
        try:
            result = _convert_case_config(
                config_path=cfg,
                output_name=output_name,
                overwrite=overwrite,
                point_format=point_format,
                cell_format=cell_format,
                line_element=line_element,
            )
            status = "wrote" if result.wrote else "exists"
            print(f"[OK] {result.case_name}: {result.msh_path} -> {result.out_path} ({status})")
            n_ok += 1
        except Exception as e:
            msg = str(e).strip()
            if not msg:
                msg = repr(e)
            print(f"[ERR] {cfg.parent.name}: {type(e).__name__}: {msg}")
            if show_traceback:
                print(traceback.format_exc().rstrip())
            n_err += 1
    return n_ok, n_err


def _resolve_config_arg(path: Path, cwd: Path) -> Path:
    p = path.expanduser()
    if not p.is_absolute():
        p = (cwd / p).resolve()
    else:
        p = p.resolve()
    return p


def _collect_configs(base_dir: Path, selected_cases: set[str] | None, config_args: list[Path], cwd: Path) -> list[Path]:
    if config_args:
        cfgs: list[Path] = []
        for c in config_args:
            cp = _resolve_config_arg(c, cwd)
            if cp.is_dir():
                cp = cp / "config.yaml"
            if not cp.exists():
                raise FileNotFoundError(f"Config path does not exist: {cp}")
            if not cp.is_file():
                raise ValueError(f"Config path is not a file: {cp}")
            cfgs.append(cp)
        return sorted(set(cfgs))
    return _find_case_configs(base_dir, selected_cases)


def main() -> int:
    ap = argparse.ArgumentParser(description="Create COMSOL-importable BDF meshes from COMSOLValidation configs.")
    ap.add_argument(
        "--base-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory containing COMSOLValidation case subdirectories.",
    )
    ap.add_argument(
        "-c",
        "--config",
        action="append",
        type=Path,
        default=[],
        help="Config path (or case directory). Repeat for multiple. If set, --base-dir/--case are ignored.",
    )
    ap.add_argument(
        "--case",
        action="append",
        default=[],
        help="Case folder name to process. Repeat for multiple cases. Default: all cases with config.yaml.",
    )
    ap.add_argument(
        "--output-name",
        default="mesh.bdf",
        help="Output file name written next to each input mesh.",
    )
    ap.add_argument(
        "--point-format",
        choices=("fixed-small", "fixed-large", "free"),
        default="free",
        help="Requested NASTRAN point format (default: free).",
    )
    ap.add_argument(
        "--cell-format",
        choices=("fixed-small", "fixed-large", "free"),
        default="free",
        help="Requested NASTRAN cell format (default: free).",
    )
    ap.add_argument(
        "--line-element",
        choices=("crod", "cbar"),
        default="crod",
        help="NASTRAN card for 2-node boundary line elements (default: crod).",
    )
    ow = ap.add_mutually_exclusive_group()
    ow.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        help="Overwrite existing outputs (default behavior).",
    )
    ow.add_argument(
        "--no-overwrite",
        dest="overwrite",
        action="store_false",
        help="Skip cases where the output already exists.",
    )
    ap.add_argument(
        "--traceback",
        action="store_true",
        help="Print full traceback for conversion errors.",
    )
    ap.set_defaults(overwrite=True)
    args = ap.parse_args()

    if args.point_format != "free" or args.cell_format != "free":
        print(
            "[INFO] Using robust free-format BDF writer; "
            f"--point-format={args.point_format} and --cell-format={args.cell_format} are ignored."
        )

    base_dir = args.base_dir.resolve()
    if not base_dir.exists() or not base_dir.is_dir():
        raise NotADirectoryError(base_dir)

    selected_cases = set(args.case) if args.case else None
    configs = _collect_configs(
        base_dir=base_dir,
        selected_cases=selected_cases,
        config_args=args.config,
        cwd=Path.cwd(),
    )
    if not configs:
        print(f"[INFO] No case configs found under {base_dir}")
        return 0

    n_ok, n_err = _iter_cases(
        configs=configs,
        output_name=args.output_name,
        overwrite=args.overwrite,
        show_traceback=args.traceback,
        point_format=args.point_format,
        cell_format=args.cell_format,
        line_element=args.line_element,
    )
    print(f"[DONE] Converted {n_ok}/{len(configs)} cases (errors: {n_err}).")
    return 1 if n_err > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
