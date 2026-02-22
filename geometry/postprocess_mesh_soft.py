#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import yaml

"""
Soft autogen script.

Hard rules:
- NO voltages
- NO fieldcage network
- NO epsilon_r auto-assignment
- Dirichlet boundaries NEVER get a 'value'
- All paths are resolved relative to the directory of config.yaml
"""

# ============================================================
# CLI
# ============================================================
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-c",
        "--config",
        required=True,
        type=Path,
        help="Path to config.yaml",
    )
    return ap.parse_args()


# ============================================================
# YAML helpers
# ============================================================
def yaml_load(path: Path) -> dict:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def yaml_dump(path: Path, data: dict) -> None:
    path.write_text(
        yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )


# ============================================================
# phys_map helpers
# ============================================================
def extract_physical_names(msh_path: Path, out_txt: Path) -> None:
    in_block = False
    with msh_path.open("r", encoding="utf-8", errors="replace") as f_in, out_txt.open(
        "w", encoding="utf-8"
    ) as f_out:
        for line in f_in:
            s = line.strip()
            if s == "$PhysicalNames":
                in_block = True
                continue
            if s == "$EndPhysicalNames":
                break
            if in_block:
                f_out.write(line)


def parse_phys_map(path: Path):
    with path.open(encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]

    n_total = int(lines[0])
    entries = lines[1:]

    regions = {}
    boundaries = {}

    for line in entries:
        dim_s, tag_s, name_s = line.split(maxsplit=2)
        dim = int(dim_s)
        tag = int(tag_s)
        name = name_s.strip('"')

        if dim == 1:
            boundaries[name] = {"bdr_id": tag}
        elif dim == 2:
            regions[name] = {"attr_id": tag}
        else:
            raise ValueError(f"Unexpected dim {dim}")

    print(f"[parse] phys_map entries={n_total}, regions={len(regions)}, boundaries={len(boundaries)}")
    return regions, boundaries


# ============================================================
# autoappend helpers
# ============================================================
def mesh_autoappend_enabled(cfg: dict) -> bool:
    mesh = cfg.get("mesh", {})
    return isinstance(mesh, dict) and bool(mesh.get("autoappend", False))


def resolve_and_update_mesh_path(cfg: dict, msh_path: Path) -> None:
    mesh_cfg = cfg.get("mesh")
    if not isinstance(mesh_cfg, dict):
        return
    if msh_path.exists():
        mesh_cfg["path"] = str(msh_path.resolve())


def apply_autogen_sections_inplace(cfg: dict, autogen: dict) -> dict:
    out = dict(cfg)
    out["materials"] = autogen["materials"]
    out["boundaries"] = autogen["boundaries"]
    return out


# ============================================================
# main
# ============================================================
def main():
    args = parse_args()

    config_path = args.config.resolve()
    base_dir = config_path.parent

    # All paths relative to config directory
    mesh_dir = base_dir / "mesh"
    med_path = mesh_dir / "mesh.med"
    msh_path = mesh_dir / "mesh22.msh"
    phys_map_path = mesh_dir / "phys_map.txt"
    output_path = mesh_dir / "config_autogen.yaml"

    cfg = yaml_load(config_path)

    # ------------------------------------------------------------
    # 0) Mesh preprocessing
    # ------------------------------------------------------------
    mesh_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["gmsh", "-0", str(med_path), "-format", "msh2", "-o", str(msh_path)],
        check=True,
    )

    extract_physical_names(msh_path, phys_map_path)

    # ------------------------------------------------------------
    # 1) Parse phys_map
    # ------------------------------------------------------------
    regions_raw, boundaries_raw = parse_phys_map(phys_map_path)

    # ------------------------------------------------------------
    # 2) Materials (attr_id only)
    # ------------------------------------------------------------
    materials = {
        name: {"attr_id": meta["attr_id"]}
        for name, meta in regions_raw.items()
    }

    # ------------------------------------------------------------
    # 3) Boundaries
    # ------------------------------------------------------------
    existing_boundaries = cfg.get("boundaries", {})
    if not isinstance(existing_boundaries, dict):
        existing_boundaries = {}

    boundaries_out = {}

    for name, meta in boundaries_raw.items():
        bdr_id = meta["bdr_id"]

        # default: dirichlet, no value
        out = {"bdr_id": bdr_id, "type": "dirichlet"}

        # preserve existing non-dirichlet BCs
        ex = existing_boundaries.get(name)
        if isinstance(ex, dict) and ex.get("type") != "dirichlet":
            out = dict(ex)
            out["bdr_id"] = bdr_id

        boundaries_out[name] = out

    blank = [n for n, v in boundaries_out.items() if v["type"] == "dirichlet"]
    print(f"[info] Dirichlet boundaries with blank values: {len(blank)}")

    # ------------------------------------------------------------
    # 4) Write config_autogen.yaml
    # ------------------------------------------------------------
    with output_path.open("w", encoding="utf-8") as out:
        out.write("materials:\n")
        for m in sorted(materials):
            out.write(f"  {m}:\n")
            out.write(f"    attr_id: {materials[m]['attr_id']}\n")

        out.write("\nboundaries:\n")
        for name in sorted(boundaries_out):
            v = boundaries_out[name]
            out.write(f"  {name}:\n")
            out.write(f"    bdr_id: {v['bdr_id']}\n")
            out.write(f"    type: {v['type']}\n")
            for k in sorted(v):
                if k in ("bdr_id", "type"):
                    continue
                out.write(f"    {k}: {v[k]}\n")

    print(f"[done] Wrote {output_path}")

    # ------------------------------------------------------------
    # 5) Optional autoappend
    # ------------------------------------------------------------
    if mesh_autoappend_enabled(cfg):
        resolve_and_update_mesh_path(cfg, msh_path)
        new_cfg = apply_autogen_sections_inplace(
            cfg,
            {"materials": materials, "boundaries": boundaries_out},
        )
        yaml_dump(config_path, new_cfg)
        print(f"[autoappend] Updated {config_path}")


if __name__ == "__main__":
    main()