#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Set
import pprint

try:
    import yaml  # type: ignore
except Exception as ex:
    raise SystemExit(
        "Missing dependency: PyYAML. Install with: pip install pyyaml\n"
        f"Import error: {ex!r}"
    )


# ---------------------- IO helpers ----------------------

def read_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def write_text(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

def resolve_path(base_dir: str, rel_or_abs: str) -> str:
    if os.path.isabs(rel_or_abs):
        return rel_or_abs
    return os.path.normpath(os.path.join(base_dir, rel_or_abs))

def py_repr(obj: Any) -> str:
    return pprint.pformat(
        obj,
        width=88,        # keeps lines reasonable
        compact=False,   # force multiline for nested dicts
        sort_dicts=True  # deterministic output
    )

def angle_from_filename(fn: str) -> Optional[float]:
    base = os.path.splitext(os.path.basename(fn))[0]
    if not base.isdigit():
        return None
    return float(int(base))


# ---------------------- Geometry extraction ----------------------

def pick_chain_from_instances(instances: Any) -> Optional[List[Any]]:
    if not isinstance(instances, list):
        return None
    for inst in instances:
        if not isinstance(inst, dict):
            continue
        pts_lists = inst.get("pts_lists")
        if not isinstance(pts_lists, list):
            continue
        for chain in pts_lists:
            if isinstance(chain, list) and len(chain) > 0:
                return chain
    return None

def pick_chain_from_json_payload(d: Dict[str, Any]) -> Optional[List[Any]]:
    if "instances" in d:
        c = pick_chain_from_instances(d.get("instances"))
        if c:
            return c

    pts_lists = d.get("pts_lists")
    if isinstance(pts_lists, list):
        for chain in pts_lists:
            if isinstance(chain, list) and chain:
                return chain

    pts = d.get("pts")
    if isinstance(pts, list) and pts:
        return pts

    outer = d.get("outer")
    if isinstance(outer, list) and outer:
        return outer

    return None


# ---------------------- Indexing on disk ----------------------

@dataclass
class Entry:
    bc_folder: str
    subcomponent: str
    angle: float
    json_path: str

def infer_group_from_bc_folder(bc: str) -> str:
    return "ptfe" if bc.strip() == "NA" else "electrode"

def index_pts_root(pts_root: str) -> List[Entry]:
    out: List[Entry] = []
    if not os.path.isdir(pts_root):
        return out

    for bc in sorted(os.listdir(pts_root)):
        bc_path = os.path.join(pts_root, bc)
        if not os.path.isdir(bc_path):
            continue

        for sub in sorted(os.listdir(bc_path)):
            sub_path = os.path.join(bc_path, sub)
            if not os.path.isdir(sub_path):
                continue

            for fn in sorted(os.listdir(sub_path)):
                if not fn.lower().endswith(".json"):
                    continue
                ang = angle_from_filename(fn)
                if ang is None:
                    continue
                fp = os.path.join(sub_path, fn)
                out.append(Entry(bc_folder=bc, subcomponent=sub, angle=ang, json_path=fp))

    return out


# ---------------------- Angle overwrites ----------------------

def normalize_angle_overwrites(cfg: Dict[str, Any]) -> Dict[str, Dict[float, float]]:
    ao = cfg.get("angle_overwrites") or {}
    out: Dict[str, Dict[float, float]] = {}
    if not isinstance(ao, dict):
        return out
    for sub, mapping in ao.items():
        sub_s = str(sub).strip()
        if not sub_s or not isinstance(mapping, dict):
            continue
        m2: Dict[float, float] = {}
        for k, v in mapping.items():
            try:
                kk = float(k)
                vv = float(v)
            except Exception:
                continue
            m2[kk] = vv
        if m2:
            out[sub_s] = m2
    return out


# ---------------------- External components ----------------------

def apply_offset_instr(pts: List[Any], ox: float, oy: float) -> List[Any]:
    out: List[Any] = []
    for ins in pts or []:
        if not isinstance(ins, list) or not ins:
            out.append(ins)
            continue

        if isinstance(ins[0], str):
            t = ins[0]
            try:
                if t == "line" and len(ins) >= 3:
                    out.append(["line", float(ins[1]) + ox, float(ins[2]) + oy])
                elif t == "arc" and len(ins) >= 5:
                    x, y, cx, cy = float(ins[1]), float(ins[2]), float(ins[3]), float(ins[4])
                    if len(ins) >= 6:
                        out.append(["arc", x + ox, y + oy, cx + ox, cy + oy, bool(ins[5])])
                    else:
                        out.append(["arc", x + ox, y + oy, cx + ox, cy + oy])
                else:
                    out.append(ins)
            except Exception:
                out.append(ins)
        else:
            if len(ins) >= 2:
                try:
                    out.append([float(ins[0]) + ox, float(ins[1]) + oy] + ins[2:])
                except Exception:
                    out.append(ins)
            else:
                out.append(ins)
    return out
def build_group_lookup(cfg: Dict[str, Any]) -> Dict[str, str]:
    """
    cfg.groups:
      Electrode: [Anode, Bell, ...]
      PTFE:      [NA]
      xenon:     [LXe0, GXe0]
    Returns mapping of component-name -> group-key (lowercase).
    """
    groups = cfg.get("groups") or {}
    if not isinstance(groups, dict):
        return {}

    out: Dict[str, str] = {}
    for gname, items in groups.items():
        if not isinstance(items, list):
            continue
        gkey = str(gname).strip().lower()
        for it in items:
            name = str(it).strip()
            if name:
                out[name] = gkey
    return out

def load_external_component_per_angle(project_root: str, cfg_dir: str, rec: Dict[str, Any], angles: List[float]) -> Tuple[str, str, Dict[float, List[Any]]]:
    name = str(rec.get("name") or "").strip()
    jpath = str(rec.get("json_path") or "").strip()
    if not name or not jpath:
        raise ValueError(f"Bad extra_components entry: {rec!r}")

    group = str(rec.get("group") or "").strip().lower()
    if not group:
        if name in ("LXe", "GXe"):
            group = "xenon"
        elif name == "PMT":
            group = "electrode"
        else:
            group = "ptfe"

    ox = float(rec.get("offset_x") or 0.0)
    oy = float(rec.get("offset_y") or 0.0)

    # IMPORTANT: resolve json_path relative to project root (sibling to pts_overwrites/out)
    full = resolve_path(project_root, jpath)
    payload = read_json(full)

    def pts_from_obj(obj: Any) -> List[Any]:
        if not isinstance(obj, dict):
            return []
        if isinstance(obj.get("pts"), list):
            return obj["pts"]
        if isinstance(obj.get("outer"), list):
            return obj["outer"]
        return []

    per: Dict[float, List[Any]] = {}

    if isinstance(payload, dict) and "by_angle_deg" in payload:
        by = payload.get("by_angle_deg") or {}
        defaults = payload.get("defaults") or {}
        dpts = apply_offset_instr(pts_from_obj(defaults), ox, oy)

        for a in angles:
            key = str(int(round(float(a))))
            if isinstance(by, dict) and key in by:
                per[a] = apply_offset_instr(pts_from_obj(by[key]), ox, oy)
            else:
                per[a] = dpts
    else:
        pts0: List[Any] = []
        if isinstance(payload, dict):
            pts0 = pts_from_obj(payload)
        elif isinstance(payload, list):
            pts0 = payload
        pts0 = apply_offset_instr(pts0, ox, oy)
        for a in angles:
            per[a] = pts0

    return name, group, per


# ---------------------- Codegen per angle ----------------------

def build_module_text_for_angle(
    cfg_path: str,
    angle: float,
    ptfe_sketches: Dict[str, Dict[str, Any]],
    electrode_sketches: Dict[str, Dict[str, Any]],
    xenon_sketches: Dict[str, Dict[str, Any]],
) -> str:
    lines: List[str] = []
    lines.append("# Auto-generated. Do not edit by hand.")
    lines.append(f"# Config: {cfg_path}")
    lines.append(f"# Slice angle: {angle:.1f} deg")
    lines.append("")
    lines.append("import math")
    lines.append("")
    lines.append(f"SLICE_ANGLE_DEG = {angle:.1f}")
    lines.append("")
    lines.append(f"_PTFE = {py_repr(ptfe_sketches)}")
    lines.append(f"_ELECTRODE = {py_repr(electrode_sketches)}")
    lines.append(f"_XENON = {py_repr(xenon_sketches)}")
    lines.append("")
    lines.append("def build_sketch_dicts(shrinkage_factor: float = 1.0):")
    lines.append("    ptfe_sketches = dict(_PTFE)")
    lines.append("    electrode_sketches = dict(_ELECTRODE)")
    lines.append("    electrode_sketches['shrinkage_factor'] = shrinkage_factor")
    lines.append("    xenon_sketches = dict(_XENON)")
    lines.append("    manual_mapping = {}")
    lines.append("    return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping")
    lines.append("")
    return "\n".join(lines)


# ---------------------- Main ----------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg", default="./pts_overwrites/XENT-TPC_20250428_overwrites.yaml")
    ap.add_argument("--out-dir", default=None, help="Defaults to <project_root>/out/slices_py")
    args = ap.parse_args(argv)

    cfg_path = os.path.abspath(args.cfg)
    cfg_dir = os.path.dirname(cfg_path)

    # project root is parent of pts_overwrites (and sibling to out)
    project_root = os.path.abspath(os.path.join(cfg_dir, os.pardir))

    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise SystemExit("Config YAML must be a mapping/object.")

    inputs = cfg.get("inputs") or {}

    group_lookup = build_group_lookup(cfg)
    # IMPORTANT: resolve pts_root relative to project root
    pts_root = resolve_path(project_root, str(inputs.get("pts_root") or "./out/pts_out"))

    out_dir = args.out_dir
    if out_dir is None:
        out_dir = os.path.join(project_root, "out", "slices_py")
    else:
        out_dir = resolve_path(project_root, out_dir)

    entries = index_pts_root(pts_root)
    if not entries:
        raise SystemExit(f"No slice JSON files discovered under: {pts_root}")

    angles = sorted({e.angle for e in entries})

    lookup: Dict[Tuple[str, str, float], str] = {}
    for e in entries:
        lookup[(e.bc_folder, e.subcomponent, e.angle)] = e.json_path

    ao = normalize_angle_overwrites(cfg)

    external_per_angle: List[Tuple[str, str, Dict[float, List[Any]]]] = []
    for rec in cfg.get("extra_components") or []:
        if not isinstance(rec, dict):
            continue
        external_per_angle.append(load_external_component_per_angle(project_root, cfg_dir, rec, angles))

    os.makedirs(os.path.abspath(out_dir), exist_ok=True)

    bc_sub_pairs = sorted({(e.bc_folder, e.subcomponent) for e in entries})

    for a in angles:
        ptfe_sketches: Dict[str, Dict[str, Any]] = {}
        electrode_sketches: Dict[str, Dict[str, Any]] = {}
        xenon_sketches: Dict[str, Dict[str, Any]] = {}

        for bc, sub in bc_sub_pairs:
            src_angle = ao.get(sub, {}).get(a, a)
            path = lookup.get((bc, sub, src_angle))
            if not path:
                continue

            d = read_json(path)
            if not isinstance(d, dict):
                continue

            pts = pick_chain_from_json_payload(d)
            if not pts:
                continue

            group = group_lookup.get(bc)  # bc is the top-level folder: Anode, NA, ...
            if group is None:
                # strict is better for debugging; you can change to a default if you want
                raise KeyError(f"BC folder {bc!r} not assigned in cfg.groups")

            if group == "ptfe":
                ptfe_sketches[sub] = {"pts": pts}
            elif group == "xenon":
                xenon_sketches[sub] = {"pts": pts}
            else:  # electrode or any other future group you treat as electrode
                electrode_sketches[sub] = {"pts": pts}


        for name, _unused_group, per in external_per_angle:
            pts = per.get(a) or []
            if not pts:
                continue

            group = group_lookup.get(name)
            if group is None:
                raise KeyError(f"Extra component {name!r} not assigned in cfg.groups")

            if group == "xenon":
                xenon_sketches[name] = {"pts": pts}
            elif group == "ptfe":
                ptfe_sketches[name] = {"pts": pts}
            else:
                electrode_sketches[name] = {"pts": pts}


        mod = build_module_text_for_angle(cfg_path, a, ptfe_sketches, electrode_sketches, xenon_sketches)
        out_path = os.path.join(os.path.abspath(out_dir), f"slice_{int(a):03d}.py")
        write_text(out_path, mod)

    print(f"Wrote {len(angles)} slice modules to: {os.path.abspath(out_dir)}")
    print(f"Angles: {angles}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
