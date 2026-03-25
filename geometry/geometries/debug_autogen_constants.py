import glob
import copy
import json
import os
import yaml


EXCLUDED_PREFIXES = ()#("XENTTPCC", "XENTTPCB", "XENTTPCD", "XENTTPCE", "XENTTPCF")
ALWAYS_EXCLUDED_COMPONENTS = {
    "XENTTPCFPart5521Part5521_0",
    "XENTTPCFPart5521Part5521_1",
}
ALWAYS_EXCLUDED_PREFIXES = ("XENTTPCG", "XENTTPCH")
DEFAULT_ANGLE_DIR = "/work/geometry/createGeometryFromCAD/DXF_slices_parts/slice_022.50deg"
DEFAULT_LAYOUT_FILE = "cleanup_autogen.yaml"
TARGET_DEBUG_COMPONENT = "XENTTPCAWarmPart088Part088"


def _is_excluded_component(key: str) -> bool:
    if key.startswith(EXCLUDED_PREFIXES):
        return True
    if key.startswith(ALWAYS_EXCLUDED_PREFIXES):
        return True
    if key in ALWAYS_EXCLUDED_COMPONENTS:
        return True
    # Defensive: exclude derivative names if cleanup adds suffixes.
    for base in ALWAYS_EXCLUDED_COMPONENTS:
        if key.startswith(base + "_"):
            return True
    return False

def _is_drawable_component(data: dict) -> bool:
    if not isinstance(data, dict):
        return False
    if bool(data.get("hull", False)):
        return True
    if "pts" in data and isinstance(data["pts"], list) and len(data["pts"]) > 0:
        return True
    if "Radius" in data:
        return True
    if "Height" in data and "Width" in data:
        return True
    return False


def _safe_min_y(data: dict):
    pts = data.get("pts", [])
    ys = []
    for item in pts:
        if isinstance(item, list) and len(item) >= 3 and isinstance(item[2], (int, float)):
            ys.append(float(item[2]))
    if not ys:
        return None
    return min(ys)


def _load_layout_manifest(angle_dir: str):
    layout_path = os.environ.get("DEBUG_AUTOGEN_LAYOUT_FILE", DEFAULT_LAYOUT_FILE)
    if not os.path.isabs(layout_path):
        layout_path = os.path.join(angle_dir, layout_path)
    if not os.path.isfile(layout_path):
        return None, layout_path
    with open(layout_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Invalid cleanup autogen layout at {layout_path}")
    return data, layout_path


def _clone_named_components(merged, names):
    missing = [name for name in names if name not in merged]
    if missing:
        raise KeyError("cleanup_autogen references missing components: %s" % ", ".join(sorted(missing)))
    return {name: copy.deepcopy(merged[name]) for name in names}


def _bucketize_from_layout(merged, layout):
    materials = layout.get("materials") or {}
    electrodes_cfg = layout.get("electrodes") or {}

    used_names = set()
    xenon_sketches = {}
    electrode_sketches = {}

    supported_materials = {"GXe", "LXe", "PTFE"}
    for material_name, entry in materials.items():
        if material_name not in supported_materials:
            raise ValueError(
                "Unsupported material in cleanup_autogen: %s. Supported: GXe, LXe, PTFE" % material_name
            )
        if not isinstance(entry, dict):
            raise ValueError(f"Invalid material entry for {material_name}: {entry!r}")
        components = list(entry.get("components", []) or [])
        if entry.get("default"):
            continue
        if material_name == "PTFE":
            continue
        if not components:
            continue
        if len(components) == 1 and components[0] == material_name:
            xenon_sketches[material_name] = copy.deepcopy(merged[components[0]])
        else:
            xenon_sketches[material_name] = {
                "hull": True,
                "sub_sketches": _clone_named_components(merged, components),
            }
        used_names.update(components)

    for electrode_name, entry in electrodes_cfg.items():
        if not isinstance(entry, dict):
            raise ValueError(f"Invalid electrode entry for {electrode_name}: {entry!r}")
        components = list(entry.get("components", []) or [])
        if not components:
            continue
        electrode_sketches[electrode_name] = {
            "hull": True,
            "sub_sketches": _clone_named_components(merged, components),
        }
        used_names.update(components)

    ptfe_sketches = {}
    for name, data in merged.items():
        if name in used_names:
            continue
        ptfe_sketches[name] = copy.deepcopy(data)

    return ptfe_sketches, electrode_sketches, xenon_sketches


def _audit_part088_gaps(merged: dict):
    keys = sorted(k for k in merged.keys() if k.startswith("XENTTPCAWarmPart088Part088_"))
    if not keys:
        return

    # numeric suffix coverage
    idxs = []
    for k in keys:
        try:
            idxs.append(int(k.rsplit("_", 1)[1]))
        except Exception:
            continue
    idxs = sorted(set(idxs))
    if idxs:
        miss = [i for i in range(idxs[0], idxs[-1] + 1) if i not in idxs]
        if miss:
            print("[debug_autogen][Part088] missing index labels: %s" % ", ".join(str(i) for i in miss[:30]))
            if len(miss) > 30:
                print("[debug_autogen][Part088] ... +%d more missing labels" % (len(miss) - 30))
        else:
            print("[debug_autogen][Part088] index labels are contiguous (%d..%d)" % (idxs[0], idxs[-1]))

    # geometric spacing audit by min y
    rows = []
    for k in keys:
        y = _safe_min_y(merged[k])
        if y is None:
            continue
        rows.append((y, k))
    rows.sort()
    if len(rows) < 3:
        return

    diffs = [rows[i + 1][0] - rows[i][0] for i in range(len(rows) - 1)]
    pos = sorted(d for d in diffs if d > 0)
    if not pos:
        return
    med = pos[len(pos) // 2]
    if med <= 0:
        return

    flagged = []
    threshold = 1.5 * med
    for i, d in enumerate(diffs):
        if d > threshold:
            y0, k0 = rows[i]
            y1, k1 = rows[i + 1]
            flagged.append((k0, k1, y0, y1, d))

    if flagged:
        print(
            "[debug_autogen][Part088] geometric spacing gaps > %.3f (median step %.3f): %d"
            % (threshold, med, len(flagged))
        )
        for k0, k1, y0, y1, d in flagged[:20]:
            print(
                "  [gap] %s (min_y=%.3f) -> %s (min_y=%.3f), delta=%.3f"
                % (k0, y0, k1, y1, d)
            )
        if len(flagged) > 20:
            print("[debug_autogen][Part088] ... +%d more geometric gaps" % (len(flagged) - 20))

def build_sketch_dicts(_unused=None):
    """
    Debug-only helper.

    Returns:
      ptfe_sketches       : dict
      electrode_sketches  : dict
      xenon_sketches      : dict
      manual_mapping      : dict (always empty)

    Each sketches dict has the form:
      { comp_key : { "pts": [...] }, ... }

    In debug mode we keep assignment stable and avoid arbitrary
    multi-bucket distribution that can alter partition behavior.
    """

    merged = {}
    angle_dir = os.environ.get("DEBUG_AUTOGEN_JSON_DIR", DEFAULT_ANGLE_DIR)
    components = ["*.json"]
    excluded = []
    invalid = []

    # Load all JSON files in this angle directory
    for i in components:
        for json_path in sorted(glob.glob(os.path.join(angle_dir, i))):
            key = os.path.splitext(os.path.basename(json_path))[0]
            if _is_excluded_component(key):
                excluded.append(key)
                continue
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            if not _is_drawable_component(data):
                invalid.append(key)
                continue
            merged[key] = data

    if not merged:
        raise RuntimeError(f"No component JSON files found in {angle_dir}")

    layout, layout_path = _load_layout_manifest(angle_dir)
    if layout is None:
        # Fallback: keep all debug components in a single material bucket for stable behavior.
        ptfe_sketches = dict(merged)
        electrode_sketches = {}
        xenon_sketches = {}
    else:
        ptfe_sketches, electrode_sketches, xenon_sketches = _bucketize_from_layout(merged, layout)
        print(
            "[debug_autogen] loaded cleanup layout %s (ptfe=%d electrode=%d xenon=%d)"
            % (layout_path, len(ptfe_sketches), len(electrode_sketches), len(xenon_sketches))
        )

    manual_mapping = {}

    return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping
  
