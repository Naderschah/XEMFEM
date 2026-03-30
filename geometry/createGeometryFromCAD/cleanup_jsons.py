"""
This code will handle the merging of repeated elements
"""
import copy
import os, json
import glob
from collections import Counter
from collections import defaultdict
from pprint import pformat
import shutil
import yaml
import sys
def load_yaml_config(cfg_path):
    """
    Loads config.yaml and normalizes/validates the expected structure.

    Features:
      - Expands ~ and environment variables in all string values
      - Resolves relative *paths* relative to the config file's directory
      - Does NOT path-resolve regex strings (e.g. override_periodicity[].match)
      - Ensures optional sections exist with sane defaults:
          override_periodicity: []
          replace_these: []
          ignore_components: []
          replacements: []
          materials: {}
          electrodes: {}
          manual_additions: {}
          slice_input_scale: 1.0
          slice_offset_y: 0.0
          manual_offset_y: 0.0
          shape_match_tol: 1e-6
          periodicity_tol: 1e-4
          grouping_mode: periodic
          slice_alignment: {}
          prune_stale_outputs: True

    Expected minimal config:
      paths:
        in_dir: "..."
        out_dir: "..."
    or batch mode:
      paths:
        in_root: "..."
        out_root: "..."
        slice_glob: "slice_*"
    """
    import os, yaml

    cfg_path = os.path.abspath(os.path.expanduser(os.path.expandvars(cfg_path)))
    if not os.path.isfile(cfg_path):
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    with open(cfg_path, "r", encoding="utf-8") as f:
        try:
            cfg = yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in {cfg_path}: {e}") from e

    if not isinstance(cfg, dict):
        raise ValueError(f"Top-level YAML must be a mapping/object in {cfg_path}")

    base_dir = os.path.dirname(cfg_path)

    def _expand_string(s: str) -> str:
        return os.path.expanduser(os.path.expandvars(s))

    def _resolve_path(p: str) -> str:
        p = _expand_string(p)
        return p if os.path.isabs(p) else os.path.abspath(os.path.join(base_dir, p))

    # Expand env/~ everywhere, but do NOT path-resolve arbitrary strings (regex!)
    def _walk_expand(v):
        if isinstance(v, str):
            return _expand_string(v)
        if isinstance(v, list):
            return [_walk_expand(x) for x in v]
        if isinstance(v, dict):
            return {k: _walk_expand(val) for k, val in v.items()}
        return v

    cfg = _walk_expand(cfg)

    # --- validation / defaults ---
    paths = cfg.get("paths")
    if not isinstance(paths, dict):
        raise ValueError("config.yaml must contain 'paths:' as a mapping/object")

    has_single = ("in_dir" in paths and "out_dir" in paths)
    has_batch = ("in_root" in paths and "out_root" in paths)
    if not has_single and not has_batch:
        raise ValueError(
            "config.yaml paths must contain either "
            "('in_dir' and 'out_dir') or ('in_root' and 'out_root')"
        )

    # Resolve actual path fields relative to config directory
    if "in_dir" in paths and paths["in_dir"] is not None:
        if not isinstance(paths["in_dir"], str) or not paths["in_dir"].strip():
            raise ValueError("config.yaml paths.in_dir must be a non-empty string")
        paths["in_dir"] = _resolve_path(paths["in_dir"])
    if "out_dir" in paths and paths["out_dir"] is not None:
        if not isinstance(paths["out_dir"], str) or not paths["out_dir"].strip():
            raise ValueError("config.yaml paths.out_dir must be a non-empty string")
        paths["out_dir"] = _resolve_path(paths["out_dir"])
    if "in_root" in paths and paths["in_root"] is not None:
        if not isinstance(paths["in_root"], str) or not paths["in_root"].strip():
            raise ValueError("config.yaml paths.in_root must be a non-empty string")
        paths["in_root"] = _resolve_path(paths["in_root"])
    if "out_root" in paths and paths["out_root"] is not None:
        if not isinstance(paths["out_root"], str) or not paths["out_root"].strip():
            raise ValueError("config.yaml paths.out_root must be a non-empty string")
        paths["out_root"] = _resolve_path(paths["out_root"])
    if "slice_glob" in paths and paths["slice_glob"] is not None:
        if not isinstance(paths["slice_glob"], str) or not paths["slice_glob"].strip():
            raise ValueError("config.yaml paths.slice_glob must be a non-empty string")

    if "extras_dir" in paths and paths["extras_dir"] is not None:
        if not isinstance(paths["extras_dir"], str):
            raise ValueError("config.yaml paths.extras_dir must be a string (or omitted)")
        paths["extras_dir"] = _resolve_path(paths["extras_dir"])

    # Optional sections with defaults
    cfg.setdefault("override_periodicity", [])
    cfg.setdefault("replace_these", [])
    cfg.setdefault("ignore_components", [])
    cfg.setdefault("replacements", [])
    cfg.setdefault("materials", {})
    cfg.setdefault("electrodes", {})
    cfg.setdefault("manual_additions", {})
    cfg.setdefault("slice_input_scale", 1.0)
    cfg.setdefault("slice_offset_y", 0.0)
    cfg.setdefault("manual_offset_y", 0.0)
    cfg.setdefault("shape_match_tol", 1e-6)
    cfg.setdefault("periodicity_tol", 1e-4)
    cfg.setdefault("grouping_mode", "periodic")
    cfg.setdefault("slice_alignment", {})
    cfg.setdefault("prune_stale_outputs", True)

    if not isinstance(cfg["override_periodicity"], list):
        raise ValueError("override_periodicity must be a list")
    if not isinstance(cfg["replace_these"], list):
        raise ValueError("replace_these must be a list")
    if not isinstance(cfg["ignore_components"], list):
        raise ValueError("ignore_components must be a list")
    if not isinstance(cfg["replacements"], list):
        raise ValueError("replacements must be a list")
    if not isinstance(cfg["materials"], dict):
        raise ValueError("materials must be a mapping/object")
    if not isinstance(cfg["electrodes"], dict):
        raise ValueError("electrodes must be a mapping/object")
    if not isinstance(cfg["manual_additions"], dict):
        raise ValueError("manual_additions must be a mapping/object")
    if not isinstance(cfg["slice_input_scale"], (int, float)):
        raise ValueError("slice_input_scale must be a number")
    if not isinstance(cfg["slice_offset_y"], (int, float)):
        raise ValueError("slice_offset_y must be a number")
    if not isinstance(cfg["manual_offset_y"], (int, float)):
        raise ValueError("manual_offset_y must be a number")
    if not isinstance(cfg["shape_match_tol"], (int, float)):
        raise ValueError("shape_match_tol must be a number")
    if not isinstance(cfg["periodicity_tol"], (int, float)):
        raise ValueError("periodicity_tol must be a number")
    if not isinstance(cfg["grouping_mode"], str):
        raise ValueError("grouping_mode must be a string")
    if cfg["grouping_mode"] not in ("periodic", "none"):
        raise ValueError("grouping_mode must be 'periodic' or 'none'")
    if not isinstance(cfg["slice_alignment"], dict):
        raise ValueError("slice_alignment must be a mapping/object")
    if not isinstance(cfg["prune_stale_outputs"], bool):
        raise ValueError("prune_stale_outputs must be a boolean")

    if cfg["slice_alignment"]:
        align = cfg["slice_alignment"]
        ref = align.get("reference_component")
        if ref is not None and (not isinstance(ref, str) or not ref.strip()):
            raise ValueError("slice_alignment.reference_component must be a non-empty string")
        target = align.get("reference_top_y", "-liquid_level")
        if not isinstance(target, str):
            raise ValueError("slice_alignment.reference_top_y must be a string")
        default_ll = align.get("default_liquid_level", 0.004)
        if not isinstance(default_ll, (int, float)):
            raise ValueError("slice_alignment.default_liquid_level must be a number")

    # Resolve only known path fields inside sections (leave regex 'match' untouched)
    for rule in cfg.get("replace_these", []) or []:
        if not isinstance(rule, dict):
            continue
        src = rule.get("source")
        if isinstance(src, dict) and isinstance(src.get("path"), str) and src["path"].strip():
            src["path"] = _resolve_path(src["path"])

    for rule in cfg.get("replacements", []) or []:
        if not isinstance(rule, dict):
            continue
        src = rule.get("source")
        if isinstance(src, dict) and isinstance(src.get("path"), str) and src["path"].strip():
            src["path"] = _resolve_path(src["path"])

    manual = cfg.get("manual_additions") or {}
    if isinstance(manual, dict):
        copy_rules = manual.get("copy") or []
        if isinstance(copy_rules, list):
            for rule in copy_rules:
                if not isinstance(rule, dict):
                    continue
                if "offset_y" in rule and not isinstance(rule["offset_y"], (int, float)):
                    raise ValueError("manual_additions.copy[].offset_y must be a number")
                if isinstance(rule.get("from"), str) and rule["from"].strip():
                    rule["from"] = _resolve_path(rule["from"])
                if isinstance(rule.get("from_glob"), str) and rule["from_glob"].strip():
                    rule["from_glob"] = _resolve_path(rule["from_glob"])
                if isinstance(rule.get("to_dir"), str) and rule["to_dir"].strip():
                    rule["to_dir"] = _resolve_path(rule["to_dir"])

    def _merge_slice_filters(existing, new_filters):
        vals = []
        for src in (existing, new_filters):
            if not src:
                continue
            if isinstance(src, str):
                vals.append(src)
            elif isinstance(src, list):
                vals.extend([str(v) for v in src if v is not None and str(v).strip()])
            else:
                vals.append(str(src))
        return vals

    def _normalize_slice_variant_key(section_name, group_name, variant_key):
        if not isinstance(variant_key, str) or not variant_key.strip():
            raise ValueError(f"{section_name}.{group_name} variant keys must be non-empty strings")
        variant_key = variant_key.strip()
        if variant_key.startswith("^"):
            return [variant_key]
        if not variant_key.startswith("slice_"):
            raise ValueError(
                f"{section_name}.{group_name}.{variant_key} is not a valid slice selector. "
                "Use a slice_* shorthand or an explicit regex starting with ^"
            )
        escaped = re.escape(variant_key)
        if variant_key.endswith("deg"):
            return [f"^{escaped}$"]
        return [f"^{escaped}(?:deg)?$"]

    def _normalize_named_component_entry(section_name, group_name, entry, slice_filters=None, priority=0, source_name=None):
        if isinstance(entry, list):
            clean_components = [str(v) for v in entry if isinstance(v, str) and v.strip()]
            return {
                "components": list(dict.fromkeys(clean_components)),
                "default": False,
                "slices": _merge_slice_filters(None, slice_filters),
                "__priority": int(priority),
                "__source_name": source_name or group_name,
            }
        if not isinstance(entry, dict):
            raise ValueError(f"{section_name}.{group_name} must be a mapping or list")
        entry = copy.deepcopy(entry)
        if "slice" in entry and "slices" not in entry:
            entry["slices"] = _merge_slice_filters(None, entry["slice"])
            entry.pop("slice", None)
        components = entry.get("components", [])
        if components is None:
            components = []
        if not isinstance(components, list):
            raise ValueError(f"{section_name}.{group_name}.components must be a list")
        clean_components = []
        for component_name in components:
            if not isinstance(component_name, str) or not component_name.strip():
                raise ValueError(
                    f"{section_name}.{group_name}.components contains invalid name: {component_name!r}"
                )
            clean_components.append(component_name)
        return {
            "components": list(dict.fromkeys(clean_components)),
            "default": bool(entry.get("default", False)),
            "slices": _merge_slice_filters(slice_filters, entry.get("slices")),
            "__priority": int(priority),
            "__source_name": source_name or group_name,
        }

    def _normalize_named_component_section(section_name):
        normalized = {}
        reserved_keys = {"components", "default", "slice", "slices"}
        for group_name, entry in (cfg.get(section_name) or {}).items():
            if not isinstance(group_name, str) or not group_name.strip():
                raise ValueError(f"{section_name} keys must be non-empty strings")
            if isinstance(entry, list):
                normalized[group_name] = [
                    _normalize_named_component_entry(section_name, group_name, entry)
                ]
                continue
            if not isinstance(entry, dict):
                raise ValueError(f"{section_name}.{group_name} must be a mapping or list")
            variant_keys = [k for k in entry.keys() if k not in reserved_keys]
            if not variant_keys:
                normalized[group_name] = [
                    _normalize_named_component_entry(section_name, group_name, entry)
                ]
                continue

            variants = []
            base_entry = {k: copy.deepcopy(v) for k, v in entry.items() if k in reserved_keys}
            if base_entry:
                variants.append(
                    _normalize_named_component_entry(
                        section_name,
                        group_name,
                        base_entry,
                        priority=0,
                        source_name=f"{group_name} (default)",
                    )
                )

            for variant_key in variant_keys:
                variant_entry = entry[variant_key]
                if isinstance(variant_entry, dict) and ("slice" in variant_entry or "slices" in variant_entry):
                    raise ValueError(
                        f"{section_name}.{group_name}.{variant_key} must not define slice/slices; "
                        "the nested key already selects the slice"
                    )
                variants.append(
                    _normalize_named_component_entry(
                        section_name,
                        group_name,
                        variant_entry,
                        slice_filters=_normalize_slice_variant_key(section_name, group_name, variant_key),
                        priority=1,
                        source_name=f"{group_name}.{variant_key}",
                    )
                )
            normalized[group_name] = variants
        cfg[section_name] = normalized

    _normalize_named_component_section("materials")
    _normalize_named_component_section("electrodes")

    return cfg

def _scale_pts_instruction(instr, scale):
  if not isinstance(instr, list) or not instr:
    return instr
  tag = instr[0]
  out = copy.deepcopy(instr)
  if len(out) > 1 and isinstance(out[1], (int, float)):
    out[1] = float(out[1]) * scale
  if len(out) > 2 and isinstance(out[2], (int, float)):
    out[2] = float(out[2]) * scale
  if tag == "arc":
    if len(out) > 3 and isinstance(out[3], (int, float)):
      out[3] = float(out[3]) * scale
    if len(out) > 4 and isinstance(out[4], (int, float)):
      out[4] = float(out[4]) * scale
  elif tag == "ellipse":
    if len(out) > 3 and isinstance(out[3], (int, float)):
      out[3] = float(out[3]) * scale
    if len(out) > 4 and isinstance(out[4], (int, float)):
      out[4] = float(out[4]) * scale
    if len(out) > 5 and isinstance(out[5], (int, float)):
      out[5] = float(out[5]) * scale
    if len(out) > 6 and isinstance(out[6], (int, float)):
      out[6] = float(out[6]) * scale
  elif tag == "spline" and len(out) > 3 and isinstance(out[3], dict):
    payload = copy.deepcopy(out[3])
    poles = payload.get("poles")
    if isinstance(poles, list):
      payload["poles"] = [
        [float(p[0]) * scale, float(p[1]) * scale]
        for p in poles
        if isinstance(p, (list, tuple)) and len(p) >= 2
      ]
    out[3] = payload
  return out


def _scale_slice_component(data, scale):
  if scale == 1.0:
    return data
  if not isinstance(data, dict):
    return data
  out = copy.deepcopy(data)
  pts = out.get("pts")
  if isinstance(pts, list):
    out["pts"] = [_scale_pts_instruction(item, scale) for item in pts]
  for key in ("HorizontalPitch", "VerticalPitch"):
    if isinstance(out.get(key), (int, float)):
      out[key] = float(out[key]) * scale
  return out


def _translate_pts_instruction(instr, dx=0.0, dy=0.0):
  if not isinstance(instr, list) or not instr:
    return instr
  tag = instr[0]
  out = copy.deepcopy(instr)
  if len(out) > 1 and isinstance(out[1], (int, float)):
    out[1] = float(out[1]) + dx
  if len(out) > 2 and isinstance(out[2], (int, float)):
    out[2] = float(out[2]) + dy
  if tag == "arc":
    if len(out) > 3 and isinstance(out[3], (int, float)):
      out[3] = float(out[3]) + dx
    if len(out) > 4 and isinstance(out[4], (int, float)):
      out[4] = float(out[4]) + dy
  elif tag == "ellipse":
    if len(out) > 3 and isinstance(out[3], (int, float)):
      out[3] = float(out[3]) + dx
    if len(out) > 4 and isinstance(out[4], (int, float)):
      out[4] = float(out[4]) + dy
  elif tag == "spline" and len(out) > 3 and isinstance(out[3], dict):
    payload = copy.deepcopy(out[3])
    poles = payload.get("poles")
    if isinstance(poles, list):
      new_poles = []
      for p in poles:
        if isinstance(p, (list, tuple)) and len(p) >= 2:
          new_poles.append([float(p[0]) + dx, float(p[1]) + dy])
        else:
          new_poles.append(p)
      payload["poles"] = new_poles
    out[3] = payload
  return out


def _translate_component(data, dx=0.0, dy=0.0):
  if (dx == 0.0 and dy == 0.0) or not isinstance(data, dict):
    return data
  out = copy.deepcopy(data)
  pts = out.get("pts")
  if isinstance(pts, list):
    out["pts"] = [_translate_pts_instruction(item, dx, dy) for item in pts]
  if isinstance(out.get("RadialPosition"), (int, float)):
    out["RadialPosition"] = float(out["RadialPosition"]) + dx
  if isinstance(out.get("VerticalPosition"), (int, float)):
    out["VerticalPosition"] = float(out["VerticalPosition"]) + dy
  sub = out.get("sub_sketches")
  if isinstance(sub, dict):
    out["sub_sketches"] = {
      key: _translate_component(value, dx, dy)
      for key, value in sub.items()
    }
  elif isinstance(sub, list):
    out["sub_sketches"] = [_translate_component(value, dx, dy) for value in sub]
  return out


def _offset_component_y(data, offset_y):
  return _translate_component(data, 0.0, offset_y)

def _component_points(data):
  points = []
  if not isinstance(data, dict):
    return points

  pts = data.get("pts")
  if isinstance(pts, list):
    for row in pts:
      if isinstance(row, list) and len(row) > 2 and isinstance(row[1], (int, float)) and isinstance(row[2], (int, float)):
        points.append((float(row[1]), float(row[2])))
      if isinstance(row, list) and row and row[0] == "spline" and len(row) > 3 and isinstance(row[3], dict):
        poles = row[3].get("poles")
        if isinstance(poles, list):
          for p in poles:
            if isinstance(p, (list, tuple)) and len(p) >= 2 and isinstance(p[0], (int, float)) and isinstance(p[1], (int, float)):
              points.append((float(p[0]), float(p[1])))

  if isinstance(data.get("RadialPosition"), (int, float)) and isinstance(data.get("VerticalPosition"), (int, float)):
    points.append((float(data["RadialPosition"]), float(data["VerticalPosition"])))

  sub = data.get("sub_sketches")
  if isinstance(sub, dict):
    for child in sub.values():
      points.extend(_component_points(child))
  elif isinstance(sub, list):
    for child in sub:
      points.extend(_component_points(child))

  return points


def _component_anchor(data):
  if not isinstance(data, dict):
    return 0.0, 0.0
  points = _component_points(data)
  if points:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), min(ys)
  x = float(data.get("RadialPosition", 0.0)) if isinstance(data.get("RadialPosition"), (int, float)) else 0.0
  y = float(data.get("VerticalPosition", 0.0)) if isinstance(data.get("VerticalPosition"), (int, float)) else 0.0
  return x, y


def _rebuild_file_meta(all_files):
  return {
      name: {
          "n_points": len(all_files[name].get("pts", [])),
          "min_x": _component_anchor(all_files[name])[0],
          "min_y": _component_anchor(all_files[name])[1],
      }
      for name in all_files.keys()
  }


def load_and_prepare_jsons(path, slice_input_scale=1.0, slice_offset_y=0.0):
  """
  Takes the path to one slice directory
  returns all objects in a json and a few precomputed things
  """
  all_files = {}

  # Deterministic ordering avoids run-to-run differences in which target names
  # become grouping anchors.
  for fname in sorted(os.listdir(path)):
    if not fname.lower().endswith(".json"):
      continue

    # your naming rule: drop the last two extensions (e.g. NAME.abc.json -> NAME)
    key = ".".join(fname.split(".")[:-2]) if len(fname.split(".")) >= 3 else os.path.splitext(fname)[0]
    fpath = os.path.join(path, fname)

    with open(fpath, "r", encoding="utf-8") as f:
      payload = _scale_slice_component(json.load(f), slice_input_scale)
      all_files[key] = _offset_component_y(payload, slice_offset_y)

  file_meta = _rebuild_file_meta(all_files)
  return all_files, file_meta


def get_slice_name_from_dir(path):
  return os.path.basename(os.path.normpath(os.path.abspath(path)))


def resolve_slice_jobs(paths):
  """
  Returns a list of (in_dir, out_dir) pairs.
  Supports:
    - single-slice mode via paths.in_dir + paths.out_dir
    - batch mode via paths.in_root + paths.out_root (+ optional paths.slice_glob)
  """
  in_root = paths.get("in_root")
  out_root = paths.get("out_root")
  if in_root and out_root:
    slice_glob = paths.get("slice_glob", "slice_*")
    pattern = os.path.join(in_root, slice_glob)
    in_dirs = sorted(
      d for d in glob.glob(pattern)
      if os.path.isdir(d)
    )
    if not in_dirs:
      raise FileNotFoundError(f"No slice directories matched {pattern}")
    return [
      (in_dir, os.path.join(out_root, get_slice_name_from_dir(in_dir)))
      for in_dir in in_dirs
    ]

  return [(paths["in_dir"], paths["out_dir"])]


def _rule_slices(rule):
  if not isinstance(rule, dict):
    return []
  slices = rule.get("slices")
  if slices is None:
    slices = rule.get("slice")
  if slices is None:
    return []
  if isinstance(slices, str):
    return [slices]
  if isinstance(slices, list):
    return [str(v) for v in slices if v is not None and str(v).strip()]
  return [str(slices)]


def rule_applies_to_slice(rule, slice_name):
  pats = _rule_slices(rule)
  if not pats:
    return True
  for pat in pats:
    if re.search(pat, slice_name):
      return True
  return False


def applicable_named_component_section(cfg, key, slice_name):
  section = cfg.get(key) or {}
  out = {}
  for name, entries in section.items():
    if isinstance(entries, dict):
      entries = [entries]
    if not isinstance(entries, list):
      continue
    matches = []
    for entry in entries:
      if not isinstance(entry, dict):
        continue
      if rule_applies_to_slice(entry, slice_name):
        matches.append(entry)
    if not matches:
      continue
    best_priority = max(int(entry.get("__priority", 0)) for entry in matches)
    best_matches = [entry for entry in matches if int(entry.get("__priority", 0)) == best_priority]
    if len(best_matches) > 1:
      source_names = [str(entry.get("__source_name", name)) for entry in best_matches]
      raise ValueError(
        "multiple %s variants matched slice %s for %s: %s"
        % (key, slice_name, name, ", ".join(sorted(source_names)))
      )
    chosen = copy.deepcopy(best_matches[0])
    chosen.pop("__priority", None)
    chosen.pop("__source_name", None)
    out[name] = chosen
  return out


def apply_ignore_components(all_files, file_meta, ignore_rules, slice_name):
  if not ignore_rules:
    return all_files, file_meta, []
  removed = []
  for rule in ignore_rules:
    if not isinstance(rule, dict):
      continue
    if not rule_applies_to_slice(rule, slice_name):
      continue
    pat = rule.get("match")
    if not pat:
      continue
    rx = re.compile(pat)
    matched = [name for name in list(all_files.keys()) if rx.search(name)]
    for name in matched:
      all_files.pop(name, None)
      file_meta.pop(name, None)
      removed.append(name)
  return all_files, file_meta, removed


def _load_component_from_source(source, all_files, in_dir, slice_input_scale, slice_offset_y):
  if isinstance(source, str):
    if source in all_files:
      return copy.deepcopy(all_files[source])
    payload = _scale_slice_component(
      _load_json_file(os.path.join(in_dir, "%s.json" % source)),
      slice_input_scale,
    )
    return _offset_component_y(payload, slice_offset_y)
  if isinstance(source, dict):
    src_path = source.get("path")
    if not src_path:
      raise ValueError(f"replacement source dict missing path: {source!r}")
    return copy.deepcopy(_load_json_file(src_path))
  raise ValueError(f"unsupported replacement source: {source!r}")


def apply_replacements(all_files, file_meta, replacement_rules, slice_name, in_dir, slice_input_scale, slice_offset_y):
  if not replacement_rules:
    return all_files, file_meta, []
  applied = []
  for rule in replacement_rules:
    if not isinstance(rule, dict):
      continue
    if not rule_applies_to_slice(rule, slice_name):
      continue
    targets = []
    if isinstance(rule.get("target"), str) and rule["target"].strip():
      targets.append(rule["target"])
    if isinstance(rule.get("targets"), list):
      targets.extend([str(v) for v in rule["targets"] if v is not None and str(v).strip()])
    targets = list(dict.fromkeys(targets))
    if not targets:
      raise ValueError(f"replacement rule missing target/targets: {rule!r}")
    missing = [name for name in targets if name not in all_files]
    if missing:
      raise KeyError("replacement targets missing from slice %s: %s" % (slice_name, ", ".join(sorted(missing))))
    source = rule.get("source")
    if source is None:
      raise ValueError(f"replacement rule missing source: {rule!r}")
    output_name = rule.get("output_name") or targets[0]
    source_component = _load_component_from_source(source, all_files, in_dir, slice_input_scale, slice_offset_y)
    target_min_x = min(file_meta[name]["min_x"] for name in targets)
    target_min_y = min(file_meta[name]["min_y"] for name in targets)
    source_min_x, source_min_y = _component_anchor(source_component)
    replacement_component = _translate_component(
      source_component,
      dx=target_min_x - source_min_x,
      dy=target_min_y - source_min_y,
    )
    for name in targets:
      all_files.pop(name, None)
      file_meta.pop(name, None)
    all_files[output_name] = replacement_component
    file_meta[output_name] = {
      "n_points": len(replacement_component.get("pts", [])),
      "min_x": _component_anchor(replacement_component)[0],
      "min_y": _component_anchor(replacement_component)[1],
    }
    applied.append({
      "output_name": output_name,
      "targets": targets,
      "source": source if isinstance(source, str) else source.get("path"),
    })
  return all_files, file_meta, applied

def get_candidates_pre(file_meta, target):
  """
  Finds all candidates based on n_points
  returned as list of strings
  """
  target_family = _group_family_key(target)
  return [
    i
    for i in file_meta.keys()
    if (_group_family_key(i) == target_family)
      and (file_meta[target]['n_points'] == file_meta[i]['n_points'])
      and (((file_meta[target]['min_x'] - file_meta[i]['min_x']) != 0.0) != 
            ((file_meta[target]['min_y'] - file_meta[i]['min_y']) != 0.0))
  ]

def _group_family_key(name):
  """
  Grouping is restricted to the component stem prior to the final subindex.
  Example:
    XENTTPCBWarmPart2260Part2260_1 -> XENTTPCBWarmPart2260Part2260
  """
  if not isinstance(name, str):
    return name
  m = re.match(r"^(.*)_\d+$", name)
  if m:
    return m.group(1)
  return name

def _quantize(v: float, tol: float) -> int:
    # maps v to an integer grid index; values within ~tol fall into same bucket
    return int(round(v / tol))

def _normalize_pts_instruction(instr, min_x: float, min_y: float, tol: float):
    if not isinstance(instr, list) or not instr:
        return instr

    tag = instr[0]
    if tag == "line":
        return (
            "line",
            _quantize(instr[1] - min_x, tol),
            _quantize(instr[2] - min_y, tol),
        )

    if tag == "arc":
        base = [
            "arc",
            _quantize(instr[1] - min_x, tol),
            _quantize(instr[2] - min_y, tol),
            _quantize(instr[3] - min_x, tol),
            _quantize(instr[4] - min_y, tol),
        ]
        if len(instr) > 5:
            base.append(instr[5])
        return tuple(base)

    if tag == "ellipse":
        base = ["ellipse"]
        numeric_pairs = {
            1: min_x,
            2: min_y,
            3: min_x,
            4: min_y,
            5: min_x,
            6: min_y,
        }
        for idx, value in enumerate(instr[1:], start=1):
            if isinstance(value, (int, float)) and idx in numeric_pairs:
                base.append(_quantize(value - numeric_pairs[idx], tol))
            else:
                base.append(value)
        return tuple(base)

    if tag == "spline":
        payload = instr[3] if len(instr) > 3 and isinstance(instr[3], dict) else {}
        poles = tuple(
            (
                _quantize(p[0] - min_x, tol),
                _quantize(p[1] - min_y, tol),
            )
            for p in payload.get("poles", [])
            if isinstance(p, (list, tuple)) and len(p) >= 2
        )
        weights = tuple(
            _quantize(float(w), tol) if isinstance(w, (int, float)) else w
            for w in payload.get("weights", [])
        )
        knots = tuple(
            _quantize(float(k), tol) if isinstance(k, (int, float)) else k
            for k in payload.get("knots", [])
        )
        mults = tuple(payload.get("mults", []))
        return (
            "spline",
            _quantize(instr[1] - min_x, tol),
            _quantize(instr[2] - min_y, tol),
            (
                ("degree", payload.get("degree")),
                ("flags", payload.get("flags")),
                ("knots", knots),
                ("mults", mults),
                ("periodic", payload.get("periodic")),
                ("poles", poles),
                ("weights", weights),
            ),
        )

    out = [tag]
    for value in instr[1:]:
        if isinstance(value, (int, float)):
            out.append(_quantize(float(value), tol))
        else:
            out.append(value)
    return tuple(out)

def _normalize_pts(pts, min_x: float, min_y: float, tol: float):
    out = []
    for j in pts:
        out.append(_normalize_pts_instruction(j, min_x, min_y, tol))
    return out

def filter_candidates(all_files, candidates, file_meta, target, tol=1e-6):
    target_keys = _normalize_pts(
        all_files[target]["pts"],
        file_meta[target]["min_x"],
        file_meta[target]["min_y"],
        tol,
    )
    target_counter = Counter(target_keys)

    kept = []
    for i in candidates:
        cand_keys = _normalize_pts(
            all_files[i]["pts"],
            file_meta[i]["min_x"],
            file_meta[i]["min_y"],
            tol,
        )
        if Counter(cand_keys) == target_counter:
            kept.append(i)
    return kept

def _quantize(v: float, tol: float) -> int:
    return int(round(v / tol)) if tol > 0 else int(v)

def _find_axis_periodic_groups(file_meta, candidates, axis: str,
                               tol_pos=1e-6, tol_step=1e-6, min_len=3):
  key = f"min_{axis}"

  # keep (coord, candidate_name)
  items = sorted((file_meta[name][key], name) for name in candidates)
  if len(items) < 2:
    return []

  coords = [c for c, _ in items]
  names  = [n for _, n in items]
  diffs  = [coords[k+1] - coords[k] for k in range(len(coords) - 1)]

  groups = []
  start = 0

  while start < len(names) - 1:
    if abs(diffs[start]) <= tol_pos:
      start += 1
      continue

    step_bin = _quantize(diffs[start], tol_step)
    end = start

    while (
        end < len(diffs)
        and abs(diffs[end]) > tol_pos
        and _quantize(diffs[end], tol_step) == step_bin
    ):
      end += 1

    if end - start + 1 >= min_len:
      groups.append({
        "axis": axis,
        "step": diffs[start],
        "members": [
          {
            "name": names[k],
            key: coords[k],
            "index": k,
          }
          for k in range(start, end + 1)
        ],
      })

    start = end

  return groups


def find_periodic_groups(file_meta, candidates,
                          tol_pos=1e-6, tol_step=1e-6, min_len=3):
  return {
      "x": _find_axis_periodic_groups(file_meta, candidates, "x",
                                      tol_pos, tol_step, min_len),
      "y": _find_axis_periodic_groups(file_meta, candidates, "y",
                                      tol_pos, tol_step, min_len),
  }


def _best_group_and_axis(groupings):
  best_axis = None
  best_group = None
  best_len = 0
  for axis in ("x", "y"):
    for g in ((groupings or {}).get(axis) or []):
      L = len(g.get("members", []))
      if L > best_len:
        best_len = L
        best_axis = axis
        best_group = g
  return best_axis, best_group, best_len


def _canonical_group_target(groupings):
  axis, group, _ = _best_group_and_axis(groupings)
  if not group:
    return None
  return _canonical_member_for_group(group, axis)

def _canonical_member_for_group(group, axis):
  coord_key = "min_x" if axis == "x" else "min_y"
  members = list(group.get("members", []))
  if not members:
    return None
  if axis == "y":
    members = sorted(members, key=lambda m: m[coord_key], reverse=True)
  else:
    members = sorted(members, key=lambda m: m[coord_key])
  return members[0]["name"]

def _group_signature(axis, group):
  return (
    axis,
    tuple(sorted(m["name"] for m in group.get("members", []))),
  )

import os
import json
import re

def _first_matching_rule(rules, target):
  for rule in (rules or []):
    pat = rule.get("match")
    if not pat:
      continue
    if re.search(pat, target):
      return rule
  return None

def _load_json_file(path):
  with open(path, "r", encoding="utf-8") as f:
    return json.load(f)

def _resolve_override_for_target(target, groupings, override_rules):
  pre = (groupings or {}).get("_pre_override")
  if isinstance(pre, dict):
    resolved = pre.get("resolved")
    if isinstance(resolved, dict):
      axis = resolved.get("axis")
      pitch = resolved.get("pitch")
      number = resolved.get("number")
      if axis in ("x", "y") and pitch is not None and number is not None:
        return {"axis": axis, "pitch": float(pitch), "number": int(number)}

  rule = _first_matching_rule(override_rules, target)
  if not rule:
    return None

  set_cfg = rule.get("set") or {}
  axis = set_cfg.get("axis")
  pitch = set_cfg.get("pitch")
  number = set_cfg.get("number")

  if axis not in ("x", "y") or number is None:
    raise ValueError("override_periodicity rule missing axis/number for target=%s: %r" % (target, rule))

  if pitch is None:
    raise ValueError(
      "override_periodicity for target=%s requires pitch here; "
      "or run apply_override_periodicity_pre to resolve it." % target
    )

  return {"axis": axis, "pitch": float(pitch), "number": int(number)}

def make_and_write_grouped(
  target,
  groupings,
  in_dir,
  out_dir,
  tol=1e-6,
  replace_rules=None,
  override_rules=None,
  slice_input_scale=1.0,
  slice_offset_y=0.0,
):
  os.makedirs(out_dir, exist_ok=True)

  def load_by_name(name):
    payload = _scale_slice_component(
      _load_json_file(os.path.join(in_dir, "%s.json" % name)),
      slice_input_scale,
    )
    return _offset_component_y(payload, slice_offset_y)

  override = _resolve_override_for_target(target, groupings, override_rules or [])

  best_axis, best_group, best_len = _best_group_and_axis(groupings)

  if override is None:
    if best_group is None or best_len < 2:
      raise ValueError("No usable periodic grouping for target=%s" % target)
    axis = best_axis
  else:
    axis = override["axis"]

  base_name = target
  if best_group is not None:
    coord_key = "min_x" if axis == "x" else "min_y"
    members = best_group.get("members", [])
    if members:
      members_sorted = sorted(
        members,
        key=lambda m: m.get(coord_key, 0.0),
        reverse=(axis == "y"),
      )
      base_name = members_sorted[0]["name"]

  base_data = load_by_name(base_name)

  rep_rule = _first_matching_rule(replace_rules or [], target)
  if rep_rule:
    src = rep_rule.get("source") or {}
    src_path = src.get("path")
    if not src_path:
      raise ValueError("replace_these rule for target=%s missing source.path: %r" % (target, rep_rule))
    template = _load_json_file(src_path)
    base_data["pts"] = template.get("pts", [])

  pts = base_data.get("pts", [])

  if override is not None:
    pitch = float(override["pitch"])
    number = int(override["number"])
  else:
    coord_key = "min_x" if axis == "x" else "min_y"
    members = best_group["members"]
    members_sorted = sorted(members, key=lambda m: m[coord_key], reverse=(axis == "y"))

    if len(members_sorted) < 2:
      raise ValueError("Grouping too small to compute pitch for target=%s" % target)

    pitch = members_sorted[1][coord_key] - members_sorted[0][coord_key]
    number = len(members_sorted)

  out = {
    "pts": pts,
    ("HorizontalPitch" if axis == "x" else "VerticalPitch"): pitch,
    "Number": number,
  }

  out_path = os.path.join(out_dir, "%s_grouped.json" % target)
  with open(out_path, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
  return out_path

def apply_override_periodicity_pre(all_files, file_meta, override_rules, tol=1e-6):
  forced_groupings = {}
  forced_used = set()

  if not override_rules:
      return all_files, file_meta, forced_groupings, forced_used

  names = list(all_files.keys())

  for rule in override_rules:
    pattern = rule["match"]
    rx = re.compile(pattern)

    use_source = rule.get("use_source", "min_extremum")
    set_cfg = rule.get("set", {})
    axis = set_cfg.get("axis")
    if axis not in ("x", "y"):
      raise ValueError(f"override_periodicity rule axis must be 'x' or 'y': {rule}")

    coord_key = "min_x" if axis == "x" else "min_y"

    matched = [n for n in names if rx.search(n)]
    if not matched:
      raise Exception("Nothing matched reges %r for override_periodicity: %r" % (pattern, rule))
      continue

    # Resolve pitch if not provided
    pitch = set_cfg.get("pitch", None)
    if pitch is None:
      coords = sorted(file_meta[n][coord_key] for n in matched if n in file_meta)
      diffs = [coords[i + 1] - coords[i] for i in range(len(coords) - 1)]
      pos_diffs = [d for d in diffs if d > tol]
      pitch = min(pos_diffs) if pos_diffs else 0.0

    pitch = abs(pitch)
    if use_source == "max_extremum":
      pitch = -pitch

    number = set_cfg.get("number")
    if number is None:
      raise ValueError(f"override_periodicity requires 'number': {rule}")

    # --- pick ONE representative to emit ---
    # Use min_extremum => smallest coord, max_extremum => largest coord
    present = [n for n in matched if n in file_meta]
    if not present:
      continue

    rep = min(present, key=lambda n: file_meta[n][coord_key]) \
      if use_source == "min_extremum" \
      else max(present, key=lambda n: file_meta[n][coord_key])

    # Remove ALL matched from downstream, but only emit rep
    for n in matched:
      if n in all_files:
        forced_used.add(n)
        all_files.pop(n, None)
        file_meta.pop(n, None)

    forced_groupings[rep] = {
        "x": [],
        "y": [],
        "_pre_override": {
            "match": pattern,
            "use_source": use_source,
            "resolved": {"axis": axis, "pitch": pitch, "number": number},
        },
    }

  return all_files, file_meta, forced_groupings, forced_used
import glob
import os
import shutil

def apply_manual_additions(manual_additions_cfg, out_dir):
  """
  Copies extra JSON files into out_dir regardless of grouping.

  Expected config shape:
    manual_additions:
      copy:
        - from: "/abs/path/A.json"
          to_name: "A_OUT"            # optional; default = stem of from
        - from_glob: "/abs/path/*.json"
          # to_dir: "/abs/path/out"   # optional; default = out_dir
  """
  if not manual_additions_cfg:
    return []

  copy_rules = manual_additions_cfg.get("copy", [])
  if not copy_rules:
    return []

  os.makedirs(out_dir, exist_ok=True)
  written_paths = []
  default_offset_y = float(manual_additions_cfg.get("_default_offset_y", 0.0))

  def _copy_with_optional_offset(src, dst, offset_y):
    if offset_y == 0.0:
      shutil.copy2(src, dst)
      return
    payload = _load_json_file(src)
    payload = _offset_component_y(payload, offset_y)
    with open(dst, "w", encoding="utf-8") as f:
      json.dump(payload, f, ensure_ascii=False, indent=2)

  for rule in copy_rules:
    # destination directory (per rule override)
    dst_dir = rule.get("to_dir", out_dir)
    offset_y = float(rule.get("offset_y", default_offset_y))

    if "from" in rule:
      src = rule["from"]
      stem = rule.get("to_name") or os.path.splitext(os.path.basename(src))[0]
      dst = os.path.join(dst_dir, f"{stem}.json")
      os.makedirs(dst_dir, exist_ok=True)
      _copy_with_optional_offset(src, dst, offset_y)
      written_paths.append(dst)
      continue

    if "from_glob" in rule:
      pattern = rule["from_glob"]
      matches = sorted(glob.glob(pattern))
      os.makedirs(dst_dir, exist_ok=True)
      for src in matches:
        stem = os.path.splitext(os.path.basename(src))[0]
        dst = os.path.join(dst_dir, f"{stem}.json")
        _copy_with_optional_offset(src, dst, offset_y)
        written_paths.append(dst)
      continue

    raise ValueError(f"manual_additions.copy entry must have 'from' or 'from_glob': {rule}")
  return written_paths


def _load_written_components_for_dir(written_paths, out_dir):
  out_dir_abs = os.path.abspath(out_dir)
  component_store = {}
  for path in sorted(written_paths):
    abs_path = os.path.abspath(path)
    if os.path.dirname(abs_path) != out_dir_abs:
      continue
    if not abs_path.lower().endswith(".json"):
      continue
    name = os.path.splitext(os.path.basename(abs_path))[0]
    component_store[name] = _load_json_file(abs_path)
  return component_store


def _clone_named_components(component_store, names, owner_name):
  present = {}
  missing = []
  for name in names:
    if name in component_store:
      present[name] = copy.deepcopy(component_store[name])
    else:
      missing.append(name)
  if missing:
    print(
      "[cleanup_jsons warning] %s references missing components: %s"
      % (owner_name, ", ".join(sorted(missing)))
    )
  return present


def _materialize_sketch_buckets(component_store, materials, electrodes_cfg):
  used_names = set()
  xenon_sketches = {}
  electrode_sketches = {}

  supported_materials = {"GXe", "LXe", "PTFE"}
  for material_name, entry in materials.items():
    if material_name not in supported_materials:
      raise ValueError(
        "Unsupported material in cleanup config: %s. Supported: GXe, LXe, PTFE" % material_name
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
      if components[0] not in component_store:
        print(
          "[cleanup_jsons warning] material %s references missing component: %s"
          % (material_name, components[0])
        )
        continue
      xenon_sketches[material_name] = copy.deepcopy(component_store[components[0]])
      used_names.add(components[0])
    else:
      sub_sketches = _clone_named_components(
        component_store,
        components,
        f"material {material_name}",
      )
      if not sub_sketches:
        print(
          "[cleanup_jsons warning] material %s has no available components; skipping"
          % material_name
        )
        continue
      xenon_sketches[material_name] = {
        "hull": True,
        "sub_sketches": sub_sketches,
      }
      used_names.update(sub_sketches.keys())

  for electrode_name, entry in electrodes_cfg.items():
    if not isinstance(entry, dict):
      raise ValueError(f"Invalid electrode entry for {electrode_name}: {entry!r}")
    components = list(entry.get("components", []) or [])
    if not components:
      continue
    sub_sketches = _clone_named_components(
      component_store,
      components,
      f"electrode {electrode_name}",
    )
    if not sub_sketches:
      print(
        "[cleanup_jsons warning] electrode %s has no available components; skipping"
        % electrode_name
      )
      continue
    electrode_sketches[electrode_name] = {
      "hull": True,
      "sub_sketches": sub_sketches,
    }
    used_names.update(sub_sketches.keys())

  ptfe_sketches = {}
  for name, data in component_store.items():
    if name in used_names:
      continue
    ptfe_sketches[name] = copy.deepcopy(data)

  return ptfe_sketches, electrode_sketches, xenon_sketches


def write_salome_python_config(cfg, slice_name, written_paths, out_dir, manual_written_paths=None):
  component_store = _load_written_components_for_dir(written_paths, out_dir)
  materials = applicable_named_component_section(cfg, "materials", slice_name)
  electrodes_cfg = applicable_named_component_section(cfg, "electrodes", slice_name)
  slice_alignment = copy.deepcopy(cfg.get("slice_alignment") or {})
  if slice_alignment and not rule_applies_to_slice(slice_alignment, slice_name):
    slice_alignment = {}
  if slice_alignment:
    slice_alignment.setdefault("reference_top_y", "-liquid_level")
    slice_alignment.setdefault("default_liquid_level", 0.004)
  manual_component_names = sorted({
    os.path.splitext(os.path.basename(path))[0]
    for path in (manual_written_paths or [])
    if os.path.abspath(os.path.dirname(path)) == os.path.abspath(out_dir)
       and path.lower().endswith(".json")
  })

  salome_cfg_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../XEMFEM_geometries")
  os.makedirs(salome_cfg_dir, exist_ok=True)
  module_path = os.path.join(salome_cfg_dir, f"{slice_name}.py")
  lines = [
    '"""Auto-generated by cleanup_jsons.py from cleaned slice outputs."""',
    "import copy",
    "",
    "liquid_level = " + repr(float(slice_alignment.get("default_liquid_level", 0.004))),
    "",
    "_COMPONENT_STORE = " + pformat(component_store, sort_dicts=False, width=100),
    "",
    "_MATERIALS = " + pformat(materials, sort_dicts=False, width=100),
    "",
    "_ELECTRODES = " + pformat(electrodes_cfg, sort_dicts=False, width=100),
    "",
    "_MANUAL_COMPONENT_NAMES = " + pformat(manual_component_names, sort_dicts=False, width=100),
    "",
    "_SLICE_ALIGNMENT = " + pformat(slice_alignment, sort_dicts=False, width=100),
    "",
    "_MANUAL_MAPPING = {}",
    "",
    "def _translate_pts_instruction(instr, dx=0.0, dy=0.0):",
    "  if not isinstance(instr, list) or not instr:",
    "    return instr",
    "  tag = instr[0]",
    "  out = copy.deepcopy(instr)",
    "  if len(out) > 1 and isinstance(out[1], (int, float)):",
    "    out[1] = float(out[1]) + dx",
    "  if len(out) > 2 and isinstance(out[2], (int, float)):",
    "    out[2] = float(out[2]) + dy",
    "  if tag == 'arc':",
    "    if len(out) > 3 and isinstance(out[3], (int, float)):",
    "      out[3] = float(out[3]) + dx",
    "    if len(out) > 4 and isinstance(out[4], (int, float)):",
    "      out[4] = float(out[4]) + dy",
    "  elif tag == 'ellipse':",
    "    if len(out) > 3 and isinstance(out[3], (int, float)):",
    "      out[3] = float(out[3]) + dx",
    "    if len(out) > 4 and isinstance(out[4], (int, float)):",
    "      out[4] = float(out[4]) + dy",
    "  elif tag == 'spline' and len(out) > 3 and isinstance(out[3], dict):",
    "    payload = copy.deepcopy(out[3])",
    "    poles = payload.get('poles')",
    "    if isinstance(poles, list):",
    "      new_poles = []",
    "      for p in poles:",
    "        if isinstance(p, (list, tuple)) and len(p) >= 2:",
    "          new_poles.append([float(p[0]) + dx, float(p[1]) + dy])",
    "        else:",
    "          new_poles.append(p)",
    "      payload['poles'] = new_poles",
    "    out[3] = payload",
    "  return out",
    "",
    "def _translate_component(data, dx=0.0, dy=0.0):",
    "  if (dx == 0.0 and dy == 0.0) or not isinstance(data, dict):",
    "    return data",
    "  out = copy.deepcopy(data)",
    "  pts = out.get('pts')",
    "  if isinstance(pts, list):",
    "    out['pts'] = [_translate_pts_instruction(item, dx, dy) for item in pts]",
    "  if isinstance(out.get('RadialPosition'), (int, float)):",
    "    out['RadialPosition'] = float(out['RadialPosition']) + dx",
    "  if isinstance(out.get('VerticalPosition'), (int, float)):",
    "    out['VerticalPosition'] = float(out['VerticalPosition']) + dy",
    "  sub = out.get('sub_sketches')",
    "  if isinstance(sub, dict):",
    "    out['sub_sketches'] = {k: _translate_component(v, dx, dy) for k, v in sub.items()}",
    "  elif isinstance(sub, list):",
    "    out['sub_sketches'] = [_translate_component(v, dx, dy) for v in sub]",
    "  return out",
    "",
    "def _component_top_y(data):",
    "  vals = []",
    "  if not isinstance(data, dict):",
    "    return None",
    "  pts = data.get('pts')",
    "  if isinstance(pts, list):",
    "    for row in pts:",
    "      if isinstance(row, list) and len(row) > 2 and isinstance(row[2], (int, float)):",
    "        vals.append(float(row[2]))",
    "      if isinstance(row, list) and row and row[0] in ('arc', 'ellipse') and len(row) > 4 and isinstance(row[4], (int, float)):",
    "        vals.append(float(row[4]))",
    "      if isinstance(row, list) and row and row[0] == 'spline' and len(row) > 3 and isinstance(row[3], dict):",
    "        poles = row[3].get('poles')",
    "        if isinstance(poles, list):",
    "          for p in poles:",
    "            if isinstance(p, (list, tuple)) and len(p) >= 2 and isinstance(p[1], (int, float)):",
    "              vals.append(float(p[1]))",
    "  if isinstance(data.get('VerticalPosition'), (int, float)):",
    "    base = float(data['VerticalPosition'])",
    "    vals.append(base)",
    "    if isinstance(data.get('Height'), (int, float)):",
    "      vals.append(base + float(data['Height']))",
    "    if isinstance(data.get('Radius'), (int, float)):",
    "      vals.append(base + float(data['Radius']))",
    "      vals.append(base - float(data['Radius']))",
    "  sub = data.get('sub_sketches')",
    "  if isinstance(sub, dict):",
    "    for child in sub.values():",
    "      y = _component_top_y(child)",
    "      if y is not None:",
    "        vals.append(y)",
    "  elif isinstance(sub, list):",
    "    for child in sub:",
    "      y = _component_top_y(child)",
    "      if y is not None:",
    "        vals.append(y)",
    "  return max(vals) if vals else None",
    "",
    "def _clone_named_components(component_store, names, owner_name):",
    "  present = {}",
    "  missing = []",
    "  for name in names:",
    "    if name in component_store:",
    "      present[name] = copy.deepcopy(component_store[name])",
    "    else:",
    "      missing.append(name)",
    "  if missing:",
    "    print('[salome config warning] %s references missing components: %s' % (owner_name, ', '.join(sorted(missing))))",
    "  return present",
    "",
    "def _materialize_sketch_buckets(component_store, materials, electrodes_cfg):",
    "  used_names = set()",
    "  xenon_sketches = {}",
    "  electrode_sketches = {}",
    "  supported_materials = {'GXe', 'LXe', 'PTFE'}",
    "  for material_name, entry in materials.items():",
    "    if material_name not in supported_materials:",
    "      raise ValueError('Unsupported material in cleanup config: %s' % material_name)",
    "    components = list((entry or {}).get('components', []) or [])",
    "    if (entry or {}).get('default') or material_name == 'PTFE' or not components:",
    "      continue",
    "    if len(components) == 1 and components[0] == material_name:",
    "      if components[0] not in component_store:",
    "        print('[salome config warning] material %s references missing component: %s' % (material_name, components[0]))",
    "        continue",
    "      xenon_sketches[material_name] = copy.deepcopy(component_store[components[0]])",
    "      used_names.add(components[0])",
    "    else:",
    "      sub_sketches = _clone_named_components(component_store, components, 'material %s' % material_name)",
    "      if not sub_sketches:",
    "        print('[salome config warning] material %s has no available components; skipping' % material_name)",
    "        continue",
    "      xenon_sketches[material_name] = {'hull': True, 'sub_sketches': sub_sketches}",
    "      used_names.update(sub_sketches.keys())",
    "  for electrode_name, entry in electrodes_cfg.items():",
    "    components = list((entry or {}).get('components', []) or [])",
    "    if not components:",
    "      continue",
    "    sub_sketches = _clone_named_components(component_store, components, 'electrode %s' % electrode_name)",
    "    if not sub_sketches:",
    "      print('[salome config warning] electrode %s has no available components; skipping' % electrode_name)",
    "      continue",
    "    electrode_sketches[electrode_name] = {'hull': True, 'sub_sketches': sub_sketches}",
    "    used_names.update(sub_sketches.keys())",
    "  ptfe_sketches = {}",
    "  for name, data in component_store.items():",
    "    if name in used_names:",
    "      continue",
    "    ptfe_sketches[name] = copy.deepcopy(data)",
    "  return ptfe_sketches, electrode_sketches, xenon_sketches",
    "",
    "def build_sketch_dicts(shrinkage_factor, liquid_level=liquid_level):",
    "  component_store = copy.deepcopy(_COMPONENT_STORE)",
    "  ref_name = _SLICE_ALIGNMENT.get('reference_component')",
    "  if ref_name:",
    "    if ref_name not in component_store:",
    "      raise KeyError('slice alignment reference component missing: %s' % ref_name)",
    "    top_y = _component_top_y(component_store[ref_name])",
    "    if top_y is None:",
    "      raise ValueError('could not determine top y for reference component: %s' % ref_name)",
    "    target = _SLICE_ALIGNMENT.get('reference_top_y', '-liquid_level')",
    "    if target != '-liquid_level':",
    "      raise ValueError('unsupported slice alignment target: %s' % target)",
    "    dy = (-float(liquid_level)) - float(top_y)",
    "    for name in list(component_store.keys()):",
    "      if name in _MANUAL_COMPONENT_NAMES:",
    "        continue",
    "      component_store[name] = _translate_component(component_store[name], 0.0, dy)",
    "  ptfe_sketches, electrode_sketches, xenon_sketches = _materialize_sketch_buckets(component_store, _MATERIALS, _ELECTRODES)",
    "  electrode_sketches['shrinkage_factor'] = shrinkage_factor",
    "  manual_mapping = copy.deepcopy(_MANUAL_MAPPING)",
    "  return ptfe_sketches, electrode_sketches, xenon_sketches, manual_mapping",
    "",
  ]
  with open(module_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
  return module_path

def run_cleanup_for_slice(cfg, in_dir, out_dir):
  slice_name = get_slice_name_from_dir(in_dir)
  slice_input_scale = float(cfg.get("slice_input_scale", 1.0))
  shape_match_tol = float(cfg.get("shape_match_tol", 1e-6))
  periodicity_tol = float(cfg.get("periodicity_tol", 1e-4))
  grouping_mode = str(cfg.get("grouping_mode", "periodic"))
  slice_alignment = copy.deepcopy(cfg.get("slice_alignment") or {})
  if slice_alignment and not rule_applies_to_slice(slice_alignment, slice_name):
    slice_alignment = {}
  slice_offset_y = 0.0 if slice_alignment.get("reference_component") else float(cfg.get("slice_offset_y", 0.0))
  manual_offset_y = float(cfg.get("manual_offset_y", 0.0))

  all_files, file_meta = load_and_prepare_jsons(
    in_dir,
    slice_input_scale=slice_input_scale,
    slice_offset_y=slice_offset_y,
  )
  all_files, file_meta, ignored = apply_ignore_components(
    all_files,
    file_meta,
    cfg.get("ignore_components", []),
    slice_name,
  )
  all_files, file_meta, applied_replacements = apply_replacements(
    all_files,
    file_meta,
    cfg.get("replacements", []),
    slice_name,
    in_dir,
    slice_input_scale,
    slice_offset_y,
  )

  used_candidates = set()
  out_groupings = []
  no_group_targets = []
  preserved_unassigned = []
  names = list(all_files.keys())

  if grouping_mode == "periodic":
    # --- apply periodicity overrides EARLY (so they don't participate in grouping) ---
    # assumes this returns (all_files, file_meta, forced_groupings, forced_used_names)
    all_files, file_meta, forced_groupings, forced_used = apply_override_periodicity_pre(
        all_files=all_files,
        file_meta=file_meta,
        override_rules=cfg.get("override_periodicity", []),
        tol=periodicity_tol,
    )

    used_candidates = set(forced_used)
    out_groupings = [
      {"target": target, "groupings": groupings}
      for target, groupings in forced_groupings.items()
    ]
    no_group_targets = []

    names = list(all_files.keys())
    emitted_group_signatures = set()

    for target in names:
      if target in used_candidates:
          continue

      soft_candidates = get_candidates_pre(file_meta, target)
      candidates = filter_candidates(
        all_files,
        soft_candidates,
        file_meta,
        target,
        tol=shape_match_tol,
      )

      # include the target itself
      candidates = list(dict.fromkeys(candidates + [target]))

      groupings = find_periodic_groups(
          file_meta,
          candidates,
          tol_pos=periodicity_tol,
          tol_step=periodicity_tol,
          min_len=3,
      )

      found_any = bool(groupings["x"]) or bool(groupings["y"])
      if not found_any:
        no_group_targets.append(target)
        continue

      target_emitted = False
      target_in_any_group = False
      for axis in ("x", "y"):
        for g in groupings[axis]:
          member_names = {m["name"] for m in g.get("members", [])}
          if target in member_names:
            target_in_any_group = True
          canonical_target = _canonical_member_for_group(g, axis)
          if canonical_target != target:
            continue
          signature = _group_signature(axis, g)
          if signature in emitted_group_signatures:
            continue
          subgroupings = {"x": [], "y": []}
          subgroupings[axis] = [g]
          out_groupings.append({
            "target": target,
            "groupings": subgroupings,
          })
          emitted_group_signatures.add(signature)
          for m in g["members"]:
            used_candidates.add(m["name"])
          target_emitted = True

      if target_emitted:
        continue
      if not target_in_any_group:
        no_group_targets.append(target)

    accounted_for = set(no_group_targets) | set(used_candidates)
    for target in names:
      if target in accounted_for:
        continue
      no_group_targets.append(target)
      preserved_unassigned.append(target)
  else:
    no_group_targets = list(names)

  os.makedirs(out_dir, exist_ok=True)
  written_paths = set()

  # write unchanged originals where no grouping was found
  for target in no_group_targets:
    dst = os.path.join(out_dir, f"{target}.json")
    with open(dst, "w", encoding="utf-8") as f:
      json.dump(all_files[target], f, ensure_ascii=False, indent=2)
    written_paths.add(os.path.abspath(dst))

  # write new files where groupings were found (and apply replace_these + override_periodicity)
  replace_rules = cfg.get("replace_these", [])
  override_rules = cfg.get("override_periodicity", [])

  for entry in out_groupings:
    target = entry["target"]
    groupings = entry["groupings"]
    out_path = make_and_write_grouped(
        target=target,
        groupings=groupings,
        in_dir=in_dir,
        out_dir=out_dir,
        tol=periodicity_tol,
        replace_rules=replace_rules,
        override_rules=override_rules,
        slice_input_scale=slice_input_scale,
        slice_offset_y=slice_offset_y,
    )
    written_paths.add(os.path.abspath(out_path))

  # manual additions copied into out_dir regardless of grouping
  manual_additions_cfg = copy.deepcopy(cfg.get("manual_additions", {}))
  manual_additions_cfg["_default_offset_y"] = manual_offset_y
  manual_written = apply_manual_additions(manual_additions_cfg, out_dir=out_dir)
  for p in manual_written:
    written_paths.add(os.path.abspath(p))

  python_path = None
  if written_paths:
    python_path = write_salome_python_config(
      cfg=cfg,
      slice_name=slice_name,
      written_paths=written_paths,
      out_dir=out_dir,
      manual_written_paths=manual_written,
    )

  if cfg.get("prune_stale_outputs", True):
    removed = 0
    for fname in os.listdir(out_dir):
      if not fname.lower().endswith(".json"):
        continue
      fpath = os.path.abspath(os.path.join(out_dir, fname))
      if fpath in written_paths:
        continue
      os.remove(fpath)
      removed += 1
    print(
      "[cleanup_jsons] wrote=%d grouped=%d copied=%d manual=%d ignored=%d replacements=%d preserved_unassigned=%d pruned=%d slice=%s out_dir=%s python=%s slice_scale=%g slice_offset_y=%g manual_offset_y=%g grouping_mode=%s shape_match_tol=%g periodicity_tol=%g"
      % (
        len(written_paths),
        len(out_groupings),
        len(no_group_targets),
        len(manual_written),
        len(ignored),
        len(applied_replacements),
        len(preserved_unassigned),
        removed,
        slice_name,
        out_dir,
        python_path,
        slice_input_scale,
        slice_offset_y,
        manual_offset_y,
        grouping_mode,
        shape_match_tol,
        periodicity_tol,
      )
    )
  else:
    print(
      "[cleanup_jsons] wrote=%d grouped=%d copied=%d manual=%d ignored=%d replacements=%d preserved_unassigned=%d prune_stale_outputs=false slice=%s out_dir=%s python=%s slice_scale=%g slice_offset_y=%g manual_offset_y=%g grouping_mode=%s shape_match_tol=%g periodicity_tol=%g"
      % (
        len(written_paths),
        len(out_groupings),
        len(no_group_targets),
        len(manual_written),
        len(ignored),
        len(applied_replacements),
        len(preserved_unassigned),
        slice_name,
        out_dir,
        python_path,
        slice_input_scale,
        slice_offset_y,
        manual_offset_y,
        grouping_mode,
        shape_match_tol,
        periodicity_tol,
      )
    )
  return {
    "slice_name": slice_name,
    "out_dir": out_dir,
    "grouped": len(out_groupings),
    "copied": len(no_group_targets),
    "manual": len(manual_written),
    "ignored": len(ignored),
    "replacements": len(applied_replacements),
    "preserved_unassigned": len(preserved_unassigned),
    "written": len(written_paths),
  }


def main():
  cfg_path = "cleanup_config.yaml"
  if len(sys.argv) > 1:
    cfg_path = sys.argv[1]
  cfg = load_yaml_config(cfg_path)  # assumes exists

  jobs = resolve_slice_jobs(cfg["paths"])
  summaries = []
  for in_dir, out_dir in jobs:
    summaries.append(run_cleanup_for_slice(cfg, in_dir, out_dir))

  if len(summaries) > 1:
    print(
      "[cleanup_jsons_batch] slices=%d written=%d grouped=%d copied=%d manual=%d ignored=%d replacements=%d preserved_unassigned=%d"
      % (
        len(summaries),
        sum(s["written"] for s in summaries),
        sum(s["grouped"] for s in summaries),
        sum(s["copied"] for s in summaries),
        sum(s["manual"] for s in summaries),
        sum(s["ignored"] for s in summaries),
        sum(s["replacements"] for s in summaries),
        sum(s["preserved_unassigned"] for s in summaries),
      )
    )


if __name__ == "__main__":
  main()
