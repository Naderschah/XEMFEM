"""
This code will handle the merging of repeated elements


"""
import os, json
from collections import Counter
from collections import defaultdict
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
          manual_additions: {}
          prune_stale_outputs: True

    Expected minimal config:
      paths:
        in_dir: "..."
        out_dir: "..."
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

    for k in ("in_dir", "out_dir"):
        if k not in paths or not isinstance(paths[k], str) or not paths[k].strip():
            raise ValueError(f"config.yaml paths.{k} must be a non-empty string")

    # Resolve actual path fields relative to config directory
    paths["in_dir"] = _resolve_path(paths["in_dir"])
    paths["out_dir"] = _resolve_path(paths["out_dir"])

    if "extras_dir" in paths and paths["extras_dir"] is not None:
        if not isinstance(paths["extras_dir"], str):
            raise ValueError("config.yaml paths.extras_dir must be a string (or omitted)")
        paths["extras_dir"] = _resolve_path(paths["extras_dir"])

    # Optional sections with defaults
    cfg.setdefault("override_periodicity", [])
    cfg.setdefault("replace_these", [])
    cfg.setdefault("manual_additions", {})
    cfg.setdefault("prune_stale_outputs", True)

    if not isinstance(cfg["override_periodicity"], list):
        raise ValueError("override_periodicity must be a list")
    if not isinstance(cfg["replace_these"], list):
        raise ValueError("replace_these must be a list")
    if not isinstance(cfg["manual_additions"], dict):
        raise ValueError("manual_additions must be a mapping/object")
    if not isinstance(cfg["prune_stale_outputs"], bool):
        raise ValueError("prune_stale_outputs must be a boolean")

    # Resolve only known path fields inside sections (leave regex 'match' untouched)
    for rule in cfg.get("replace_these", []) or []:
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
                if isinstance(rule.get("from"), str) and rule["from"].strip():
                    rule["from"] = _resolve_path(rule["from"])
                if isinstance(rule.get("from_glob"), str) and rule["from_glob"].strip():
                    rule["from_glob"] = _resolve_path(rule["from_glob"])
                if isinstance(rule.get("to_dir"), str) and rule["to_dir"].strip():
                    rule["to_dir"] = _resolve_path(rule["to_dir"])

    return cfg

def load_and_prepare_jsons(path):
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
      all_files[key] = json.load(f)

  file_meta = {
      name: {
          "n_points": len(all_files[name]["pts"]),
          "min_x": min(row[1] for row in all_files[name]["pts"]),
          "min_y": min(row[2] for row in all_files[name]["pts"]),
      }
      for name in all_files.keys()
  }
  return all_files, file_meta

def get_candidates_pre(file_meta, target):
  """
  Finds all candidates based on n_points
  returned as list of strings
  """
  return [
    i
    for i in file_meta.keys()
    if (file_meta[target]['n_points'] == file_meta[i]['n_points'])
      and (((file_meta[target]['min_x'] - file_meta[i]['min_x']) != 0.0) != 
            ((file_meta[target]['min_y'] - file_meta[i]['min_y']) != 0.0))
  ]

def _quantize(v: float, tol: float) -> int:
    # maps v to an integer grid index; values within ~tol fall into same bucket
    return int(round(v / tol))

def _normalize_pts(pts, min_x: float, min_y: float, tol: float):
    # returns hashable keys for Counter: (label, qx, qy)
    out = []
    for j in pts:
        label = j[0]
        x = j[1] - min_x
        y = j[2] - min_y
        out.append((label, _quantize(x, tol), _quantize(y, tol)))
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
):
  os.makedirs(out_dir, exist_ok=True)

  def load_by_name(name):
    return _load_json_file(os.path.join(in_dir, "%s.json" % name))

  override = _resolve_override_for_target(target, groupings, override_rules or [])

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
      members_sorted = sorted(members, key=lambda m: m.get(coord_key, 0.0))
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
    members_sorted = sorted(members, key=lambda m: m[coord_key])

    if len(members_sorted) < 2:
      raise ValueError("Grouping too small to compute pitch for target=%s" % target)

    pitch = members_sorted[1][coord_key] - members_sorted[0][coord_key]
    number = len(members_sorted) - 1

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

  for rule in copy_rules:
    # destination directory (per rule override)
    dst_dir = rule.get("to_dir", out_dir)

    if "from" in rule:
      src = rule["from"]
      stem = rule.get("to_name") or os.path.splitext(os.path.basename(src))[0]
      dst = os.path.join(dst_dir, f"{stem}.json")
      os.makedirs(dst_dir, exist_ok=True)
      shutil.copy2(src, dst)
      written_paths.append(dst)
      continue

    if "from_glob" in rule:
      pattern = rule["from_glob"]
      matches = sorted(glob.glob(pattern))
      os.makedirs(dst_dir, exist_ok=True)
      for src in matches:
        stem = os.path.splitext(os.path.basename(src))[0]
        dst = os.path.join(dst_dir, f"{stem}.json")
        shutil.copy2(src, dst)
        written_paths.append(dst)
      continue

    raise ValueError(f"manual_additions.copy entry must have 'from' or 'from_glob': {rule}")
  return written_paths

def main():
  tol = 1e-6

  cfg_path = "cleanup_config.yaml"
  if len(sys.argv) > 1:
    cfg_path = sys.argv[1]
  cfg = load_yaml_config(cfg_path)  # assumes exists

  in_dir = cfg["paths"]["in_dir"]
  out_dir = cfg["paths"]["out_dir"]

  all_files, file_meta = load_and_prepare_jsons(in_dir)

  # --- apply periodicity overrides EARLY (so they don't participate in grouping) ---
  # assumes this returns (all_files, file_meta, forced_groupings, forced_used_names)
  all_files, file_meta, forced_groupings, forced_used = apply_override_periodicity_pre(
      all_files=all_files,
      file_meta=file_meta,
      override_rules=cfg.get("override_periodicity", []),
      tol=tol,
  )

  used_candidates = set(forced_used)
  out_groupings = dict(forced_groupings)   # targets -> groupings/override spec
  no_group_targets = []

  names = list(all_files.keys())

  for target in names:
    if target in used_candidates:
        continue

    soft_candidates = get_candidates_pre(file_meta, target)
    candidates = filter_candidates(all_files, soft_candidates, file_meta, target, tol=tol)

    # include the target itself
    candidates = list(dict.fromkeys(candidates + [target]))

    groupings = find_periodic_groups(
        file_meta,
        candidates,
        tol_pos=tol,
        tol_step=tol,
        min_len=3,
    )

    found_any = bool(groupings["x"]) or bool(groupings["y"])
    if not found_any:
      no_group_targets.append(target)
      continue

    out_groupings[target] = groupings

    # mark every object appearing in any grouping as used
    for axis in ("x", "y"):
      for g in groupings[axis]:
        for m in g["members"]:
          used_candidates.add(m["name"])

  os.makedirs(out_dir, exist_ok=True)
  written_paths = set()

  # write unchanged originals where no grouping was found
  for target in no_group_targets:
    src = os.path.join(in_dir, f"{target}.json")
    dst = os.path.join(out_dir, f"{target}.json")
    shutil.copy2(src, dst)
    written_paths.add(os.path.abspath(dst))

  # write new files where groupings were found (and apply replace_these + override_periodicity)
  replace_rules = cfg.get("replace_these", [])
  override_rules = cfg.get("override_periodicity", [])

  for target, groupings in out_groupings.items():
    out_path = make_and_write_grouped(
        target=target,
        groupings=groupings,
        in_dir=in_dir,
        out_dir=out_dir,
        tol=tol,
        replace_rules=replace_rules,
        override_rules=override_rules,
    )
    written_paths.add(os.path.abspath(out_path))

  # manual additions copied into out_dir regardless of grouping
  manual_written = apply_manual_additions(cfg.get("manual_additions", {}), out_dir=out_dir)
  for p in manual_written:
    written_paths.add(os.path.abspath(p))

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
      "[cleanup_jsons] wrote=%d grouped=%d copied=%d manual=%d pruned=%d out_dir=%s"
      % (
        len(written_paths),
        len(out_groupings),
        len(no_group_targets),
        len(manual_written),
        removed,
        out_dir,
      )
    )
  else:
    print(
      "[cleanup_jsons] wrote=%d grouped=%d copied=%d manual=%d prune_stale_outputs=false out_dir=%s"
      % (
        len(written_paths),
        len(out_groupings),
        len(no_group_targets),
        len(manual_written),
        out_dir,
      )
    )


if __name__ == "__main__":
  main()
