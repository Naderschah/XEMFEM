"""
This code will handle the merging of repeated elements


"""
import os, json
from collections import Counter
from collections import defaultdict
import shutil

def load_and_prepare_jsons(path):
  """
  Takes the path to one slice directory

  returns all objects in a json and a few precomputed things
  """
  all_files = {}
  for i in os.listdir(path):
    all_files[".".join(i.split(".")[:-2])] = json.loads(os.path.join(path, i))
  # Data for cheap comparisons
  file_meta = {
    i:{
      'n_points': len(all_files[i]['pts'])
      'min_x': min(pts_row[1] for pts_row in all_files[i]['pts'])
      'min_y': min(pts_row[2] for pts_row in all_files[i]['pts'])
    }
    for i in all_files.keys()
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

def make_and_write_grouped(target, groupings, in_dir, out_dir):
  """
  Loads ONE original JSON (chosen from an extremum member of the best periodic group),
  then writes a new JSON containing:
    {
      "pts": ...,
      "HorizontalPitch" or "VerticalPitch": <signed displacement>,
      "Number": <repetitions not including original>
    }
  """
  os.makedirs(out_dir, exist_ok=True)

  def load_json(name):
    path = os.path.join(in_dir, f"{name}.json")
    with open(path, "r", encoding="utf-8") as f:
      return json.load(f)

  best_axis = None
  best_group = None
  best_len = 0

  for axis in ("x", "y"):
    for g in groupings.get(axis, []):
      L = len(g.get("members", []))
      if L > best_len:
        best_len = L
        best_axis = axis
        best_group = g

  if best_group is None or best_len < 2:
      raise ValueError(f"No usable periodic grouping for target={target}")

  members = best_group["members"]

  coord_key = "min_x" if best_axis == "x" else "min_y"
  members_sorted = sorted(members, key=lambda m: m[coord_key])
  base = members_sorted[0]
  base_name = base["name"]

  if len(members_sorted) >= 2:
      pitch = members_sorted[1][coord_key] - members_sorted[0][coord_key]
  else:
      raise Exception("Logic Error: Only one group memeber")

  number = len(members_sorted) - 1

  data = load_json(base_name)

  out = {
      "pts": data.get("pts", []),
      ("HorizontalPitch" if best_axis == "x" else "VerticalPitch"): pitch,
      "Number": number,
  }

  out_path = os.path.join(out_dir, f"{target}_grouped.json")
  with open(out_path, "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

import re

def apply_override_periodicity_pre(all_files, file_meta, override_rules, tol=1e-6):
  """
  Pre-pass that:
    1) Finds targets matching override_periodicity regex rules
    2) Resolves missing pitch (if pitch is None) from min-coordinate spacing along the chosen axis
    3) Produces forced_groupings + forced_used_names so these objects do not participate in later grouping
    4) Removes matched objects from all_files + file_meta (preferred so they won't be considered later)

  Returns:
    all_files, file_meta, forced_groupings, forced_used_names
  """
  forced_groupings = {}
  forced_used = set()

  if not override_rules:
    return all_files, file_meta, forced_groupings, forced_used

  # Snapshot names once (we will pop from dicts)
  names = list(all_files.keys())

  for rule in override_rules:
    pattern = rule["match"]
    rx = re.compile(pattern)

    use_source = rule.get("use_source", "min_extremum")
    set_cfg = rule.get("set", {})
    axis = set_cfg.get("axis")  # "x" or "y"
    if axis not in ("x", "y"):
        raise ValueError(f"override_periodicity rule axis must be 'x' or 'y': {rule}")

    coord_key = "min_x" if axis == "x" else "min_y"

    matched = [n for n in names if rx.search(n)]
    if not matched:
      continue

    # Resolve pitch if not provided: minimal positive spacing between sorted min coords
    pitch = set_cfg.get("pitch", None)
    if pitch is None:
      coords = sorted(file_meta[n][coord_key] for n in matched if n in file_meta)
      diffs = [coords[i + 1] - coords[i] for i in range(len(coords) - 1)]
      pos_diffs = [d for d in diffs if d > tol]
      pitch_mag = min(pos_diffs) if pos_diffs else 0.0
      pitch = pitch_mag

    # Enforce direction based on extremum choice unless pitch is explicitly signed how you want.
    # Convention used here:
    #   min_extremum -> +|pitch|
    #   max_extremum -> -|pitch|
    pitch = abs(pitch)
    if use_source == "max_extremum":
      pitch = -pitch

    number = set_cfg.get("number")
    if number is None:
      raise ValueError(f"override_periodicity requires 'number' (repetitions): {rule}")

    # Create forced grouping entries and remove from all_files/file_meta
    for n in matched:
      if n not in all_files:
        continue

      forced_used.add(n)

      # Minimal "groupings" structure so downstream code can call make_and_write_grouped.
      # make_and_write_grouped is expected to apply override_rules again; we also provide
      # resolved values here as a convenience.
      forced_groupings[n] = {
          "x": [],
          "y": [],
          "_pre_override": {
              "match": pattern,
              "use_source": use_source,
              "resolved": {
                  "axis": axis,
                  "pitch": pitch,
                  "number": number,
              },
          },
      }

      # Prefer excluding these from later grouping computation
      all_files.pop(n, None)
      file_meta.pop(n, None)

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
    return

  copy_rules = manual_additions_cfg.get("copy", [])
  if not copy_rules:
    return

  os.makedirs(out_dir, exist_ok=True)

  for rule in copy_rules:
    # destination directory (per rule override)
    dst_dir = rule.get("to_dir", out_dir)

    if "from" in rule:
      src = rule["from"]
      stem = rule.get("to_name") or os.path.splitext(os.path.basename(src))[0]
      dst = os.path.join(dst_dir, f"{stem}.json")
      os.makedirs(dst_dir, exist_ok=True)
      shutil.copy2(src, dst)
      continue

    if "from_glob" in rule:
      pattern = rule["from_glob"]
      matches = sorted(glob.glob(pattern))
      os.makedirs(dst_dir, exist_ok=True)
      for src in matches:
        stem = os.path.splitext(os.path.basename(src))[0]
        dst = os.path.join(dst_dir, f"{stem}.json")
        shutil.copy2(src, dst)
      continue

    raise ValueError(f"manual_additions.copy entry must have 'from' or 'from_glob': {rule}")

def main():
  tol = 1e-6

  cfg_path = "config.yaml"  # wherever you keep it
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

  # write unchanged originals where no grouping was found
  for target in no_group_targets:
    src = os.path.join(in_dir, f"{target}.json")
    dst = os.path.join(out_dir, f"{target}.json")
    shutil.copy2(src, dst)

  # write new files where groupings were found (and apply replace_these + override_periodicity)
  replace_rules = cfg.get("replace_these", [])
  override_rules = cfg.get("override_periodicity", [])

  for target, groupings in out_groupings.items():
    make_and_write_grouped(
        target=target,
        groupings=groupings,
        in_dir=in_dir,
        out_dir=out_dir,
        tol=tol,
        replace_rules=replace_rules,
        override_rules=override_rules,
    )

  # manual additions copied into out_dir regardless of grouping
  apply_manual_additions(cfg.get("manual_additions", {}), out_dir=out_dir)