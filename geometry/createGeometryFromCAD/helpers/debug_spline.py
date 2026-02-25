"""
Inspect failing DXF SPLINE entities (ezdxf construction_tool / path conversion failures)
and print detailed diagnostics about what is structurally wrong.

Usage:
  python inspect_bad_splines.py input.dxf
Optional:
  python inspect_bad_splines.py input.dxf --write-bad bad_splines_only.dxf --json report.json

Notes:
- In ezdxf: order = degree + 1
- Common failure causes:
  * order < 2 (degree < 1)
  * order > number of control points
  * knot_count != control_point_count + order
  * knots not nondecreasing
  * weights count mismatch
"""

import argparse
import json
import math
import sys
from collections import Counter

import ezdxf
from ezdxf.path import make_path


def is_finite_number(x) -> bool:
    try:
        return x is not None and math.isfinite(float(x))
    except Exception:
        return False


def safe_len(seq):
    try:
        return len(seq)
    except Exception:
        return None


def format_first(seq, n=10):
    try:
        return list(seq)[:n]
    except Exception:
        return None


def monotone_nondecreasing(values):
    # returns (ok, first_bad_index)
    try:
        for i in range(len(values) - 1):
            if values[i] > values[i + 1]:
                return False, i
        return True, None
    except Exception:
        return None, None


def knot_multiplicities(knots, tol=0.0):
    """
    Count knot multiplicities. If tol>0, group nearly-equal knots.
    """
    if knots is None:
        return None
    if tol <= 0:
        c = Counter(knots)
        return dict(sorted(c.items(), key=lambda kv: kv[0]))
    # tolerance grouping
    groups = []
    for k in knots:
        placed = False
        for g in groups:
            if abs(g["value"] - k) <= tol:
                g["count"] += 1
                placed = True
                break
        if not placed:
            groups.append({"value": float(k), "count": 1})
    groups.sort(key=lambda g: g["value"])
    return groups


def summarize_spline(s):
    """
    Extract key SPLINE fields as best as possible without failing.
    """
    degree = getattr(s.dxf, "degree", None)
    order = (degree + 1) if isinstance(degree, int) else None

    # Control points / fit points / knots / weights:
    cp = None
    fp = None
    knots = None
    weights = None

    try:
        cp = list(s.control_points)
    except Exception:
        pass

    try:
        fp = list(s.fit_points)
    except Exception:
        pass

    # knots access differs across versions; try multiple:
    try:
        knots = list(s.knots)
    except Exception:
        try:
            knots = list(s.get_knot_values())  # not always present
        except Exception:
            knots = None

    try:
        weights = list(s.weights)
    except Exception:
        weights = None

    return {
        "handle": getattr(s.dxf, "handle", None),
        "layer": getattr(s.dxf, "layer", None),
        "degree": degree,
        "order": order,
        "flags": getattr(s.dxf, "flags", None),
        "control_points_count": safe_len(cp),
        "fit_points_count": safe_len(fp),
        "knot_count": safe_len(knots),
        "weights_count": safe_len(weights),
        "control_points_first": format_first(cp, 5),
        "fit_points_first": format_first(fp, 5),
        "knots_first": format_first(knots, 12),
        "weights_first": format_first(weights, 12),
        "knots": knots,
        "weights": weights,
    }


def validate_spline_struct(info):
    """
    Run common NURBS/B-spline consistency checks and return a list of issues.
    """
    issues = []
    degree = info["degree"]
    order = info["order"]
    n_ctrl = info["control_points_count"]
    n_knots = info["knot_count"]
    n_w = info["weights_count"]
    knots = info["knots"]
    weights = info["weights"]

    # Basic degree/order checks
    if not isinstance(degree, int):
        issues.append("degree is missing or not an int")
    else:
        if degree < 1:
            issues.append(f"degree={degree} < 1 -> order={order} invalid (order must be >= 2)")
    if isinstance(order, int) and order < 2:
        issues.append(f"order={order} < 2")

    # Control points vs order
    if isinstance(order, int) and isinstance(n_ctrl, int):
        if n_ctrl <= 0:
            issues.append("no control points")
        if order > n_ctrl:
            issues.append(f"order={order} > control_points_count={n_ctrl} (too few control points for degree)")

    # Knot vector checks
    if knots is None:
        issues.append("knot vector missing/unreadable")
    else:
        if isinstance(order, int) and isinstance(n_ctrl, int) and isinstance(n_knots, int):
            expected = n_ctrl + order
            if n_knots != expected:
                issues.append(f"knot_count={n_knots} != control_points_count+order={expected} (common export corruption)")

        ok, idx = monotone_nondecreasing(knots)
        if ok is False:
            issues.append(f"knot vector is not nondecreasing (first decrease at index {idx}: {knots[idx]} > {knots[idx+1]})")
        elif ok is None:
            issues.append("knot vector monotonicity check failed")

        # finiteness
        bad_knots = [i for i, k in enumerate(knots) if not is_finite_number(k)]
        if bad_knots:
            issues.append(f"knot vector contains non-finite values at indices {bad_knots[:10]}")

        # typical clamped endpoints: multiplicity at ends should be == order
        # (not mandatory for all splines, but very common in CAD)
        if isinstance(order, int) and len(knots) > 0:
            start = knots[0]
            end = knots[-1]
            start_mult = sum(1 for k in knots if k == start)
            end_mult = sum(1 for k in knots if k == end)
            if start_mult < min(order, len(knots)):
                issues.append(f"start knot multiplicity {start_mult} < order {order} (may be unclamped or invalid)")
            if end_mult < min(order, len(knots)):
                issues.append(f"end knot multiplicity {end_mult} < order {order} (may be unclamped or invalid)")

    # Weights checks (NURBS/rational). If weights absent, that’s usually OK (all weights=1).
    if weights is not None:
        if isinstance(n_ctrl, int) and isinstance(n_w, int) and n_w != n_ctrl:
            issues.append(f"weights_count={n_w} != control_points_count={n_ctrl}")
        bad_w = [i for i, w in enumerate(weights) if (not is_finite_number(w) or float(w) <= 0.0)]
        if bad_w:
            issues.append(f"weights contain non-finite or non-positive values at indices {bad_w[:10]}")

    return issues


def try_failure_reason(s):
    """
    Try the same operations that typically fail in plotting:
      - s.construction_tool()
      - make_path(s)
    Return exception info if any.
    """
    for name, fn in [
        ("construction_tool", lambda: s.construction_tool()),
        ("make_path", lambda: make_path(s)),
    ]:
        try:
            fn()
        except Exception as e:
            return {
                "failed_at": name,
                "exception_type": type(e).__name__,
                "exception_msg": str(e),
            }
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dxf", help="Input DXF file")
    ap.add_argument("--write-bad", help="Write a DXF containing only failing SPLINE entities")
    ap.add_argument("--json", help="Write a JSON report")
    ap.add_argument("--max", type=int, default=0, help="Max number of failing splines to print (0 = no limit)")
    args = ap.parse_args()

    doc = ezdxf.readfile(args.dxf)
    msp = doc.modelspace()

    bad_infos = []
    total = 0
    bad = 0

    for s in msp.query("SPLINE"):
        total += 1

        fail = try_failure_reason(s)
        if not fail:
            continue

        bad += 1
        info = summarize_spline(s)
        info["failure"] = fail
        info["issues"] = validate_spline_struct(info)
        # add multiplicities (use exact equality; add tol version if needed)
        if info["knots"] is not None:
            info["knot_multiplicities"] = knot_multiplicities(info["knots"], tol=0.0)

        bad_infos.append(info)

        if args.max and bad > args.max:
            break

    # Print a readable report
    print(f"Total SPLINE entities: {total}")
    print(f"Failing SPLINE entities: {bad}")
    for i, info in enumerate(bad_infos, start=1):
        print("\n" + "=" * 80)
        print(f"[{i}] handle={info['handle']} layer={info['layer']}")
        print(f"failed_at={info['failure']['failed_at']}: {info['failure']['exception_type']}: {info['failure']['exception_msg']}")
        print(f"degree={info['degree']} order={info['order']} flags={info['flags']}")
        print(f"control_points={info['control_points_count']} fit_points={info['fit_points_count']} knots={info['knot_count']} weights={info['weights_count']}")
        print(f"knots_first={info['knots_first']}")
        print(f"weights_first={info['weights_first']}")
        if info["issues"]:
            print("issues:")
            for msg in info["issues"]:
                print(f"  - {msg}")
        else:
            print("issues: (none detected by basic checks)")
        # optional: show end multiplicities compactly
        km = info.get("knot_multiplicities")
        if isinstance(km, dict) and len(km) <= 20:
            print(f"knot multiplicities (all): {km}")
        elif km is not None:
            print(f"knot multiplicities: {len(km)} unique values (not printed)")

    # Write DXF with only bad splines
    if args.write_bad:
        out = ezdxf.new(doc.dxfversion)
        out_msp = out.modelspace()

        # Make layers used by bad splines so they appear reasonably
        layers = set(i["layer"] for i in bad_infos if i.get("layer"))
        for layer in layers:
            if layer not in out.layers:
                out.layers.new(layer)

        for info in bad_infos:
            h = info["handle"]
            # re-fetch entity by handle from original doc:
            try:
                ent = doc.entitydb.get(h)
                if ent is None:
                    continue
                out_msp.add_entity(ent.copy())
            except Exception:
                # fallback: skip if copy fails
                continue

        out.saveas(args.write_bad)
        print(f"\nWrote failing splines to: {args.write_bad}")

    # Write JSON report
    if args.json:
        # remove potentially huge arrays if you want a smaller file:
        compact = []
        for info in bad_infos:
            c = dict(info)
            # keep full knots/weights if you want deep analysis; comment these out to shrink:
            # c.pop("knots", None)
            # c.pop("weights", None)
            compact.append(c)

        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "input": args.dxf,
                    "total_splines": total,
                    "failing_splines": bad,
                    "items": compact,
                },
                f,
                indent=2,
                ensure_ascii=False,
            )
        print(f"Wrote JSON report to: {args.json}")


if __name__ == "__main__":
    main()