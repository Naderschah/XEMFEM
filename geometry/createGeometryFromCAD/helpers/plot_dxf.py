angle = '360'

DXF_GLOB = "/work/geometry/createGeometryFromCAD/DXF_slices_parts/slice_{}.00deg/XENT_TPC_B_Warm*.dxf".format(angle)
OUT_PNG = "PMT_{}.png".format(angle)
# pip install ezdxf matplotlib

import os
import glob
import ezdxf
import matplotlib.pyplot as plt

from ezdxf.recover import readfile as recover_readfile
from ezdxf.addons.drawing import RenderContext, Frontend
from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
from ezdxf.addons.drawing.config import Configuration, ColorPolicy, BackgroundPolicy
from ezdxf.bbox import extents


cfg = Configuration(
    color_policy=ColorPolicy.BLACK,
    background_policy=BackgroundPolicy.WHITE,
)

def spline_is_valid(e) -> bool:
    """Conservative validity checks to avoid ezdxf Basis() 'invalid order'."""
    try:
        # order = degree + 1
        degree = int(getattr(e.dxf, "degree", 3))
        order = degree + 1
        if order < 2:
            return False

        # control points
        cpc = int(getattr(e.dxf, "n_control_points", 0))
        # Some DXFs may not set n_control_points reliably; fall back to length
        if cpc <= 0:
            try:
                cpc = len(e.control_points)  # type: ignore[attr-defined]
            except Exception:
                cpc = 0
        if cpc <= 0:
            return False

        # order must not exceed count
        if order > cpc:
            return False

        # knot count should be count + order (for clamped B-spline; ezdxf expects this)
        try:
            knots = list(e.knots)  # type: ignore[attr-defined]
            if len(knots) != cpc + order:
                return False
        except Exception:
            # If knots are inaccessible, let it pass (recover might have fixed it)
            pass

        return True
    except Exception:
        return False

def filter_bad_entities(e) -> bool:
    # Return True to draw entity, False to skip
    if e.dxftype() == "SPLINE":
        return spline_is_valid(e)
    return True

dxf_files = sorted(glob.glob(DXF_GLOB))
if not dxf_files:
    raise FileNotFoundError(f"No DXF files matched: {DXF_GLOB}")

fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1])
ax.set_aspect("equal", adjustable="box")
ax.axis("off")

backend = MatplotlibBackend(ax)

xmin = ymin = float("inf")
xmax = ymax = float("-inf")

ok = 0
skipped = 0

for i, path in enumerate(dxf_files, 1):
    try:
        # Use recover to tolerate/repair malformed DXFs
        doc, auditor = recover_readfile(path)
        # If you want to see issues:
        # if auditor.has_errors:
        #     print(f"[{i}] AUDIT errors in {path}: {len(auditor.errors)}")

        msp = doc.modelspace()
        ctx = RenderContext(doc)
        Frontend(ctx, backend, config=cfg).draw_layout(
            msp,
            finalize=False,
            filter_func=filter_bad_entities,  # skip invalid SPLINEs
        )

        # Update global extents (also guarded)
        try:
            box = extents(msp, fast=False)
            xmin = min(xmin, box.extmin.x)
            ymin = min(ymin, box.extmin.y)
            xmax = max(xmax, box.extmax.x)
            ymax = max(ymax, box.extmax.y)
        except Exception:
            # If extents fails for this file, ignore it
            pass

        ok += 1
    except Exception as ex:
        skipped += 1
        print(f"[{i}/{len(dxf_files)}] SKIP {path}  ({type(ex).__name__}: {ex})")

# If extents never got updated (all failed), just save what was drawn
if xmin != float("inf"):
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

out_path = os.path.abspath(OUT_PNG)
fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0)
plt.close(fig)

print(f"Rendered: {ok}, skipped files: {skipped}")
print(out_path)