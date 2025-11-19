import numpy as np
from pathlib import Path

# -------- Files --------
phys_map_path = Path("phys_map.txt")
output_path = Path("config_autogen.yaml")

# -------- Fixed definitions --------
fixed_boundaries = {
    "BC_Bell":         {"type": "dirichlet", "value": 10.0},
    "BC_Gate":         {"type": "dirichlet", "value": 5.0},
    "BC_Anode":        {"type": "dirichlet", "value": 0.0},
    "BC_Cathode":      {"type": "dirichlet", "value": -5.0},
    "BC_TopScreen":    {"type": "dirichlet", "value": 0.0},
    "BC_BottomScreen": {"type": "dirichlet", "value": 0.0},
    "BC_CopperRing":   {"type": "dirichlet", "value": 0.0},
    "BC_PMTs":         {"type": "dirichlet", "value": 0.0},
}

# Lookup for materials: region name → epsilon_r
volume_materials = {
    "LXe": 1.95,
    "GXe": 1.0,
    "PTFE": 2.1,
}

# Resistive Chain
R1 = 1.0 # First 5
R2 = 1.0 # Split 
R3 = 1.0 # Middle series 
RC = 1.0 # Last to Cathode 

# -------- Parse phys_map.txt --------
with phys_map_path.open() as f:
    lines = [l.strip() for l in f if l.strip()]
entries = lines[1:]  # skip number row

materials = {}
boundaries_fixed = {}
rings = []  # (name, tag)

for line in entries:
    dim_s, tag_s, name_s = line.split(maxsplit=2)
    dim = int(dim_s)
    tag = int(tag_s)
    name = name_s.strip('"')

    if dim == 2:
        if name in volume_materials:
            materials[name] = {"attr_id": tag, "epsilon_r": volume_materials[name]}

    elif dim == 1:
        if name in fixed_boundaries:
            boundaries_fixed[name] = {
                "bdr_id": tag,
                "type": fixed_boundaries[name]["type"],
                "value": fixed_boundaries[name]["value"],
            }
        else:
            rings.append((name, tag))

# -------- Assign voltages to rings --------
N = len(rings)
if N > 1:
    voltages = np.linspace(-30, 30, N)
else:
    voltages = np.array([0.0])

boundaries_rings = {}

for (name, tag), V in zip(rings, voltages):
    boundaries_rings[name] = {
        "bdr_id": tag,
        "type": "dirichlet",
        "value": float(V),
    }

# -------- Write to config_autogen.yaml --------
with output_path.open("w") as out:
    out.write("materials:\n")
    for m, v in materials.items():
        out.write(f"  {m}:\n")
        out.write(f"    attr_id: {v['attr_id']}\n")
        out.write(f"    epsilon_r: {v['epsilon_r']}\n")

    out.write("\nboundaries:\n")
    for name, v in boundaries_fixed.items():
        out.write(f"  {name}:\n")
        out.write(f"    bdr_id: {v['bdr_id']}\n")
        out.write(f"    type: {v['type']}\n")
        out.write(f"    value: {v['value']}\n")

    for name, v in sorted(boundaries_rings.items()):
        out.write(f"  {name}:\n")
        out.write(f"    bdr_id: {v['bdr_id']}\n")
        out.write(f"    type: {v['type']}\n")
        out.write(f"    value: {v['value']}\n")

print(f"✔ Wrote {output_path} with {len(boundaries_rings)} ring boundaries.")
