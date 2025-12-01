import numpy as np
from pathlib import Path
import re

"""
Top Field Cage has maximum index 
Same with guard
"""

# ============================================================
# 0. USER CONFIGURATION
# ============================================================

phys_map_path = Path("phys_map.txt")
output_path   = Path("config_autogen.yaml")

# Volumes in phys_map → epsilon_r
volume_materials = {
    "LXeGroup":   1.95,
    "GXeGroup":   1.0,
    "PTFE_Group": 2.1,
}


grounded = {"type": "dirichlet", "value": 0}
# Fixed boundary conditions (explicit Dirichlet values)
# Names MUST match phys_map.txt
fixed_boundaries = {
    # SR2 config https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenonnt:analysis:available_data_for_sr2#dates
    "BC_TopScreen":    {"type": "dirichlet", "value": -100},
    "BC_Anode":        {"type": "dirichlet", "value": 4850},
    "BC_Gate":         {"type": "dirichlet", "value": +300},
    "FieldCage_71":    {"type": "dirichlet", "value": +650.},
    "BC_Cathode":      {"type": "dirichlet", "value": -2750},
    "BC_BottomScreen": {"type": "dirichlet", "value": -2750},

    "BC_PMT":          {"type": "dirichlet", "value": -1300},
    "BC_Bell":         grounded,
    "BC_CopperRing":   grounded,
    "BC_Cryostat":     grounded
}
fixed_boundaries = {
    # SR0 config From Francescos thesis
    "BC_TopScreen":    {"type": "dirichlet", "value": -900},
    "BC_Anode":        {"type": "dirichlet", "value": 4900},
    "BC_Gate":         {"type": "dirichlet", "value": 300},
    "FieldCage_71":    {"type": "dirichlet", "value": 650.},
    "BC_Cathode":      {"type": "dirichlet", "value": -2750},
    "BC_BottomScreen": {"type": "dirichlet", "value": -1300},

    "BC_PMT":          {"type": "dirichlet", "value": -1300},
    "BC_Bell":         grounded,
    "BC_CopperRing":   grounded,
    "BC_Cryostat":     grounded
}
fixed_boundaries = {
    # Design config From Francescos thesis
    "BC_TopScreen":    {"type": "dirichlet", "value": -1500},
    "BC_Anode":        {"type": "dirichlet", "value": 6500},
    "BC_Gate":         {"type": "dirichlet", "value": -1000},
    "FieldCage_71":    {"type": "dirichlet", "value": -950.},
    "BC_Cathode":      {"type": "dirichlet", "value": -30000},
    "BC_BottomScreen": {"type": "dirichlet", "value": -1500},

    "BC_PMT":          {"type": "dirichlet", "value": -1300},
    "BC_Bell":         grounded,
    "BC_CopperRing":   grounded,
    "BC_Cryostat":     grounded
}

# Resistor values (only ratios matter for voltages)
# These will be written into fieldcage_network.resistors
R1 = 1.25e9   # between first 5 FieldCage_ and also last 3 after merge
R2 = 2.5e9  # split / merge resistors
R3 = 5e9    # middle region resistors in both branches
R_C = 7.e9  # from last FieldCage_ to Cathode


# ============================================================
# 1. PARSE phys_map.txt
# ============================================================

_idx_re = re.compile(r".*_(\d+)$")

def extract_index(name: str, prefix: str):
    """Extract integer index from names like 'FieldCage_12' or 'FieldCageGuard_3'."""
    if not name.startswith(prefix):
        return None
    m = _idx_re.match(name)
    if not m:
        return None
    return int(m.group(1))


def parse_phys_map(path: Path):
    with path.open() as f:
        lines = [l.strip() for l in f if l.strip()]

    n_total = int(lines[0])
    entries = lines[1:]

    regions = {}    # dim==2
    boundaries = {} # dim==1

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
            raise ValueError(f"Unexpected dim {dim} in phys_map line: {line}")

    print(f"[parse] total entries in phys_map: {n_total}")
    print(f"[parse] boundaries: {len(boundaries)}, regions: {len(regions)}")

    return regions, boundaries


regions_raw, boundaries_raw = parse_phys_map(phys_map_path)


# ============================================================
# 2. COMPUTE VOLTAGES FOR FieldCage_* AND FieldCageGuard_*
# ============================================================

def solve_fieldcage_network(fieldcage_names, guard_names, V_top, V_cathode):
    """
    Build and solve the resistor network:

      - First 5 FieldCage_ in series with R1
      - At FieldCage_5, split via R2:
          * to next FieldCage_ (main chain)
          * to FieldCageGuard_1 (guard chain)
      - Both chains continue in series with R3
      - At the bottom, guard chain re-merges to main chain via R2
      - Last 3 FieldCage_ after merge use R1
      - Last FieldCage_ connects to Cathode via R_C

    V_top:     potential at FieldCage_1
    V_cathode: potential at Cathode (from BC_Cathode)

    Returns dict: name -> voltage for all FieldCage_* and FieldCageGuard_*.
    """

    if not fieldcage_names:
        return {}

    # Sort by index
    fieldcage_names = sorted(fieldcage_names, key=lambda n: extract_index(n, "FieldCage_"), reverse=True)
    guard_names     = sorted(guard_names,     key=lambda n: extract_index(n, "FieldCageGuard_"), reverse=True)

    fc     = fieldcage_names
    guards = guard_names
    n_fc   = len(fc)
    n_guard = len(guards)

    # Fixed potentials in the network
    fixed_V = {
        fc[0]:      V_top,     # first field cage
        "Cathode":  V_cathode,
    }

    # Unknown nodes: all other FieldCage_* and all FieldCageGuard_*
    unknown_nodes = list(fc[1:]) + list(guards)
    N = len(unknown_nodes)

    if N == 0:
        return {fc[0]: V_top}

    node_index = {name: i for i, name in enumerate(unknown_nodes)}
    G = np.zeros((N, N), dtype=float)
    b = np.zeros(N, dtype=float)

    edges = []  # for consistency checks

    def add_resistor(a, bnode, R):
        g = 1.0 / R
        edges.append((a, bnode, R))

        a_unknown = a in node_index
        b_unknown = bnode in node_index
        a_fixed   = a in fixed_V
        b_fixed   = bnode in fixed_V

        if a_unknown and b_unknown:
            ia = node_index[a]
            ib = node_index[bnode]
            G[ia, ia] += g
            G[ib, ib] += g
            G[ia, ib] -= g
            G[ib, ia] -= g
        elif a_unknown and b_fixed:
            ia = node_index[a]
            G[ia, ia] += g
            b[ia] += g * fixed_V[bnode]
        elif b_unknown and a_fixed:
            ib = node_index[bnode]
            G[ib, ib] += g
            b[ib] += g * fixed_V[a]
        else:
            # both fixed or both unknown-without-equation: nothing to add to G,b
            pass

    # If there are not enough rings or no guards, fall back to simple series chain
    if n_fc < 8 or n_guard == 0:
        for i in range(n_fc - 1):
            add_resistor(fc[i], fc[i + 1], R1)
        add_resistor(fc[-1], "Cathode", R_C)
    else:
        # First 5 FieldCage_ with R1: FC1-2-3-4-5
        for i in range(0, 4):
            add_resistor(fc[i], fc[i + 1], R1)

        # Split at FC_5 (index 4) via R2:
        #   FC_5 --R2--> FC_6
        #   FC_5 --R2--> Guard_1
        if n_fc > 5:
            add_resistor(fc[4], fc[5], R2)
        if n_guard > 0:
            add_resistor(fc[4], guards[0], R2)

        # Middle sections with R3
        merge_index = n_fc - 3  # 4th-from-last fieldcage is merge node

        # Main branch: FC_6 ... FC_merge
        for i in range(5, merge_index):
            add_resistor(fc[i], fc[i + 1], R3)

        # Guard branch: Guard_1 ... Guard_last
        for j in range(0, n_guard - 1):
            add_resistor(guards[j], guards[j + 1], R3)

        # Merge: Guard_last --R2--> FC_merge
        if n_guard > 0:
            add_resistor(guards[-1], fc[merge_index], R2)

        # Last FieldCage_ (including merge node) up to last ring with R1
        for i in range(merge_index, n_fc - 1):
            add_resistor(fc[i], fc[i + 1], R1)

        # Last FieldCage_ to Cathode via R_C
        add_resistor(fc[-1], "Cathode", R_C)

    # Solve the linear system G * V_unknown = b
    V_unknown = np.linalg.solve(G, b)

    voltages = {}
    voltages.update(fixed_V)
    for name, idx in node_index.items():
        voltages[name] = float(V_unknown[idx])

    # -------- Self-consistency checks --------
    V1 = voltages[fc[0]]
    Vc = fixed_boundaries["BC_Cathode"]["value"]
    drop_direct = V1 - Vc

    # Path 1: straight FieldCage chain FC1 → ... → FC_last → Cathode
    drop_fc_chain = 0.0
    for a, bnode in zip(fc, fc[1:]):
        drop_fc_chain += voltages[a] - voltages[bnode]
    drop_fc_chain += voltages[fc[-1]] - Vc

    print(f"[check] ΔV(FC1→Cathode) (direct) = {drop_direct:.6g}")
    print(f"[check] ΔV via FieldCage chain    = {drop_fc_chain:.6g} "
          f"(Δ = {drop_fc_chain - drop_direct:+.3e})")

    # Path 2: use guard branch (if present)
    if n_guard > 0 and n_fc >= 8:
        merge_index = n_fc - 4
        merge_name = fc[merge_index]
        path_guard = (
            [fc[0], fc[1], fc[2], fc[3], fc[4]] +
            guards +
            [merge_name] +
            fc[merge_index + 1:] +
            ["Cathode"]
        )
        drop_guard = 0.0
        for a, bnode in zip(path_guard, path_guard[1:]):
            Va = voltages[a] if a in voltages else Vc
            Vb = voltages[bnode] if bnode in voltages else Vc
            drop_guard += Va - Vb
        print(f"[check] ΔV via guard branch      = {drop_guard:.6g} "
              f"(Δ = {drop_guard - drop_direct:+.3e})")

    return voltages


# Collect FieldCage_*/FieldCageGuard_* names from phys_map
fieldcage_names = [name for name in boundaries_raw if name.startswith("FieldCage_")]
guard_names     = [name for name in boundaries_raw if name.startswith("FieldCageGuard_")]

print(f"[network] FieldCage_*: {len(fieldcage_names)}, FieldCageGuard_*: {len(guard_names)}")

# Determine which FieldCage_* is the first in the chain (highest index)
top_fc_name = max(fieldcage_names, key=lambda n: extract_index(n, "FieldCage_"))

fc_voltages = solve_fieldcage_network(
    fieldcage_names,
    guard_names,
    V_top=fixed_boundaries[top_fc_name]["value"],
    V_cathode=fixed_boundaries["BC_Cathode"]["value"],
)


# ============================================================
# 3. BUILD MATERIALS + BOUNDARIES DATA
# ============================================================

# 3a. Materials: find any region without epsilon_r
materials = {}
missing_material_eps = []
for name, meta in regions_raw.items():
    if name in volume_materials:
        materials[name] = {
            "attr_id": meta["attr_id"],
            "epsilon_r": volume_materials[name],
        }
    else:
        missing_material_eps.append(name)

if missing_material_eps:
    print("[warning] No epsilon_r specified for materials:")
    for name in sorted(missing_material_eps):
        print(f"  - {name}")

# 3b. Boundaries: fixed + rings
boundaries_fixed_out = {}
for name, bc in fixed_boundaries.items():
    if name not in boundaries_raw:
        print(f"[warning] Fixed boundary '{name}' not found in phys_map; skipping.")
        continue
    tag = boundaries_raw[name]["bdr_id"]
    boundaries_fixed_out[name] = {
        "bdr_id": tag,
        "type": bc["type"],
        "value": bc["value"],
    }

boundaries_rings = {}
for name in fieldcage_names:
    tag = boundaries_raw[name]["bdr_id"]
    boundaries_rings[name] = {
        "bdr_id": tag,
        "type": "dirichlet",
        "value": float(fc_voltages[name]),
    }
for name in guard_names:
    tag = boundaries_raw[name]["bdr_id"]
    boundaries_rings[name] = {
        "bdr_id": tag,
        "type": "dirichlet",
        "value": float(fc_voltages[name]),
    }

assigned_boundary_names = set(boundaries_fixed_out) | set(boundaries_rings)
unassigned_boundaries = sorted(
    name for name in boundaries_raw if name not in assigned_boundary_names
)

if unassigned_boundaries:
    print("[warning] The following boundaries have no assigned voltage:")
    for name in unassigned_boundaries:
        tag = boundaries_raw[name]["bdr_id"]
        print(f"  - {name} (bdr_id={tag})")
else:
    print("[info] All boundaries in phys_map have an assigned voltage.")


# ============================================================
# 4. BUILD FIELDCAGE NETWORK DESCRIPTION FOR CONFIG
# ============================================================

# Sort fieldcage / guard names by index
fieldcage_names_sorted = sorted(fieldcage_names, key=lambda n: extract_index(n, "FieldCage_"), reverse=True)
guard_names_sorted     = sorted(guard_names,     key=lambda n: extract_index(n, "FieldCageGuard_"), reverse=True)

n_fc   = len(fieldcage_names_sorted)
n_guard = len(guard_names_sorted)

nodes_desc = []
edges_desc = []

if n_fc > 0:
    top_fc_name = max(fieldcage_names, key=lambda n: extract_index(n, "FieldCage_"))
    for name in fieldcage_names_sorted:
        nodes_desc.append({
            "name": name,
            "boundary": name,
            "fixed": (name == top_fc_name),
        })

# Guard rings: all unknown
for name in guard_names_sorted:
    nodes_desc.append({
        "name": name,
        "boundary": name,
        "fixed": False,
    })

# Cathode node: fixed from BC_Cathode
nodes_desc.append({
    "name": "Cathode",
    "boundary": "BC_Cathode",
    "fixed": True,
})

# Build edges according to the same topology as solve_fieldcage_network
if n_fc > 0:
    fc = fieldcage_names_sorted
    guards = guard_names_sorted

    if n_fc < 8 or n_guard == 0:
        # Simple series chain: all R1, last to Cathode via R_C
        for i in range(n_fc - 1):
            edges_desc.append({"n1": fc[i], "n2": fc[i+1], "R": "R1"})
        edges_desc.append({"n1": fc[-1], "n2": "Cathode", "R": "R_C"})
    else:
        # First 5 with R1: FC1–2–3–4–5
        for i in range(0, 4):
            edges_desc.append({"n1": fc[i], "n2": fc[i+1], "R": "R1"})

        # Split at FC_5 via R2 for guard R3 for ring
        if n_fc > 5:
            edges_desc.append({"n1": fc[4], "n2": fc[5], "R": "R3"})
        if n_guard > 0:
            edges_desc.append({"n1": fc[4], "n2": guards[0], "R": "R2"})

        # Middle with R3
        merge_index = n_fc - 3  # 4th from last
        # Main branch: FC_6 ... FC_merge with R3
        for i in range(5, merge_index):
            edges_desc.append({"n1": fc[i], "n2": fc[i+1], "R": "R3"})
        # Guard branch: Guard_1 ... Guard_last with R3
        for j in range(0, n_guard - 1):
            edges_desc.append({"n1": guards[j], "n2": guards[j+1], "R": "R3"})
        # Merge: Guard_last --R2--> FC_merge
        if n_guard > 0:
            edges_desc.append({"n1": guards[-1], "n2": fc[merge_index], "R": "R2"})
        # Last 3 with R1: FC_merge..FC_last
        for i in range(merge_index, n_fc - 1):
            edges_desc.append({"n1": fc[i], "n2": fc[i+1], "R": "R1"})
        # Last FC to Cathode via R_C
        edges_desc.append({"n1": fc[-1], "n2": "Cathode", "R": "R_C"})


# ============================================================
# 5. WRITE config_autogen.yaml
# ============================================================

with output_path.open("w") as out:
    # materials
    out.write("materials:\n")
    for m in sorted(materials):
        v = materials[m]
        out.write(f"  {m}:\n")
        out.write(f"    attr_id: {v['attr_id']}\n")
        out.write(f"    epsilon_r: {v['epsilon_r']}\n")

    # boundaries
    out.write("\nboundaries:\n")

    # fixed boundaries first (Anode, Gate, Cathode, etc. - last is field cage 1)
    for name in sorted(boundaries_fixed_out):
        v = boundaries_fixed_out[name]
        out.write(f"  {name}:\n")
        out.write(f"    bdr_id: {v['bdr_id']}\n")
        out.write(f"    type: {v['type']}\n")
        out.write(f"    value: {v['value']}\n")

    # then all FieldCage_*/FieldCageGuard_* rings
    def sort_key_ring(n):
        # sort FieldCage_N..1, then FieldCageGuard_M..1 in index order (descending)
        if n.startswith("FieldCageGuard_"):
            idx = extract_index(n, "FieldCageGuard_") or 0
            return (1, -idx)
        else:
            idx = extract_index(n, "FieldCage_") or 0
            return (0, -idx)

    for name in sorted(boundaries_rings, key=sort_key_ring):
        if name in boundaries_fixed_out:
            continue   # skip FieldCage_1 once
        v = boundaries_rings[name]
        out.write(f"  {name}:\n")
        out.write(f"    bdr_id: {v['bdr_id']}\n")
        out.write(f"    type: {v['type']}\n")
        out.write(f"    value: {v['value']}\n")

    # fieldcage_network section
    out.write("\nfieldcage_network:\n")
    out.write("  enabled: true\n")

    # nodes
    out.write("  nodes:\n")
    for nd in nodes_desc:
        out.write(f"    - name: \"{nd['name']}\"\n")
        out.write(f"      boundary: \"{nd['boundary']}\"\n")
        out.write(f"      fixed: {str(nd['fixed']).lower()}\n")

    # resistor values
    out.write("  resistors:\n")
    out.write(f"    R1: {R1}\n")
    out.write(f"    R2: {R2}\n")
    out.write(f"    R3: {R3}\n")
    out.write(f"    R_C: {R_C}\n")

    # edges
    out.write("  edges:\n")
    for e in edges_desc:
        out.write(f"    - n1: \"{e['n1']}\"\n")
        out.write(f"      n2: \"{e['n2']}\"\n")
        out.write(f"      R: \"{e['R']}\"\n")

print(f"[done] Wrote {output_path} with {len(boundaries_rings)} FieldCage/guard boundaries "
      f"and {len(nodes_desc)} nodes / {len(edges_desc)} edges in fieldcage_network.")
