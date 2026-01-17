def get_bc():
    # False = (top=max index). True = (top=min index).
    grounded = {"type": "dirichlet", "value": 0}
    fixed_boundaries = {
        "PTFE_Wall_Charge": {"type": "neumann", "value": 0 },#"depth_dependent": True, 
                               # "z_bot": -0.0008 - 1.5008, "z_top": -0.0008, 
                               # # For now taken from https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenonnt:ftoschi:electric_field_matching_xenonnt&s[]=wall&s[]=charge
                               # "value_bot": -1, "value_top": -5},
        "BC_TopScreeningElectrode":    grounded,
        "BC_AnodeElectrode":           grounded,
        "BC_GateElectrode":            grounded,

        # Intended TOP-ring voltage (will be rebound to actual top ring name in phys_map)
        "BC_TopFieldShaping":          grounded,

        "BC_CathodeElectrode":         grounded,
        "BC_BottomScreeningElectrode": grounded,

        "BC_AllPMTs":      grounded,
        "BC_Bell":         grounded,
        "BC_CopperRing":   grounded,
        "BC_Cryostat":     grounded,
    }
    return fixed_boundaries
