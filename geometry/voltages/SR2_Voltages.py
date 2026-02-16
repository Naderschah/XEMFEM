# https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenonnt:analysis:available_data_for_sr2#dates
def get_bc():
    grounded = {"type": "dirichlet", "value": 0}
    fixed_boundaries = {
        "PTFE_Wall_Charge": {"type": "neumann", "depth_dependent": True, 
                                # For now taken from https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenonnt:ftoschi:electric_field_matching_xenonnt&s[]=wall&s[]=charge
                                "z_bot": -0.0008 - 1.5008, "z_top": -0.0008, "value_bot": -0.1e-6, "value_top": -0.5e-6},
        "BC_r0":                        {"type": "neumann", "value": 0},
        "BC_TopScreeningElectrode":    {"type": "dirichlet", "value": -100},
        "BC_AnodeElectrode":           {"type": "dirichlet", "value": 4850},
        "BC_GateElectrode":            {"type": "dirichlet", "value": +300},

        # Intended TOP-ring voltage (will be rebound to actual top ring name in phys_map)
        "BC_TopFieldShaping":          {"type": "dirichlet", "value": 650},

        "BC_CathodeElectrode":         {"type": "dirichlet", "value": -2750},
        "BC_BottomScreeningElectrode": {"type": "dirichlet", "value": -2750},

        "BC_AllPMTs":                  {"type": "dirichlet", "value": -1300},
        "BC_Bell":                     grounded,
        "BC_CopperRing":               grounded,
        "BC_Cryostat":                 grounded,
    }
    return fixed_boundaries