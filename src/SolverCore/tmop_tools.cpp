#include "tmop_tools.h"

bool ApplyTMOPStep(mfem::ParMesh &,
                   const Config &cfg)
{
    if (!cfg.mesh.tmop.enable)
    {
        return false;
    }

    // TMOP integration is wired into the solve loop and config parser first.
    // The implementation will be filled in next against the MFEM TMOP API.
    return false;
}
