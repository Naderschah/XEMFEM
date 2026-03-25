#pragma once

#include "Config.h"
#include "mfem.hpp"

bool ApplyTMOPStep(mfem::ParMesh &pmesh,
                   const Config &cfg);
