#pragma once

#include "Config.h"
#include "trace_fieldlines.h"

// ------------------ CIV Computation Methods ------------------
// Sample n random points and trace (requires too many samples)
double ComputeCIV_RandomSample(const Config           &cfg,
                               const SimulationResult &result);
double ComputeCIV_FixedGrid(const Config            &cfg,
                            const SimulationResult &result);
double ComputeCIV_AdaptiveGrid(const Config            &cfg,
                               const SimulationResult &result);
double ComputeCIV_RowSweep(const Config            &cfg,
                           const SimulationResult &result);