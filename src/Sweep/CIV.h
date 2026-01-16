#pragma once

#include "Config.h"
#include "trace_fieldlines.h"

// ------------------ CIV Computation Methods ------------------
// Sample n random points and trace (requires too many samples)
double ComputeCIV_RandomSample(const Config           &cfg,
                               const SimulationResult &result);
// Sample n columns and find the highest charge insensitive element
// Everything below is charge insensitive
double ComputeCIV_ColumnSweep(const SimulationResult &sim,
                              const Config           &cfg);
// TODO comment
double ComputeCIV_Marching(const SimulationResult &sim, const Config &cfg);

double compute_civ(const Config &cfg, const SimulationResult &result);