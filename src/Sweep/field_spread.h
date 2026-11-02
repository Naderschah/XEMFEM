#pragma once

#include "solver_api.h"
#include "Config.h"

struct FieldHomogeneityMetrics
{
    double integration_volume = std::numeric_limits<double>::quiet_NaN();
    double mean_magnitude = std::numeric_limits<double>::quiet_NaN();
    double lower_quantile = std::numeric_limits<double>::quiet_NaN();
    double upper_quantile = std::numeric_limits<double>::quiet_NaN();
    double relative_rms_magnitude = std::numeric_limits<double>::quiet_NaN();
    double relative_rms_radial = std::numeric_limits<double>::quiet_NaN();
    double relative_rms_axial = std::numeric_limits<double>::quiet_NaN();
    double mean_axial = std::numeric_limits<double>::quiet_NaN();
    double mean_abs_axial = std::numeric_limits<double>::quiet_NaN();
};

double computeFieldSpreadMetric(
    const Config &cfg,
    const SimulationResult &result,
    FieldHomogeneityMetrics *homogeneity = nullptr);
