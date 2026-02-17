#pragma once

#include <iosfwd> 
#include <string>
#include "mfem.hpp"
#include "mfem/linalg/mumps.hpp"
#include "boundary_conditions.h"
#include "Config.h"
using namespace mfem;

std::unique_ptr<mfem::ParGridFunction> SolvePoisson(ParFiniteElementSpace &pfes, const BoundaryConditionGroups BCs, const std::shared_ptr<const Config>& cfg, std::string dir_overwrite = "");

// Residual file logging instead of terminal dump
class ResidualFileMonitor : public mfem::IterativeSolverMonitor
{
public:
    // The stream is not owned; caller must keep it alive.
    explicit ResidualFileMonitor(std::ostream &os);

    // Called by MFEM on each residual evaluation.
    void MonitorResidual(int it,
                         mfem::real_t norm,
                         const mfem::Vector &r,
                         bool final) override;

    // This is far more verbose and just discarded
    void MonitorSolution(int it,
                         mfem::real_t norm,
                         const mfem::Vector &x,
                         bool final) override;

private:
    std::ostream &out_;
};
