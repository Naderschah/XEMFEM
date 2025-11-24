#pragma once

#include <iosfwd> 

#include "mfem.hpp"
#include "boundary_conditions.h"
using namespace mfem;

struct Config; // forward declaration - still used?
std::unique_ptr<mfem::ParGridFunction> SolvePoisson(ParFiniteElementSpace &pfes, const mfem::Array<int> &dirichlet_attr, const std::shared_ptr<const Config>& cfg);

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
