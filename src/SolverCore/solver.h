#pragma once

#include "mfem.hpp"
#include "boundary_conditions.h"
using namespace mfem;

struct Config; // forward declaration - still used?
std::unique_ptr<mfem::GridFunction> SolvePoisson(mfem::FiniteElementSpace &fespace, const mfem::Array<int> &dirichlet_attr, const std::shared_ptr<const Config>& cfg);
