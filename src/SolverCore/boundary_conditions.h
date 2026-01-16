#pragma once

#include "mfem.hpp"
#include "boundary_conditions.h"
#include "Config.h"
using namespace mfem;

// Struct holding all boundary attributes for BC's 
struct BoundaryConditionGroups
{
    // True boundary (external)
    Array<int> dirichlet_attr;
    Array<int> neumann_attr;

    bool has_dirichlet() const { return !dirichlet_attr.IsEmpty(); }
    bool has_neumann() const { return !neumann_attr.IsEmpty(); }
};


BoundaryConditionGroups GetBoundaryConditionGroups(const mfem::Mesh *mesh, const std::shared_ptr<const Config> &cfg);
void ApplyDirichletValues(GridFunction &V, const Array<int> &dirichlet_attr, const std::shared_ptr<const Config>&);
void ApplyNeumannValues(ParLinearForm &b, const Array<int> &neumann_attr, const std::shared_ptr<const Config>& cfg, mfem::Coefficient &w, std::vector<std::unique_ptr<mfem::Coefficient>> &owned_coeffs, std::vector<mfem::Array<int>> &owned_markers);
