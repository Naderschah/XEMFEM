#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "mfem.hpp"
#include "boundary_conditions.h"
using namespace mfem;

// Returns an Array of Dirichlet boundary attributes (top and bottom plates)
Array<int> GetDirichletAttributes(Mesh *mesh);

void ApplyDirichletValues(GridFunction &V, const Array<int> &dirichlet_attr);

#endif
