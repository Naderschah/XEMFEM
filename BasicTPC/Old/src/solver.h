#ifndef SOLVER_H
#define SOLVER_H

#include "mfem.hpp"
#include "constants.h"
#include "boundary_conditions.h"
using namespace mfem;

GridFunction SolvePoisson(FiniteElementSpace &fespace, const Array<int> &dirichlet_attr);


#endif
