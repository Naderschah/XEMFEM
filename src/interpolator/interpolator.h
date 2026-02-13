#pragma once 

#include "Config.h"
#include "solver_api.h"

#include <vector>
#include <string>
#include <hdf5.h>
#include <limits>
#include "mfem.hpp"
#include <mfem/fem/gslib.hpp>
#include "path_handler.h"

using namespace mfem;

struct GridSample
{
  int dim;
  int Nx, Ny, Nz;
  Vector origin;   // size dim
  Vector spacing;  // size dim
  // output arrays sized (Nx*Ny*Nz):
  std::vector<double> Ex, Ey, Ez;  // for dim==3
  std::vector<double> Emag;
  std::vector<uint8_t> valid;      // 1 if inside (or accepted), 0 otherwise
};

void SampleEFieldOnCartesianGrid(ParMesh &pmesh,
                                const ParGridFunction &V,
                                const int Nx, const int Ny, const int Nz,
                                GridSample &out,
                                Config cfg,
                                const bool H1_project,
                                const bool accept_surface_projection);

int do_interpolate(Config cfg);