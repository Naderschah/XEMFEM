#pragma once
#ifndef COMPUTE_ELECTRIC_FIELD_H
#define COMPUTE_ELECTRIC_FIELD_H

#include "mfem.hpp"
#include <memory>
#include <string>

class ElectricFieldPostprocessor
{
public:
  explicit ElectricFieldPostprocessor(mfem::ParFiniteElementSpace &V_h1);

  std::unique_ptr<mfem::ParGridFunction> MakeE() const;     // vector field (vdim = dim)
  std::unique_ptr<mfem::ParGridFunction> MakeEmag() const;  // scalar field

  void ComputeElectricField(const mfem::ParGridFunction &V,
                            mfem::ParGridFunction &E_out,
                            double scale = -1.0) const;

  void ComputeFieldMagnitude(const mfem::ParGridFunction &E,
                             mfem::ParGridFunction &Emag) const;

  void SaveComponents(const mfem::ParGridFunction &E,
                      const std::string &prefix) const;

  void LoadE(mfem::ParGridFunction &E_out,
             const std::string &prefix);

  int Dimension() const { return dim_; }

private:
  // NOTE: these must be non-const pointers for MFEM APIs that expect non-const.
  mfem::ParMesh                 &mesh_;
  mfem::ParFiniteElementSpace  *fes_h1_;   // space of V (non-const ptr)
  int                            dim_;
  int                            p_;

  // Owned FE spaces (non-const pointers returned by get()).
  std::unique_ptr<mfem::ParFiniteElementSpace> sfes_L2_; // scalar L2
  std::unique_ptr<mfem::ParFiniteElementSpace> vfes_L2_; // vector L2 (vdim = dim)

  std::unique_ptr<mfem::ParFiniteElementSpace> sfes_L2_mag_; // scalar L2(p=0)

  std::unique_ptr<mfem::DiscreteLinearOperator> gradOp_; // H1 → L2^dim
};

// -------- Optional simple wrappers (keep your old call style) --------
void InitFieldPostprocessor(mfem::ParFiniteElementSpace &V_h1);
std::unique_ptr<mfem::ParGridFunction> CreateE();     // vector container
std::unique_ptr<mfem::ParGridFunction> CreateEmag();  // scalar container
void ComputeElectricField(mfem::ParGridFunction &V, mfem::ParGridFunction &E, double scale = -1.0);
void ComputeFieldMagnitude(const mfem::ParGridFunction &E, mfem::ParGridFunction &Emag);
void SaveEComponents(const mfem::ParGridFunction &E, const std::string &prefix);

#endif // COMPUTE_ELECTRIC_FIELD_H
