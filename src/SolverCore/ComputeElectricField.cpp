#include "ComputeElectricField.h"
#include <fstream>
#include <algorithm>

using namespace mfem;

namespace {
// ||Vector GridFunction|| coefficient for projection
class NormOfVectorGF : public Coefficient
{
public:
  explicit NormOfVectorGF(const GridFunction &E, int vdim)
    : E_(&E), val_(vdim), Ecf_(E_) {}

  double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
  {
    Ecf_.Eval(val_, T, ip);
    return val_.Norml2();
  }
private:
  const GridFunction *E_;
  mutable Vector      val_;
  VectorGridFunctionCoefficient Ecf_;
};

// Extract one component of a vector GridFunction as a scalar Coefficient
class VectorComponentOfGF : public Coefficient
{
public:
  VectorComponentOfGF(const GridFunction &E, int comp, int vdim)
    : E_(&E), comp_(comp), tmp_(vdim), Ecf_(E_) {}

  double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
  {
    Ecf_.Eval(tmp_, T, ip);
    MFEM_ASSERT(comp_ >= 0 && comp_ < tmp_.Size(), "Invalid component index");
    return tmp_[comp_];
  }
private:
  const GridFunction *E_;
  int                 comp_;
  mutable Vector      tmp_;
  VectorGridFunctionCoefficient Ecf_;
};
} // namespace

// -------------------- ElectricFieldPostprocessor ----------------------------

ElectricFieldPostprocessor::ElectricFieldPostprocessor(FiniteElementSpace &V_h1,
                                                       bool smooth_output)
  : mesh_(*V_h1.GetMesh()),
    fes_h1_(&V_h1),               // non-const pointer
    dim_(mesh_.Dimension()),
    p_(V_h1.GetFE(0)->GetOrder()),
    smooth_(smooth_output)
{
  // L2 spaces (order p-1 for ∇ of H1(p), clamped at 0)
  auto *fecL2 = new L2_FECollection(std::max(p_-1, 0), dim_);
  sfes_L2_.reset(new FiniteElementSpace(&mesh_, fecL2));          // scalar L2
  vfes_L2_.reset(new FiniteElementSpace(&mesh_, fecL2, dim_));    // vector L2 (vdim = dim)

  // Discrete gradient H1 -> L2^dim (expects non-const fes pointers)
  gradOp_.reset(new DiscreteLinearOperator(fes_h1_, vfes_L2_.get()));
  gradOp_->AddDomainInterpolator(new GradientInterpolator);
  gradOp_->Assemble();
  gradOp_->Finalize();

  if (smooth_)
  {
    auto *fech1 = new H1_FECollection(std::max(p_-1, 0), dim_);
    sfes_H1_.reset(new FiniteElementSpace(&mesh_, fech1));          // scalar H1
    vfes_H1_.reset(new FiniteElementSpace(&mesh_, fech1, dim_));    // vector H1
  }
}

std::unique_ptr<GridFunction> ElectricFieldPostprocessor::MakeE() const
{
  if (!smooth_) { return std::make_unique<GridFunction>(vfes_L2_.get()); }
  return std::make_unique<GridFunction>(vfes_H1_.get());
}

std::unique_ptr<GridFunction> ElectricFieldPostprocessor::MakeEmag() const
{
  if (!smooth_) { return std::make_unique<GridFunction>(sfes_L2_.get()); }
  return std::make_unique<GridFunction>(sfes_H1_.get());
}

void ElectricFieldPostprocessor::ComputeElectricField(const GridFunction &V,
                                                      GridFunction &E_out,
                                                      double scale) const
{
  if (!smooth_)
  {
    MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
      "ComputeElectricField: E_out must be a vector L2 GridFunction (vdim=dim).");

    gradOp_->Mult(V, E_out); // E_out = ∇V
    E_out *= scale;          // typically -1.0 for E = -∇V
  }
  else
  {
    GridFunction E_l2(vfes_L2_.get());
    gradOp_->Mult(V, E_l2);
    E_l2 *= scale;

    MFEM_VERIFY(E_out.FESpace() == vfes_H1_.get(),
      "ComputeElectricField(smooth=true): E_out must be vector H1 (vdim=dim).");

    VectorGridFunctionCoefficient Ec(&E_l2);
    E_out.ProjectCoefficient(Ec);
  }
}

void ElectricFieldPostprocessor::ComputeFieldMagnitude(const GridFunction &E,
                                                       GridFunction &Emag) const
{
  NormOfVectorGF normE(E, dim_);
  Emag.ProjectCoefficient(normE);
}

void ElectricFieldPostprocessor::SaveComponents(const GridFunction &E,
                                                const std::string &prefix) const
{
  // Use non-const pointer type here (GridFunction expects FiniteElementSpace*)
  FiniteElementSpace *scalar_space =
      (!smooth_ ? sfes_L2_.get() : sfes_H1_.get());

  for (int d = 0; d < dim_; ++d)
  {
    VectorComponentOfGF comp_cf(E, d, dim_);
    GridFunction comp_gf(scalar_space);   // OK now: non-const FES*
    comp_gf.ProjectCoefficient(comp_cf);

    std::string fname = prefix + (d == 0 ? "_ex.gf" : (d == 1 ? "_ey.gf" : "_ez.gf"));
    std::ofstream ofs(fname);
    comp_gf.Save(ofs);
  }
}

// -------------------- Optional simple wrappers ------------------------------
namespace { static std::unique_ptr<ElectricFieldPostprocessor> g_post; }

void InitFieldPostprocessor(FiniteElementSpace &V_h1, bool smooth_output)
{
  g_post = std::make_unique<ElectricFieldPostprocessor>(V_h1, smooth_output);
}

std::unique_ptr<GridFunction> CreateE()
{
  MFEM_VERIFY(g_post, "CreateE: call InitFieldPostprocessor() first.");
  return g_post->MakeE();
}

std::unique_ptr<GridFunction> CreateEmag()
{
  MFEM_VERIFY(g_post, "CreateEmag: call InitFieldPostprocessor() first.");
  return g_post->MakeEmag();
}

void ComputeElectricField(GridFunction &V, GridFunction &E, double scale)
{
  MFEM_VERIFY(g_post, "ComputeElectricField: call InitFieldPostprocessor() first.");
  g_post->ComputeElectricField(V, E, scale);
}

void ComputeFieldMagnitude(const GridFunction &E, GridFunction &Emag)
{
  MFEM_VERIFY(g_post, "ComputeFieldMagnitude: call InitFieldPostprocessor() first.");
  g_post->ComputeFieldMagnitude(E, Emag);
}

void SaveEComponents(const GridFunction &E, const std::string &prefix)
{
  MFEM_VERIFY(g_post, "SaveEComponents: call InitFieldPostprocessor() first.");
  g_post->SaveComponents(E, prefix);
}
