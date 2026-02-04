#include "ComputeElectricField.h"
#include <fstream>
#include <algorithm>

using namespace mfem;

namespace {
// ||Vector GridFunction|| coefficient for projection
class NormOfVectorGF : public Coefficient
{
public:
  explicit NormOfVectorGF(const ParGridFunction &E, int vdim)
    : E_(&E), val_(vdim), Ecf_(E_) {}

  double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
  {
    Ecf_.Eval(val_, T, ip);
    return val_.Norml2();
  }
private:
  const ParGridFunction *E_;
  mutable Vector        val_;
  VectorGridFunctionCoefficient Ecf_;
};

// Extract one component of a vector GridFunction as a scalar Coefficient
class VectorComponentOfGF : public Coefficient
{
public:
  VectorComponentOfGF(const ParGridFunction &E, int comp, int vdim)
    : E_(&E), comp_(comp), tmp_(vdim), Ecf_(E_) {}

  double Eval(ElementTransformation &T, const IntegrationPoint &ip) override
  {
    Ecf_.Eval(tmp_, T, ip);
    MFEM_ASSERT(comp_ >= 0 && comp_ < tmp_.Size(), "Invalid component index");
    return tmp_[comp_];
  }
private:
  const ParGridFunction *E_;
  int                 comp_;
  mutable Vector      tmp_;
  VectorGridFunctionCoefficient Ecf_;
};
} // namespace

// -------------------- ElectricFieldPostprocessor ----------------------------

ElectricFieldPostprocessor::ElectricFieldPostprocessor(ParFiniteElementSpace &V_h1)
  : mesh_(dynamic_cast<mfem::ParMesh&>(*V_h1.GetMesh())),
    fes_h1_(&V_h1),               // non-const pointer
    dim_(mesh_.Dimension()),
    p_(V_h1.GetFE(0)->GetOrder())
{
  // Vector L2 space for E: order p-1 (clamped at 0)
  auto *fecL2_vec = new L2_FECollection(std::max(p_ - 1, 0), dim_);
  vfes_L2_.reset(new ParFiniteElementSpace(&mesh_, fecL2_vec, dim_)); // vector L2 (vdim = dim)

  // Scalar L2 space for components Ex/Ey/Ez: same order as vector L2
  sfes_L2_.reset(new ParFiniteElementSpace(&mesh_, fecL2_vec));       // scalar L2(p-1)

  // Scalar L2(P0) space for |E| (Emag): piecewise constant, one DOF per element
  auto *fecL2_p0 = new L2_FECollection(0, dim_);
  sfes_L2_mag_.reset(new ParFiniteElementSpace(&mesh_, fecL2_p0));    // scalar L2(P0)

  // Discrete gradient H1 -> L2^dim (expects non-const fes pointers)
  gradOp_.reset(new DiscreteLinearOperator(fes_h1_, vfes_L2_.get()));
  gradOp_->AddDomainInterpolator(new GradientInterpolator);
  gradOp_->Assemble();
  gradOp_->Finalize();
}

std::unique_ptr<ParGridFunction> ElectricFieldPostprocessor::MakeE() const
{
  return std::make_unique<ParGridFunction>(vfes_L2_.get()); 
}

std::unique_ptr<ParGridFunction> ElectricFieldPostprocessor::MakeEmag() const
{
  return std::make_unique<ParGridFunction>(sfes_L2_mag_.get());
}

void ElectricFieldPostprocessor::ComputeElectricField(const ParGridFunction &V,
                                                      ParGridFunction &E_out,
                                                      double scale) const
{
  MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
    "ComputeElectricField: E_out must be a vector L2 GridFunction (vdim=dim).");

  gradOp_->Mult(V, E_out); // E_out = ∇V
  E_out *= scale;          // -1.0 for E = -∇V
  
}

void ElectricFieldPostprocessor::ComputeFieldMagnitude(const ParGridFunction &E,
                                                       ParGridFunction &Emag) const
{
  NormOfVectorGF normE(E, dim_);
  Emag.ProjectCoefficient(normE); // L2(P0) - no negative values
}

void ElectricFieldPostprocessor::SaveComponents(const ParGridFunction &E,
                                                const std::string &prefix) const
{
  // Components Ex/Ey/Ez stored in scalar L2(p-1)
  ParFiniteElementSpace *scalar_space = sfes_L2_.get();

  for (int d = 0; d < dim_; ++d)
  {
    VectorComponentOfGF comp_cf(E, d, dim_);
    ParGridFunction comp_gf(scalar_space);
    comp_gf.ProjectCoefficient(comp_cf);

    std::string fname = prefix + (d == 0 ? "_ex.gf"
                                         : (d == 1 ? "_ey.gf" : "_ez.gf"));
    comp_gf.SaveAsSerial(fname.c_str(), /*precision=*/16, /*save_rank=*/0);
  }
}

void ElectricFieldPostprocessor::LoadE(mfem::ParGridFunction &E_out,
                                       const std::string &prefix)
{
    MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
                "LoadE: E_out must be a vector L2 GridFunction (vdim = dim).");

    const std::string fname = prefix + ".gf"; // e.g. ".../E.gf"
    std::ifstream ifs(fname);
    if (!ifs)
    {
        mfem::mfem_error(("LoadE: cannot open '" + fname + "'").c_str());
    }

    E_out.Load(ifs);

    // Optional sanity checks (cheap)
    MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
                "LoadE: loaded data does not match expected vector L2 space.");
}

// -------------------- Optional simple wrappers ------------------------------
namespace { static std::unique_ptr<ElectricFieldPostprocessor> g_post; }

void InitFieldPostprocessor(ParFiniteElementSpace &V_h1)
{
  g_post = std::make_unique<ElectricFieldPostprocessor>(V_h1);
}

std::unique_ptr<ParGridFunction> CreateE()
{
  MFEM_VERIFY(g_post, "CreateE: call InitFieldPostprocessor() first.");
  return g_post->MakeE();
}

std::unique_ptr<ParGridFunction> CreateEmag()
{
  MFEM_VERIFY(g_post, "CreateEmag: call InitFieldPostprocessor() first.");
  return g_post->MakeEmag();
}

void ComputeElectricField(ParGridFunction &V, ParGridFunction &E, double scale)
{
  MFEM_VERIFY(g_post, "ComputeElectricField: call InitFieldPostprocessor() first.");
  g_post->ComputeElectricField(V, E, scale);
}

void ComputeFieldMagnitude(const ParGridFunction &E, ParGridFunction &Emag)
{
  MFEM_VERIFY(g_post, "ComputeFieldMagnitude: call InitFieldPostprocessor() first.");
  g_post->ComputeFieldMagnitude(E, Emag);
}

void SaveEComponents(const ParGridFunction &E, const std::string &prefix)
{
  MFEM_VERIFY(g_post, "SaveEComponents: call InitFieldPostprocessor() first.");
  g_post->SaveComponents(E, prefix);
}
