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

ElectricFieldPostprocessor::ElectricFieldPostprocessor(FiniteElementSpace &V_h1)
  : mesh_(*V_h1.GetMesh()),
    fes_h1_(&V_h1),               // non-const pointer
    dim_(mesh_.Dimension()),
    p_(V_h1.GetFE(0)->GetOrder())
{
  // Vector L2 space for E: order p-1 (clamped at 0)
  auto *fecL2_vec = new L2_FECollection(std::max(p_ - 1, 0), dim_);
  vfes_L2_.reset(new FiniteElementSpace(&mesh_, fecL2_vec, dim_)); // vector L2 (vdim = dim)

  // Scalar L2 space for components Ex/Ey/Ez: same order as vector L2
  sfes_L2_.reset(new FiniteElementSpace(&mesh_, fecL2_vec));       // scalar L2(p-1)

  // Scalar L2(P0) space for |E| (Emag): piecewise constant, one DOF per element
  auto *fecL2_p0 = new L2_FECollection(0, dim_);
  sfes_L2_mag_.reset(new FiniteElementSpace(&mesh_, fecL2_p0));    // scalar L2(P0)

  // Discrete gradient H1 -> L2^dim (expects non-const fes pointers)
  gradOp_.reset(new DiscreteLinearOperator(fes_h1_, vfes_L2_.get()));
  gradOp_->AddDomainInterpolator(new GradientInterpolator);
  gradOp_->Assemble();
  gradOp_->Finalize();
}

std::unique_ptr<GridFunction> ElectricFieldPostprocessor::MakeE() const
{
  return std::make_unique<GridFunction>(vfes_L2_.get()); 
}

std::unique_ptr<GridFunction> ElectricFieldPostprocessor::MakeEmag() const
{
  return std::make_unique<GridFunction>(sfes_L2_mag_.get());
}

void ElectricFieldPostprocessor::ComputeElectricField(const GridFunction &V,
                                                      GridFunction &E_out,
                                                      double scale) const
{
  MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
    "ComputeElectricField: E_out must be a vector L2 GridFunction (vdim=dim).");

  gradOp_->Mult(V, E_out); // E_out = ∇V
  E_out *= scale;          // -1.0 for E = -∇V
  
}

void ElectricFieldPostprocessor::ComputeFieldMagnitude(const GridFunction &E,
                                                       GridFunction &Emag) const
{
  NormOfVectorGF normE(E, dim_);
  Emag.ProjectCoefficient(normE); // L2(P0) - no negative values
}

void ElectricFieldPostprocessor::SaveComponents(const GridFunction &E,
                                                const std::string &prefix) const
{
  // Components Ex/Ey/Ez stored in scalar L2(p-1)
  FiniteElementSpace *scalar_space = sfes_L2_.get();

  for (int d = 0; d < dim_; ++d)
  {
    VectorComponentOfGF comp_cf(E, d, dim_);
    GridFunction comp_gf(scalar_space);
    comp_gf.ProjectCoefficient(comp_cf);

    std::string fname = prefix + (d == 0 ? "_ex.gf"
                                         : (d == 1 ? "_ey.gf" : "_ez.gf"));
    std::ofstream ofs(fname);
    comp_gf.Save(ofs);
  }
}
void ElectricFieldPostprocessor::LoadE(GridFunction &E_out,
                                       const std::string &prefix)
{
    MFEM_VERIFY(E_out.FESpace() == vfes_L2_.get(),
                "LoadE: E_out must be a vector L2 GridFunction (vdim = dim).");

    FiniteElementSpace *scalar_space = sfes_L2_.get();

    // Load scalar components Ex/Ey/Ez as GridFunctions on sfes_L2_
    std::unique_ptr<GridFunction> comp_gf[3];

    for (int d = 0; d < dim_; ++d)
    {
        std::string fname = prefix + (d == 0 ? "_ex.gf"
                                 : (d == 1 ? "_ey.gf" : "_ez.gf"));

        std::ifstream ifs(fname);
        if (!ifs)
        {
            mfem::mfem_error(("LoadE: cannot open component file '" + fname + "'").c_str());
        }

        comp_gf[d] = std::make_unique<GridFunction>(scalar_space);
        comp_gf[d]->Load(ifs);
    }

    // Vector coefficient constructed from scalar component grid functions
    class VectorFromScalarGFs : public VectorCoefficient
    {
    public:
        const GridFunction *ex_;
        const GridFunction *ey_;
        const GridFunction *ez_;
        int dim_;

        VectorFromScalarGFs(int dim,
                            const GridFunction *ex,
                            const GridFunction *ey,
                            const GridFunction *ez)
            : VectorCoefficient(dim),
              ex_(ex), ey_(ey), ez_(ez), dim_(dim)
        {}

        virtual void Eval(Vector &V,
                          ElementTransformation &T,
                          const IntegrationPoint &ip) override
        {
            // Build temporary scalar coefficients on the fly
            if (dim_ >= 1 && ex_)
            {
                GridFunctionCoefficient cf_ex(ex_);
                V(0) = cf_ex.Eval(T, ip);
            }
            if (dim_ >= 2 && ey_)
            {
                GridFunctionCoefficient cf_ey(ey_);
                V(1) = cf_ey.Eval(T, ip);
            }
            if (dim_ >= 3 && ez_)
            {
                GridFunctionCoefficient cf_ez(ez_);
                V(2) = cf_ez.Eval(T, ip);
            }
        }
    };

    VectorFromScalarGFs vec_cf(
        dim_,
        comp_gf[0] ? comp_gf[0].get() : nullptr,
        (dim_ > 1 && comp_gf[1]) ? comp_gf[1].get() : nullptr,
        (dim_ > 2 && comp_gf[2]) ? comp_gf[2].get() : nullptr);

    // Project the vector coefficient into the vector L2 space of E_out
    E_out.ProjectCoefficient(vec_cf);
}
// -------------------- Optional simple wrappers ------------------------------
namespace { static std::unique_ptr<ElectricFieldPostprocessor> g_post; }

void InitFieldPostprocessor(FiniteElementSpace &V_h1)
{
  g_post = std::make_unique<ElectricFieldPostprocessor>(V_h1);
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
