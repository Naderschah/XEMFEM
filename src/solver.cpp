#include "solver.h"

// Internal Helper for axisymmetric
inline std::unique_ptr<mfem::Coefficient>
MakeAxisymWeightCoeff(bool axisymmetric, int radial_idx)
{
  if (!axisymmetric) {
    return std::make_unique<mfem::ConstantCoefficient>(1.0);
  }

  // FunctionCoefficient expects physical coords X
  auto w_of_X = [radial_idx](const mfem::Vector &X) -> double {
    double r = X(radial_idx);
    if (r < 0.0) r = 0.0;               // guard tiny negatives from roundoff
    return 2.0 * M_PI * r;
  };
  return std::make_unique<mfem::FunctionCoefficient>(w_of_X);
}

// Internal Helper for distributed
inline bool IsDistributed(const mfem::FiniteElementSpace &fes)
{
  return dynamic_cast<const mfem::ParFiniteElementSpace*>(&fes) != nullptr;
}

inline MPI_Comm GetComm(const mfem::FiniteElementSpace &fes)
{
#ifdef MFEM_USE_MPI
  if (auto pfes = dynamic_cast<const mfem::ParFiniteElementSpace*>(&fes))
    return pfes->GetParMesh()->GetComm();
#endif
  return MPI_COMM_SELF; // safe fallback for serial path
}

// Build ε(x) that is piecewise-constant over element attributes (volume tags)
static PWConstCoefficient BuildEpsilonPWConst(const Mesh &mesh, const std::shared_ptr<const Config>& cfg)
{
    // attributes are 1-based; we need a vector sized by the max attribute id
    const int max_attr = mesh.attributes.Max();  // e.g. 2004 in your case
    Vector eps_by_attr(max_attr);
    eps_by_attr = 0.0;

    // helper that writes using 1-based attributes
    auto set_eps = [&](int attr, double val)
    {
        MFEM_VERIFY(attr >= 1 && attr <= max_attr,
                    "Volume attribute id out of range.");
        eps_by_attr(attr - 1) = val;
    };

    // Set all materials epsilon from config 
    for (const auto& [name, mat] : cfg->materials)
    {
      // Skip materials that have no explicit attr_id (e.g. "Default")
      if (mat.id > 0)  set_eps(mat.id, mat.epsilon_r);
      else if (name == "Default")
      {
        // Fill in any remaining attributes that were not explicitly set
        for (int a = 1; a <= max_attr; ++a)
        {
          // TODO Verify the logic + check if we want to actually do this 
          if (eps_by_attr(a - 1) == 0.0)  eps_by_attr(a - 1) = mat.epsilon_r;
        }
      }
    }
    // ctor available in your MFEM: PWConstCoefficient(Vector &c)
    return PWConstCoefficient(eps_by_attr);
}

std::unique_ptr<mfem::GridFunction> SolvePoisson(mfem::FiniteElementSpace &fespace,
                                                const mfem::Array<int> &dirichlet_attr,
                                                const std::shared_ptr<const Config>& cfg)
{
  using namespace mfem;

  const Mesh &mesh = *fespace.GetMesh();
  PWConstCoefficient epsilon_pw = BuildEpsilonPWConst(mesh, cfg);

  // common handles (filled inside the single if/else)
  std::unique_ptr<GridFunction> V;        // GridFunction or ParGridFunction
  std::unique_ptr<BilinearForm> a;        // BilinearForm or ParBilinearForm
  std::unique_ptr<LinearForm>   b;        // LinearForm   or ParLinearForm
  OperatorHandle                A;        // SparseMatrix or HypreParMatrix
  Vector                        X, B;     // works for both
  std::function<std::unique_ptr<Solver>(OperatorHandle&)> make_prec;
  MPI_Comm                      comm = MPI_COMM_SELF;

  const bool par = (dynamic_cast<ParFiniteElementSpace*>(&fespace) != nullptr);

  if (par) {
    // ---------- parallel concrete types ----------
    auto &pfes = static_cast<ParFiniteElementSpace&>(fespace);
    V = std::make_unique<ParGridFunction>(&pfes);
    a = std::make_unique<ParBilinearForm>(&pfes);
    b = std::make_unique<ParLinearForm>(&pfes);
    comm = pfes.GetParMesh()->GetComm();

    // build the proper preconditioner later, after A exists
    make_prec = [=](OperatorHandle &Ah) -> std::unique_ptr<Solver>
    {
      auto *Ap = Ah.As<HypreParMatrix>();
      return std::make_unique<HypreBoomerAMG>(*Ap);
    };
  } else {
    // ---------- serial concrete types ----------
    V = std::make_unique<GridFunction>(&fespace);
    a = std::make_unique<BilinearForm>(&fespace);
    b = std::make_unique<LinearForm>(&fespace);
    comm = MPI_COMM_SELF;

    make_prec = [=](OperatorHandle &Ah) -> std::unique_ptr<Solver>
    {
      auto *As = Ah.As<SparseMatrix>();
      return std::make_unique<DSmoother>(*As);
    };
  }

  // ---------- shared body ----------
  *V = 0.0;
  ApplyDirichletValues(*V, dirichlet_attr, cfg);

  auto w = MakeAxisymWeightCoeff(cfg->solver.axisymmetric, 0);
  mfem::ProductCoefficient weps(*w, epsilon_pw);
  a->AddDomainIntegrator(new DiffusionIntegrator(weps));
  a->Assemble();                 // Finalize() not needed with OperatorHandle path

  b->Assemble();

  Array<int> ess_tdof;
  fespace.GetEssentialTrueDofs(dirichlet_attr, ess_tdof);

  a->FormLinearSystem(ess_tdof, *V, *b, A, X, B);   // fills A, X, B for both modes

  auto P = make_prec(A);          // DSmoother or HypreBoomerAMG, already decided

  CGSolver cg(comm);
  cg.SetOperator(*A.Ptr());       // OperatorHandle -> Operator&
  cg.SetPreconditioner(*P);
  cg.SetRelTol(cfg->solver.rtol);
  cg.SetAbsTol(cfg->solver.atol);
  cg.SetMaxIter(cfg->solver.maxiter);
  cg.SetPrintLevel(cfg->solver.printlevel);
  cg.Mult(B, X);

  a->RecoverFEMSolution(X, *b, *V);

  return V;
}
