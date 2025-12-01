#include "solver.h"

#include <filesystem>
#include <iomanip>
#include "HYPRE_parcsr_ls.h"

// Writing residuals to file
ResidualFileMonitor::ResidualFileMonitor(std::ostream &os)
    : out_(os)
{
}

void ResidualFileMonitor::MonitorResidual(int it,
                                          mfem::real_t norm,
                                          const mfem::Vector &,
                                          bool final)
{
    out_ << it << " " << std::scientific << norm;
    if (final)
    {
        out_ << " final";
    }
    out_ << std::endl;
}
void ResidualFileMonitor::MonitorSolution(int,
                                          mfem::real_t,
                                          const mfem::Vector&,
                                          bool)
{
    // no-op
}





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

std::unique_ptr<mfem::ParGridFunction> SolvePoisson(ParFiniteElementSpace &pfes,
                                                const mfem::Array<int> &dirichlet_attr,
                                                const std::shared_ptr<const Config>& cfg)
{
  using namespace mfem;

  if (cfg->debug.debug) {std::cout << "[DEBUG] In Poisson Solver" << std::endl;}

  const std::filesystem::path directory(cfg->save_path);
  const std::filesystem::path log_path = directory / "solver.log";

  const Mesh &mesh = *pfes.GetMesh();
  PWConstCoefficient epsilon_pw = BuildEpsilonPWConst(mesh, cfg);
  // TODO add check that epsilon and dirichlet are full and buitl 
  // common handles (filled inside the single if/else)
  std::unique_ptr<ParGridFunction> V;     // GridFunction or ParGridFunction
  std::unique_ptr<ParBilinearForm> a;        // BilinearForm or ParBilinearForm
  std::unique_ptr<ParLinearForm>   b;        // LinearForm   or ParLinearForm
  OperatorHandle                A;        // SparseMatrix or HypreParMatrix
  Vector                        X, B;     // works for both
  std::function<std::unique_ptr<Solver>(OperatorHandle&)> make_prec;

  V = std::make_unique<ParGridFunction>(&pfes);
  a = std::make_unique<ParBilinearForm>(&pfes);
  b = std::make_unique<ParLinearForm>(&pfes);

  MPI_Comm comm = pfes.GetParMesh()->GetComm();

  // build the proper preconditioner later, after A exists
  make_prec = [=](OperatorHandle &Ah) -> std::unique_ptr<Solver>
  {
    auto *Ap = Ah.As<HypreParMatrix>();
    auto prec = std::make_unique<HypreBoomerAMG>(*Ap);

    // Make hypre actually print something (level > 0)
    prec->SetPrintLevel(cfg->solver.printlevel);

    // Get underlying hypre solver object and set its print file
    {
      HYPRE_Solver hsolver = (HYPRE_Solver)(*prec);
      HYPRE_BoomerAMGSetPrintFileName(hsolver, log_path.c_str());
    }

    return prec;
  };

  *V = 0.0;
  ApplyDirichletValues(*V, dirichlet_attr, cfg);

  auto w = MakeAxisymWeightCoeff(cfg->solver.axisymmetric, 0);
  mfem::ProductCoefficient weps(*w, epsilon_pw);
  a->AddDomainIntegrator(new DiffusionIntegrator(weps));
  a->Assemble();                

  b->Assemble();

  Array<int> ess_tdof;
  pfes.GetEssentialTrueDofs(dirichlet_attr, ess_tdof);

  a->FormLinearSystem(ess_tdof, *V, *b, A, X, B);  

  auto P = make_prec(A);
  // It will error on setup, for ADS/AMS this happens in l1 row norm producing singular coarse-grid matrices causing the errors
  // but this is handled fine in those cases, so I am hoping it is here as well TODO Really should figure out the origin of this 
  // Future debug note: This only started happening when I compiled hypre inside the Dockerfile, up until that point this didnt 
  // happen, but OpenMP didnt work either (ie always 1 thread never more)  
  if (auto *hypre_prec = dynamic_cast<mfem::HypreSolver*>(P.get()))
  {   
    if (cfg->debug.printHypreWarnings) hypre_prec->SetErrorMode(mfem::HypreSolver::WARN_HYPRE_ERRORS);
    else hypre_prec->SetErrorMode(mfem::HypreSolver::IGNORE_HYPRE_ERRORS);
  }

  CGSolver cg(comm);
  cg.SetOperator(*A.Ptr());     
  cg.SetPreconditioner(*P);
  cg.SetRelTol(cfg->solver.rtol);
  cg.SetAbsTol(cfg->solver.atol);
  cg.SetMaxIter(cfg->solver.maxiter);
  cg.SetPrintLevel(cfg->solver.printlevel);
  // Print to file TODO need to homogenixze my saving customs
  // already  have a save path but need to pass a per run save path
  std::ofstream it_log(log_path);
  ResidualFileMonitor monitor(it_log);
  cg.SetMonitor(monitor);

  cg.Mult(B, X);

  a->RecoverFEMSolution(X, *b, *V);


  return V;
}
