#include "solver.h"

GridFunction SolvePoisson(FiniteElementSpace &fespace, const Array<int> &dirichlet_attr)
{
    GridFunction V(&fespace);
    V = 0.0;
    ApplyDirichletValues(V, dirichlet_attr);

    BilinearForm a(&fespace);
    ConstantCoefficient epsilon_coeff(epsilon0);
    a.AddDomainIntegrator(new DiffusionIntegrator(epsilon_coeff));
    a.Assemble();

    LinearForm b(&fespace);
    b.Assemble();

    Array<int> ess_tdof_list;
    fespace.GetEssentialTrueDofs(dirichlet_attr, ess_tdof_list);
    a.EliminateEssentialBC(ess_tdof_list, V, b);
    a.Finalize();

    
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, V, b, A, X, B);

    // Solve with CG
    CGSolver cg;
    
    cg.SetRelTol(0);
    cg.SetAbsTol(tol);
    cg.SetMaxIter(maxiter);
    cg.SetPrintLevel(printlevel);

    DSmoother prec(A);
    cg.SetPreconditioner(prec);
    cg.SetOperator(A);
    cg.Mult(B, X);

    a.RecoverFEMSolution(X, b, V);

    return V;   // safe: fespace alive in main
}
