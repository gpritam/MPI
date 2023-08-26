#include "HypreSolvers.h"

//_______________________________________________________________________________
// This is an example (not recommended) of how we can modify things about AMG that
// affect the solve phase based on how FlexGMRES is doing. For another preconditioner 
// it may make sense to modify the tolerance.
//_______________________________________________________________________________
int hypre_FlexGMRESModifyPCAMGExample ( void *precond_data,
                                        int iterations,
                                        double rel_residual_norm )
{
    HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data,(rel_residual_norm > 0.1 ? 10 : 1));
    
    return 0;
}

//________________________________________________________________________________________________
// Solve Ax = rhsHypre by Hypre
// 
// solver_id = 0 : AMG
// solver_id = 1 : Conjugate Gradient method with preconditioner
// solver_id = 2 : Conjugate Gradient method with Parasails preconditioner
// solver_id = 3 : Conjugate Gradient method with AMG preconditioner
// solver_id = 4 : Flexible GMRES with AMG preconditioner
//________________________________________________________________________________________________
void HypreSolve ( HYPRE_ParCSRMatrix AHypre, 
                  HYPRE_ParVector rhsHypre, 
                  HYPRE_ParVector solutionHypre, 
                  int &iterations, 
                  double &residual, 
                  const int solver_id, 
                  const double tolerance,
                  const int MaximumCycle,
                  const bool PrintSolutionInformation )
{
    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    HYPRE_Solver solver, preconditioner;
    
    if (solver_id == 0)
    {
        HYPRE_BoomerAMGCreate(&solver);
        
        // Set some parameters (See Reference Manual for more parameters)
        if (PrintSolutionInformation)
            HYPRE_BoomerAMGSetPrintLevel(solver, 3);     // print solve info + parameters
        
        HYPRE_BoomerAMGSetOldDefault(solver);            // Falgout coarsening with modified classical interpolaiton
        HYPRE_BoomerAMGSetRelaxType(solver, 3);          // G-S/Jacobi hybrid relaxation
        HYPRE_BoomerAMGSetRelaxOrder(solver, 1);         // uses C/F relaxation
        HYPRE_BoomerAMGSetNumSweeps(solver, 1);          // Sweeeps on each level
        HYPRE_BoomerAMGSetMaxLevels(solver, 35);         // maximum number of levels
        HYPRE_BoomerAMGSetTol(solver, tolerance);        // convergence tolerance
        HYPRE_BoomerAMGSetMaxIter(solver, MaximumCycle); // Maximum number of cycle
        
        // Now setup and solve
        HYPRE_BoomerAMGSetup(solver, AHypre, rhsHypre, solutionHypre);
        HYPRE_BoomerAMGSolve(solver, AHypre, rhsHypre, solutionHypre);
        
        HYPRE_BoomerAMGGetNumIterations(solver, &iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual);
        
        HYPRE_BoomerAMGDestroy(solver);
    }
    else if (solver_id == 1)
    {
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        
        // Set some parameters (See Reference Manual for more parameters)
        HYPRE_PCGSetMaxIter(solver, 1000);     // maximum iterations
        HYPRE_PCGSetTol(solver, tolerance);    // convergence tolerance
        HYPRE_PCGSetTwoNorm(solver, 1);        // use the two norm as the stopping criteria
        
        if (PrintSolutionInformation)
            HYPRE_PCGSetPrintLevel(solver, 2);     // prints out the iteration info
        
        HYPRE_PCGSetLogging(solver, 1);            // needed to get run info later
        
        // Now setup and solve
        HYPRE_ParCSRPCGSetup(solver, AHypre, rhsHypre, solutionHypre);
        HYPRE_ParCSRPCGSolve(solver, AHypre, rhsHypre, solutionHypre);
        
        HYPRE_PCGGetNumIterations(solver, &iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &residual);
        
        HYPRE_ParCSRPCGDestroy(solver);
    }
    else if (solver_id == 2)
    {
        int      sai_max_levels = 1, sai_sym = 1;
        double   sai_threshold = 0.1, sai_filter = 0.05;
        
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        
        // Set some parameters (See Reference Manual for more parameters)
        HYPRE_PCGSetMaxIter(solver, 1000);  // maximum iterations
        HYPRE_PCGSetTol(solver, tolerance); // convergence tolerance
        HYPRE_PCGSetTwoNorm(solver, 1);     // use the two norm as the stopping criteria
        
        if (PrintSolutionInformation)
            HYPRE_PCGSetPrintLevel(solver, 2); // print solve info
        
        HYPRE_PCGSetLogging(solver, 1);        // needed to get run info later
        
        // Now set up the ParaSails preconditioner and specify any parameters
        HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &preconditioner);
        
        // Set some parameters (See Reference Manual for more parameters)
        HYPRE_ParaSailsSetParams(preconditioner, 0.1, 1);
        HYPRE_ParaSailsSetFilter(preconditioner, 0.05);
        HYPRE_ParaSailsSetSym(preconditioner, 1);
        HYPRE_ParaSailsSetLogging(preconditioner, 3);
        
        // Set the PCG preconditioner
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,(HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, preconditioner);
        
        // Now setup and solve!
        HYPRE_ParCSRPCGSetup(solver, AHypre, rhsHypre, solutionHypre);
        HYPRE_ParCSRPCGSolve(solver, AHypre, rhsHypre, solutionHypre);
        
        HYPRE_PCGGetNumIterations(solver, &iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &residual);
        
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_ParaSailsDestroy(preconditioner);
    }
    else if (solver_id == 3)
    {
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        
        // Set some parameters (See Reference Manual for more parameters)
        HYPRE_PCGSetMaxIter(solver, 1000);     // maximum iterations
        HYPRE_PCGSetTol(solver, tolerance);    // convergence tolerance
        HYPRE_PCGSetTwoNorm(solver, 1);        // use the two norm as the stopping criteria
        
        if (PrintSolutionInformation)
            HYPRE_PCGSetPrintLevel(solver, 2); // print solve info
        
        HYPRE_PCGSetLogging(solver, 1);        // needed to get run info later
        
        // Now set up the AMG preconditioner and specify any parameters
        HYPRE_BoomerAMGCreate(&preconditioner);
        
        if (PrintSolutionInformation)
            HYPRE_BoomerAMGSetPrintLevel(preconditioner, 1);     // print amg solutionHypre info
        
        HYPRE_BoomerAMGSetCoarsenType(preconditioner, 6);        // 
        HYPRE_BoomerAMGSetOldDefault(preconditioner);            // 
        HYPRE_BoomerAMGSetRelaxType(preconditioner, 6);          // Sym G.S./Jacobi hybrid
        HYPRE_BoomerAMGSetNumSweeps(preconditioner, 1);          // 
        HYPRE_BoomerAMGSetTol(preconditioner, 0.0);              // convergence tolerance
        HYPRE_BoomerAMGSetMaxIter(preconditioner, 1);            // do only one iteration
        
        // Set the PCG preconditioner
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
        
        // Now setup and solve!
        HYPRE_ParCSRPCGSetup(solver, AHypre, rhsHypre, solutionHypre);
        HYPRE_ParCSRPCGSolve(solver, AHypre, rhsHypre, solutionHypre);
        
        HYPRE_PCGGetNumIterations(solver, &iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &residual);
        
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(preconditioner);
    }
    else if (solver_id == 4)
    {
        int restart = 30, modify = 1;
        
        HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);
        
        // Set some parameters (See Reference Manual for more parameters)
        HYPRE_FlexGMRESSetKDim(solver, restart);     // 
        HYPRE_FlexGMRESSetMaxIter(solver, 1000);     // max iterations
        HYPRE_FlexGMRESSetTol(solver, tolerance);    // convergence tolerance
        
        if (PrintSolutionInformation)
            HYPRE_FlexGMRESSetPrintLevel(solver, 2);     // print solve info
        
        HYPRE_FlexGMRESSetLogging(solver, 1);            // needed to get run info later
        
        // Now set up the AMG preconditioner and specify any parameters
        HYPRE_BoomerAMGCreate(&preconditioner);
        
        if (PrintSolutionInformation)
            HYPRE_BoomerAMGSetPrintLevel(preconditioner, 1);      // print amg solutionHypre info
        
        HYPRE_BoomerAMGSetCoarsenType(preconditioner, 6);         // 
        HYPRE_BoomerAMGSetOldDefault(preconditioner);             // 
        HYPRE_BoomerAMGSetRelaxType(preconditioner, 6);           // Sym G.S./Jacobi hybrid
        HYPRE_BoomerAMGSetNumSweeps(preconditioner, 1);           // 
        HYPRE_BoomerAMGSetTol(preconditioner, 0.0);               // convergence tolerance zero
        HYPRE_BoomerAMGSetMaxIter(preconditioner, 1);             // do only one iteration
        
        // Set the FlexGMRES preconditioner
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, preconditioner);
        
        if (modify)
            HYPRE_FlexGMRESSetModifyPC( solver,(HYPRE_PtrToModifyPCFcn) hypre_FlexGMRESModifyPCAMGExample);
        
        // Now setup and solve
        HYPRE_ParCSRFlexGMRESSetup(solver, AHypre, rhsHypre, solutionHypre);
        HYPRE_ParCSRFlexGMRESSolve(solver, AHypre, rhsHypre, solutionHypre);
        
        HYPRE_FlexGMRESGetNumIterations(solver, &iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &residual);
        
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(preconditioner);
    }
    else
    {
        if (rank == 0)
        {
            std::cout << "Invalid solver_id specified!" << std::endl;
            
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    
    if ( (rank == 0) && (PrintSolutionInformation) )
    {
        std::cout << std::endl << "Iterations = " << iterations << std::endl;
        std::cout << "Final Relative Residual Norm = " << residual << std::endl << std::endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

//_______________________________________________________________________________
// Plot a square matrix, A.
// A : N x N
// 
// This function produces .tec file in the same directory from where this 
// function is called.
//_______________________________________________________________________________
void PlotHypreMatrix ( HYPRE_IJMatrix AHypre, 
                       const int N )
{
    double Color;
    int returnsize = 1, TwoColor;
    
    std::ofstream TecplotWrite("MatrixVisualization.tec", std::ios::out);
    TecplotWrite.flags( std::ios::dec | std::ios::scientific );
    TecplotWrite.precision(8);
    
    if ( !TecplotWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    TecplotWrite << "TITLE = \"Matrix visualization\"\nVariables = \"X\",\"Y\",\"TwoColor\",\"Color\"" << std::endl;
    TecplotWrite << "Zone N = " << 4*N*N << ", E = " << N*N << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
    
    for (int i{}; i < N; i++)
    {
        for (int j{}; j < N; j++)
        {
            HYPRE_IJMatrixGetValues(AHypre, 1, &returnsize, &i, &j, &Color);
            
            TwoColor =  (Color != 0.0 ? 1 : 0);
            
            TecplotWrite << (j-0.5)/N << "\t" << (N-i-0.5)/N << "\t" << TwoColor << "\t" << Color << std::endl;
            TecplotWrite << (j+0.5)/N << "\t" << (N-i-0.5)/N << "\t" << TwoColor << "\t" << Color << std::endl;
            TecplotWrite << (j+0.5)/N << "\t" << (N-i+0.5)/N << "\t" << TwoColor << "\t" << Color << std::endl;
            TecplotWrite << (j-0.5)/N << "\t" << (N-i+0.5)/N << "\t" << TwoColor << "\t" << Color << std::endl;
        }
    }
    
    for (int i{}; i < N*N; i++)
        TecplotWrite << 4*i+1 << "\t" << 4*i+2 << "\t" << 4*i+3 << "\t" << 4*i+4 << std::endl;
    
    TecplotWrite.close();
}
