//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 1D Helmholtz equation with periodic boundary condition.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 20.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=TestPeriodic1D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Periodic1D_MPI.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x )
{
    return sin(8.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x )
{
    return -(64.0*PI*PI/(Lx*Lx) + Lambda)*sin(8.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void ConstructOperator ( )
{
    for (int i{}; i < (N-2); ++i)
        for (int j{}; j < (N-2); ++j)
            Operator[i][j] = Dxx[i+1][j+1] - (i == j ? Lambda : 0.0);
    
    LUDecomposition(Operator,npivot,determinant,N-2);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
int main ( int argc, char *argv[] )
{
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::cout.flags( std::ios::dec );
    std::cout.precision(8);    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Create plan
    CreatePDEPlan(512,8,0.0,1.0);
    
    // Solve first PDE
    FieldVariable V[Ns];
    
    for (int ii{}; ii < Ns; ++ii)
        Allocate(V[ii].phi,N);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int ii{}; ii < Ns; ++ii)
        for (int i{}; i < N; ++i)
            rhsPDE[ii][i] = RHS(x[ii*(N-1)+i]);
    
    SolvePDE(V,ConstructOperator);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (Rank == 0)
    {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << std::endl << "Total time elapsed : " << duration.count() << " microseconds." << std::endl << std::endl;
    }
    
    WriteFile(V,ExactFunction);
    
    ComputeError(V,ExactFunction);
    
    // Deatroy plan    
    for (int ii{}; ii < Ns; ++ii)
        Deallocate(V[ii].phi,N);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
