//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 2D Helmholtz equation with periodic boundary condition.
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
// make TARGET=TestPeriodic2D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Periodic2D_MPI.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x, 
                       const double y )
{
    return sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x, 
             const double y )
{
    return -(64.0*PI*PI/(Lx*Lx) + 36.0*PI*PI/(Ly*Ly) + Lambda)*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void ConstructOperator ( const int iw )
{
    for (int i{}; i < (N-2); ++i)
        for (int j{}; j < (N-2); ++j)
            Operator[i][j] = Dyy[i+1][j+1] - (i == j ? (Lambda + 4.0*PI*PI*iw*iw/(Lx*Lx)) : 0.0);
    
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
    CreatePDEPlan(512,512,10,0.0,1.0,0.0,1.0);
    CreateFFTPlan(Nx);
    
    // Solve first PDE
    FieldVariable V[Ns];
    
    for (int jj{}; jj < Ns; ++jj)
        Allocate(V[jj].phi,NX,N);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int jj{}; jj < Ns; ++jj)
        for (int j{}; j < N; ++j)
            for (int i{}; i < Nx; ++i)
                rhsPDE[jj][i][j] = RHS((xl+i*dx),y[jj*(N-1)+j]);
    
    SolvePDE(V,ConstructOperator);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (Rank == 0)
    {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    }
    
    WriteFile(V,ExactFunction);
    
    ComputeError(V,ExactFunction);
    
    // Deatroy plan    
    for (int jj{}; jj < Ns; ++jj)
        Deallocate(V[jj].phi,NX,N);
    
    DestroyFFTPlan(Nx);
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
