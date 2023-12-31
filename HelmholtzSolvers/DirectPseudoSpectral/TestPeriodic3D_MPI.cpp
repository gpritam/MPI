//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 3D Helmholtz equation with periodic boundary condition.
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
// make TARGET=TestPeriodic3D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Periodic3D_MPI.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x, 
                       const double y, 
                       const double z )
{
    return sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x, 
             const double y, 
             const double z )
{
    return -(64.0*PI*PI/(Lx*Lx) + 36.0*PI*PI/(Ly*Ly) + 16.0*PI*PI/(Lz*Lz) + Lambda)*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void ConstructOperator ( const int iw, 
                         const int jw )
{
    for (int i{}; i < (N-2); ++i)
        for (int j{}; j < (N-2); ++j)
            Operator[i][j] = Dzz[i+1][j+1] - (i == j ? (Lambda + 4.0*PI*PI*iw*iw/(Lx*Lx) + 4.0*PI*PI*jw*jw/(Ly*Ly)) : 0.0);
    
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
    
    std::cout.flags(std::ios::dec);
    std::cout.precision(8);    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Create plan
    CreatePDEPlan(64,32,250,10,0.0,1.0,0.0,1.0,0.0,1.0);
    CreateFFTPlan(Maximum(Nx,Ny));
    
    // Solve first PDE
    FieldVariable V[Ns];
    
    for (int kk{}; kk < Ns; ++kk)
        Allocate(V[kk].phi,2*Nx,NY,N);
    
    auto start = high_resolution_clock::now();
    
    Lambda = 1.0;
    
    for (int kk{}; kk < Ns; ++kk)
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                for (int k{}; k < N; ++k)
                    rhsPDE[kk][i][j][k] = RHS(xl+i*dx,yl+j*dy,z[kk*(N-1)+k]);
    
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
    
    // Destroy plan    
    for (int kk{}; kk < Ns; ++kk)
        Deallocate(V[kk].phi,2*Nx,NY,N);
    
    DestroyFFTPlan(Maximum(Nx,Ny));
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
