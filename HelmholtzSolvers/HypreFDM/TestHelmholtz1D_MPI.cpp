//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 1D Helmholtz equation.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 15.08.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//_______________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=TestHelmholtz1D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Helmholtz1D_MPI.h"

using namespace std::chrono;

//_______________________________________________________________________________
// Exact solution
//_______________________________________________________________________________
double ExactFunction ( const double x )
{
    return sin(2.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// Boundary condition at the left face
// 
// Important : Neumann boundary conditions are always in the outward normal to 
// the boundary.
//_______________________________________________________________________________
double LeftBC ( )
{
    double x = xl;
    
    return (-alphal*2.0*PI*cos(2.0*PI*(x-xl)/Lx)/Lx + betal*sin(2.0*PI*(x-xl)/Lx));
}

//_______________________________________________________________________________
// Boundary condition at the right face
//_______________________________________________________________________________
double RightBC ( )
{
    double x = xr;
    
    return (alphar*2.0*PI*cos(2.0*PI*(x-xl)/Lx)/Lx + betar*sin(2.0*PI*(x-xl)/Lx));
}

//_______________________________________________________________________________
// Define right hand side
//_______________________________________________________________________________
double RHS ( const double x )
{
    return -(4.0*PI*PI/(Lx*Lx) + Lambda)*sin(2.0*PI*(x-xl)/Lx);
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
    CreatePDEPlan(512,0.0,1.0,false);
    
    // Solve first PDE
    double *V, *W;
    
    Allocate(V,Nx+2*offset,offset);
    Allocate(W,Nx+2*offset,offset);
    
    auto start = high_resolution_clock::now();
    
    alphal = 1.0; betal  = 0.0;
    alphar = 0.0; betar  = 1.0;
    Lambda = 3.0;
    
    leftBC = LeftBC();
    rightBC = RightBC();
    
    for (int i{}; i < Nx; ++i)
        rhsPDE[i] = RHS(x[i]);
    
    SolvePDE(V);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (Rank == 0)
    {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        std::cout << std::endl << "Total time elapsed : " << duration.count() << " milliseconds." << std::endl << std::endl;
    }
    
    WriteFile(V,ExactFunction);
    ComputeError(V,ExactFunction);
    
    // Solve second PDE
    alphal = 1.0; betal  = 1.0;
    alphar = 0.0; betar  = 1.0;
    Lambda = 5.0;
    
    leftBC = LeftBC();
    rightBC = RightBC();
    
    for (int i{}; i < Nx; ++i)
        rhsPDE[i] = RHS(x[i]);
    
    SolvePDE(W);
    
    WriteFile(W,ExactFunction);
    ComputeError(W,ExactFunction);
    
    // Deatroy plan    
    Deallocate(V,Nx+2*offset,offset);
    Deallocate(W,Nx+2*offset,offset);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
