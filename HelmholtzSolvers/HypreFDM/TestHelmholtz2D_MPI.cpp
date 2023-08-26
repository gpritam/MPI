//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 2D Helmholtz equation.
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
// make TARGET=TestHelmholtz2D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Helmholtz2D_MPI.h"

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
// Boundary condition at the left face (parallel to y axis on which x=xl)
//_______________________________________________________________________________
double LeftBC ( const double y )
{
    double x = xl;
    
    return (-alphal*8.0*PI*cos(8.0*PI*(x-xl)/Lx)/Lx + betal*sin(8.0*PI*(x-xl)/Lx) )*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// Boundary condition at the right face (parallel to y axis on which x=xr)
//_______________________________________________________________________________
double RightBC ( const double y )
{
    double x = xr;
    
    return (alphar*8.0*PI*cos(8.0*PI*(x-xl)/Lx)/Lx + betar*sin(8.0*PI*(x-xl)/Lx) )*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// Boundary condition at the bottom face (parallel to x axis on which y=yl)
// 
// Important : Neumann boundary conditions are always in the outward normal to 
// the boundary.
//_______________________________________________________________________________
double BottomBC ( const double x )
{
    double y = yl;
    
    return (alphab*6.0*PI*sin(6.0*PI*(y-yl)/Ly)/Ly + betab*cos(6.0*PI*(y-yl)/Ly))*sin(8.0*PI*(x-xl)/Lx);
}

//_______________________________________________________________________________
// Boundary condition at the top face (parallel to x axis on which y=yr)
//_______________________________________________________________________________
double TopBC ( const double x )
{
    double y = yr;
    
    return (-alphat*6.0*PI*sin(6.0*PI*(y-yl)/Ly)/Ly + betat*cos(6.0*PI*(y-yl)/Ly))*sin(8.0*PI*(x-xl)/Lx);
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
    CreatePDEPlan(512,512,0.0,1.0,0.0,1.0,false,false);
    
    // Solve first PDE
    double **V;
    
    Allocate(V,Nx+2*offset,Ny+2*offset,offset,offset);
    
    auto start = high_resolution_clock::now();
    
    alphal = 1.0; betal  = 0.0;
    alphar = 0.0; betar  = 1.0;
    
    alphab = 1.0; betab  = 0.0;
    alphat = 0.0; betat  = 1.0;
    
    Lambda = 3.0;
    
    for (int i{}; i < Nx; ++i)
    {
        bottomBC[i] = BottomBC(x[i]);
        topBC[i] = TopBC(x[i]);
    }
    
    for (int j{}; j < Ny; ++j)
    {
        leftBC[j] = LeftBC(y[j]);
        rightBC[j] = RightBC(y[j]);
    }
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            rhsPDE[i][j] = RHS(x[i],y[j]);
    
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
    
    // Deatroy plan    
    Deallocate(V,Nx+2*offset,Ny+2*offset,offset,offset);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
