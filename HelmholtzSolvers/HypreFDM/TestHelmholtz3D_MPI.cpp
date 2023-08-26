//_______________________________________________________________________________
//_______________________________________________________________________________
// Test MPI-parallel Program for 3D Helmholtz equation.
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
// make TARGET=TestHelmholtz3D_MPI.cpp
// make run Nproc=4
//_______________________________________________________________________________

#include "Helmholtz3D_MPI.h"

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
// Boundary condition at the left face (y-z plane on which x=xl)
//_______________________________________________________________________________
double LeftBC ( const double y, 
                const double z )
{
    double x = xl;
    
    return (-alphal*8.0*PI*cos(8.0*PI*(x-xl)/Lx)/Lx + betal*sin(8.0*PI*(x-xl)/Lx) )*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Boundary condition at the right face (y-z plane on which x=xr)
//_______________________________________________________________________________
double RightBC ( const double y, 
                 const double z )
{
    double x = xr;
    
    return (alphar*8.0*PI*cos(8.0*PI*(x-xl)/Lx)/Lx + betar*sin(8.0*PI*(x-xl)/Lx) )*cos(6.0*PI*(y-yl)/Ly)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Boundary condition at the bottom face (x-z plane on which y=yl)
//_______________________________________________________________________________
double BottomBC ( const double x, 
                  const double z )
{
    double y = yl;
    
    return (alphab*6.0*PI*sin(6.0*PI*(y-yl)/Ly)/Ly + betab*cos(6.0*PI*(y-yl)/Ly) )*sin(8.0*PI*(x-xl)/Lx)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Boundary condition at the top face (x-z plane on which y=yr)
//_______________________________________________________________________________
double TopBC ( const double x, 
               const double z )
{
    double y = yr;
    
    return (-alphat*6.0*PI*sin(6.0*PI*(y-yl)/Ly)/Ly + betat*cos(6.0*PI*(y-yl)/Ly) )*sin(8.0*PI*(x-xl)/Lx)*cos(4.0*PI*(z-zl)/Lz);
}

//_______________________________________________________________________________
// Boundary condition at the back face (x-y plane on which z=zl)
// 
// Important : Neumann boundary conditions are always in the outward normal to 
// the boundary.
//_______________________________________________________________________________
double BackBC ( const double x, 
                const double y )
{
    double z = zl;
    
    return (alphaback*4.0*PI*sin(4.0*PI*(z-zl)/Lz)/Lz + betaback*cos(4.0*PI*(z-zl)/Lz))*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
}

//_______________________________________________________________________________
// Boundary condition at the front face (x-y plane on which z=zr)
//_______________________________________________________________________________
double FrontBC ( const double x, 
                 const double y )
{
    double z = zr;
    
    return (-alphaf*4.0*PI*sin(4.0*PI*(z-zl)/Lz)/Lz + betaf*cos(4.0*PI*(z-zl)/Lz))*sin(8.0*PI*(x-xl)/Lx)*cos(6.0*PI*(y-yl)/Ly);
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
    CreatePDEPlan(64,32,250,0.0,1.0,0.0,1.0,0.0,1.0,false,false,false);
    
    // Solve first PDE
    double ***V;
    
    Allocate(V,Nx+2*offset,Ny+2*offset,Nz+2*offset,offset,offset,offset);
    
    auto start = high_resolution_clock::now();
    
    alphal = 0.0; betal  = 1.0;
    alphar = 0.0; betar  = 1.0;
    
    alphab = 0.0; betab  = 1.0;
    alphat = 0.0; betat  = 1.0;
    
    alphaback = 1.0; betaback = 0.0;
    alphaf = 1.0;      betaf    = 0.0;
    
    Lambda = 2.5;
    
    for (int j{}; j < Ny; ++j)
    {
        for (int k{}; k < Nz; ++k)
        {
            leftBC[j][k] = LeftBC(y[j],z[k]);
            rightBC[j][k] = RightBC(y[j],z[k]);
        }
    }
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            backBC[i][j] = BackBC(x[i],y[j]);
            frontBC[i][j] = FrontBC(x[i],y[j]);
        }
    }
    
    for (int i{}; i < Nx; ++i)
    {
        for (int k{}; k < Nz; ++k)
        {
            bottomBC[i][k] = BottomBC(x[i],z[k]);
            topBC[i][k] = TopBC(x[i],z[k]);
        }
    }
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                rhsPDE[i][j][k] = RHS(x[i],y[j],z[k]);
    
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
    Deallocate(V,Nx+2*offset,Ny+2*offset,Nz+2*offset,offset,offset,offset);
    
    DestroyPDEPlan();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
