//________________________________________________________________________________________________
//________________________________________________________________________________________________
// MPI-parallel program for simulating spinoidal decomposition.
// 
// Developed by: Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 17.02.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=SpinoidalDecomposition.cpp
// make run Nproc=4
//________________________________________________________________________________________________

#include "General.h"
#include "General_MPI_2D.h"

using namespace std::chrono;

#define TECPLOT
//#define VISIT

//________________________________________________________________________________________________
// Npx : Number of processors along x-direction
// Npy : Number of processors along y-direction
// 
// Nx : Total number of control volumes in x direction
// Ny : Total number of control volumes in y direction
//________________________________________________________________________________________________

int Nx = 512, Ny = 512, Npx, Npy;

const int SaveInterval = 2000, offset = 2;
const double dt = 1.0E-3, Tf = 250.0;

const double M = 1.0, kappa = 0.5, alpha = 1.0;

double **C, **D, Laplacian, mu, dCdx, dCdy, dx, dy, Lx = Nx, Ly = Ny, t = 0.0;
int TimeStep = 0;

char *s;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SetSubDomains ( int &Npx, 
                     int &Npy, 
                     int &Nx, 
                     int &Ny )
{
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    bool reverse = false;
    double ratio = Lx/Ly;
    
    Npx = size;
    Npy = 1;
    
    if (Lx < Ly)
    {
        reverse = true;
        ratio = Ly/Lx;
    }
    
    int imax = floor(sqrt(double(size)));
    double distance = Absolute(ratio-size);
    
    for (int i{2}; i <= imax; ++i)
    {
        if (size % i == 0)
        {
            if (Absolute(ratio-size/i) < distance)
            {
                distance = Absolute(ratio-size/i);
                Npx = size/i;
                Npy = i;
            }
        }
    }
    
    if (reverse)
        Swap(Npx,Npy);
    
    Nx = ceil(double(Nx)/Npx);
    Ny = ceil(double(Ny)/Npy);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Allocate ()
{
    Allocate(C,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(D,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(C,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(D,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            C[i][j] = 0.4 + 0.02*(0.5 - RandomNumber());
    
    ExchangeData(C,Nx,Ny,offset);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double ComputeEnergy ( double **C )
{
    double localenergy = 0.0, energy;
    
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            dCdx = 0.5*(C[i+1][j] - C[i-1][j])/dx;
            dCdy = 0.5*(C[i][j+1] - C[i][j-1])/dy;
            
            localenergy += alpha*C[i][j]*C[i][j]*(1.0 - C[i][j])*(1.0 - C[i][j]) + 0.5*kappa*(dCdx*dCdx + dCdy*dCdy);
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (size > 1)
        MPI_Allreduce(&localenergy, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    else
        energy = localenergy;
    
    return energy*dx*dy;
}

//________________________________________________________________________________________________
// Main code
//________________________________________________________________________________________________
int main ( int argc, char *argv[] )
{
    MPI_Init(&argc, &argv);
    
    std::cout.flags( std::ios::dec | std::ios::fixed );
    std::cout.precision(8);
    
    auto start = high_resolution_clock::now();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Begin code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SetSubDomains(Npx,Npy,Nx,Ny);
    
    SetMPIEnvironment MPI_ReadyToUse(Nx,Ny,Npx,Npy,offset);
    
    Lx = Npx*Nx;
    Ly = Npy*Ny;
    
    dx = Lx/(Npx*Nx);
    dy = Ly/(Npy*Ny);
    
    Allocate();
    
    Initialize();
    
    #ifdef VISIT
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::out);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        Output << "!NBLOCKS " << size << std::endl;
        
        Output.close();
    }
    #endif
    
    // Update concentration
    while (t < Tf)
    {
        double energy = ComputeEnergy(C);
        
        if (Rank == 0)
            std::cout << "Time step = " << TimeStep << ", Current time = " << t << ", Energy = " << energy << std::endl;
        
        // Step 1: Compute D
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                Laplacian = (C[i+1][j] - 2.0*C[i][j] + C[i-1][j])/(dx*dx) + (C[i][j+1] - 2.0*C[i][j] + C[i][j-1])/(dy*dy);
                
                mu = 2.0*alpha*C[i][j]*(1.0 - C[i][j])*(1.0 - 2.0*C[i][j]);
                
                D[i][j] = mu - kappa*Laplacian;
            }
        }
        
        ExchangeData(D,Nx,Ny,offset);
        
        // Step 2: Update C
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                Laplacian = (D[i+1][j] - 2.0*D[i][j] + D[i-1][j])/(dx*dx) + (D[i][j+1] - 2.0*D[i][j] + D[i][j-1])/(dy*dy);
                
                C[i][j] += dt*M*Laplacian;
                
                if (C[i][j] >= 0.99999) C[i][j] = 0.99999;
                if (C[i][j] <= 0.00001) C[i][j] = 0.00001;
            }
        }
        
        ExchangeData(C,Nx,Ny,offset);
        
        // Step 3: Write file
        if (TimeStep % SaveInterval == 0)
        {
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d-%04d.tec",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios :: out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"C\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    Output << (I*Nx+i)*dx/Lx << "\t" << (J*Ny+j)*dy/Ly << "\t" << C[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%d-%04d.vtk",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios::out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Spinoidal decomposition" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << (I*Nx+i)*dx/Lx << "\t" << (J*Ny+j)*dy/Ly << "\t" << 0.0 << endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS C float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << C[i][j] << std::endl;
            
            Output.close();
            
            if (Rank == 0)
            {
                std::ofstream Output("Output/Field.visit", std::ios::app);
                
                if ( !Output )
                    ErrorMessage("Output file couldnot be opened!");
                
                for (int i = 0; i < size; ++i)
                {
                    sprintf(s,"Field-%d-%04d.vtk",TimeStep/SaveInterval,i);
                    
                    Output << s << std::endl;
                }
                
                Output.close();
            }
            #endif
        }
        
        TimeStep++;
        
        t += dt;
        
        CheckNaN(C,Nx,Ny);
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    Deallocate();
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // End code
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (Rank == 0)
    {
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        
        std::cout << std::endl << "Total time elapsed : " << duration.count() << " seconds." << std::endl << std::endl;
    }
    
    // Finalize the MPI environment.
    MPI_Finalize();
    
    return 0;
}
