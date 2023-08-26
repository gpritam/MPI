//________________________________________________________________________________________________
//________________________________________________________________________________________________
// MPI-parallel program for simulating dendritic solidification.
// 
// Developed by: Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 31.01.2023
// Bangalore
//________________________________________________________________________________________________
//________________________________________________________________________________________________

//________________________________________________________________________________________________
// To run this code, issue the following commands
// 
// make TARGET=DendriticGrowth.cpp
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

const int SaveInterval = 50, offset = 2;
const double Lx = 9.0, Ly = 9.0;
const double xmin = 0.0, xmax = xmin + Lx, ymin = 0.0, ymax = ymin + Ly;
const double dt = 5.0E-5, Tf = 0.4;

const double tau = 3.0E-4, epsilonb = 1.0E-2, kappa = 1.8, delta = 2.0E-2;
const double beta = 6.0, alpha = 0.9, Gamma = 10.0, Teq = 1.0, Theta0 = 0.2, r0 = 5.0*Lx*Lx/250000.0;

double **phi, **phiold, **T, **Told, **Epsilon, **Epsilon_theta, dx, dy;
double Laplacian, Theta, phix, phiy, term0, term1, term2, term3, phixp, phixm, phiyp, phiym, m, t = 0.0;
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
    Allocate(phi,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(T,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(phiold,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Told,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Epsilon,Nx+2*offset,Ny+2*offset,offset,offset);
    Allocate(Epsilon_theta,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Allocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Deallocate ()
{
    Deallocate(phi,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(T,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(phiold,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Told,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Epsilon,Nx+2*offset,Ny+2*offset,offset,offset);
    Deallocate(Epsilon_theta,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Deallocate(s,200);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void Initialize ()
{
    double x, y, r;
    
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            x = (I*Nx+i)*dx;
            y = (J*Ny+j)*dy;
            
            r = (x-Lx/2.0)*(x-Lx/2.0) + (y-Ly/2.0)*(y-Ly/2.0);
            
            phi[i][j] = (r < r0 ? 1.0 : 0.0);
            
            T[i][j] = 0.0;
        }
    }
    
    ExchangeData(phi,Nx,Ny,offset);
    ExchangeData(T,Nx,Ny,offset);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
double PhiX ( const int i, 
              const int j )
{
    return 0.5*(phiold[i+1][j]-phiold[i-1][j])/dx;
}

//________________________________________________________________________________________________
// 
//_______________________________________________________________________________________________
double PhiY ( const int i, 
              const int j )
{
    return 0.5*(phiold[i][j+1]-phiold[i][j-1])/dy;
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
        if (Rank == 0)
            std::cout << "Time step = " << TimeStep << ", Current time = " << t << std::endl;
        
        // Step 1
        for (int i = -offset; i < Nx+offset; ++i)
        {
            for (int j = -offset; j < Ny+offset; ++j)
            {
                phiold[i][j] = phi[i][j];
                Told[i][j] = T[i][j];
            }
        }
        
        // Step 2: Compute Epsilon and (d Epsilon/d theta)
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                phix = 0.5*(phiold[i+1][j] - phiold[i-1][j])/dx;
                phiy = 0.5*(phiold[i][j+1] - phiold[i][j-1])/dy;
                
                Theta = atan2(phiy,phix);
                
                Epsilon[i][j] = epsilonb + epsilonb*delta*cos(beta*(Theta-Theta0));
                Epsilon_theta[i][j] = -epsilonb*beta*delta*sin(beta*(Theta-Theta0));
            }
        }
        
        ExchangeData(Epsilon,Nx,Ny,offset);
        ExchangeData(Epsilon_theta,Nx,Ny,offset);
        
        // Step 3: Update phase-field parameter and temperature
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                // Update phase-field parameter
                Laplacian = (phiold[i+1][j] - 2.0*phiold[i][j] + phiold[i-1][j])/(dx*dx) + (phiold[i][j+1] - 2.0*phiold[i][j] + phiold[i][j-1])/(dy*dy);
                m = alpha*atan(Gamma*(Teq-Told[i][j]))/PI;
                
                term0 = 0.5*(Epsilon[i][j+1]*Epsilon_theta[i][j+1]*PhiX(i,j+1) - Epsilon[i][j-1]*Epsilon_theta[i][j-1]*PhiX(i,j-1))/dy;
                term1 = 0.5*(Epsilon[i+1][j]*Epsilon_theta[i+1][j]*PhiY(i+1,j) - Epsilon[i-1][j]*Epsilon_theta[i-1][j]*PhiY(i-1,j))/dx;
                term2 = Epsilon[i][j]*Epsilon[i][j]*Laplacian;
                term3 = phiold[i][j]*(1.0-phiold[i][j])*(phiold[i][j]-0.5+m);
                
                phi[i][j] += dt*(term0 - term1 + term2 + term3)/tau;
                
                // Update temperature
                Laplacian = (Told[i+1][j] - 2.0*Told[i][j] + Told[i-1][j])/(dx*dx) + (Told[i][j+1] - 2.0*Told[i][j] + Told[i][j-1])/(dy*dy);
                
                T[i][j] += dt*Laplacian + kappa*(phi[i][j]-phiold[i][j]);
            }
        }
        
        ExchangeData(phi,Nx,Ny,offset);
        ExchangeData(T,Nx,Ny,offset);
        
        // Step 4: Write file
        if (TimeStep % SaveInterval == 0)
        {
            #ifdef TECPLOT
            sprintf(s,"Output/Field-%d-%04d.tec",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios :: out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "TITLE = Flow" << std::endl << "VARIABLES = \"X\", \"Y\", \"phi\", \"T\" " << std::endl;
            Output << "Zone T = U I = " << Ny << " J = " << Nx << std::endl;
            
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Ny; ++j)
                    Output << (I*Nx+i)*dx << "\t" << (J*Ny+j)*dy << "\t" << phi[i][j] << "\t" << T[i][j] << std::endl;
            
            Output.close();
            #endif
            
            #ifdef VISIT
            sprintf(s,"Output/Field-%d-%04d.vtk",TimeStep/SaveInterval,Rank);
            
            std::ofstream Output( s, std::ios :: out );
            Output.flags( std::ios::dec );
            Output.precision(10);
            
            if ( !Output )
                ErrorMessage("Output file couldnot be opened!");
            
            Output << "# vtk DataFile Version 3.1" << std::endl;
            Output << "Dendritic growth" << std::endl;
            Output << "ASCII" << std::endl;
            Output << "DATASET STRUCTURED_GRID" << std::endl;
            Output << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << std::endl; 
            Output << "POINTS " << Nx*Ny << " FLOAT" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << (I*Nx+i)*dx << "\t" << (J*Ny+j)*dy << "\t" << 0.0 << std::endl;
            
            Output << std::endl << "POINT_DATA " << Nx*Ny << std::endl;
            Output << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    Output << phi[i][j] << std::endl;
            
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
        
        CheckNaN(phi,Nx,Ny);
        
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
