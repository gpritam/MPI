//_______________________________________________________________________________
//_______________________________________________________________________________
// MPI-parallel Program for 2D Helmholtz equation (Domain decomposition along y direction).
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 20.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

#include "NonPeriodic2D_MPI.h"

struct PDE
{
    double **ul, **ur;
    
    double *zetal, *zetar, *b, *c, *rhs;
};

PDE *S;

double *br, *rhsr;

int Rank, size;

int Nx, Ny, NX, N, Ns;

double xl, xr, yl, yr, Lx, Ly, dx;

double alphal, betal, alphar, betar, Lambda;

double **Dy, **Dyy, **Operator, *y, ***rhsPDE, *leftBC, *rightBC, determinant;

int *npivot;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const int Ny, 
                     const int N, 
                     const double xl, 
                     const double xr, 
                     const double yl, 
                     const double yr )
{
    ::Nx = Nx;
    ::Ny = Ny;
    ::N  = N;
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    
    NX = 2*Nx;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    dx = Lx/Nx;
    
    Ns = ceil(((double)(Ny))/(size*N));
    
    S = new (std::nothrow) PDE [Ns];
    
    if (S == nullptr)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int jj{}; jj < Ns; ++jj)
    {
        Allocate(S[jj].ul,NX,N);
        Allocate(S[jj].ur,NX,N);
        
        Allocate(S[jj].zetal,NX);
        Allocate(S[jj].zetar,NX);
        Allocate(S[jj].b,NX);
        Allocate(S[jj].c,NX);
        Allocate(S[jj].rhs,NX);
    }
    
    Allocate(y,Ns*(N-1)+1);
    Allocate(rhsPDE,Ns,NX,N);
    
    Allocate(br,NX);
    Allocate(rhsr,NX);
    Allocate(leftBC,NX);
    Allocate(rightBC,NX);
    
    Allocate(Dy,N,N);
    Allocate(Dyy,N,N);
    
    Allocate(Operator,N-2,N-2);
    Allocate(npivot,N-2);
    
    double Lyp = Ly/size;
    
    double *ylocal;
    
    Allocate(ylocal,N);
    
    GaussLobattoLegendrePoints(ylocal,N);
    
    for (int jj{}; jj < Ns; ++jj)
    {
        double yl, yr;
        
        yl = Rank*Lyp + jj*Lyp/Ns;
        yr = yl + Lyp/Ns;
        
        for (int j{}; j < N; ++j)
            y[jj*(N-1)+j] = 0.5*(yl + yr + (Lyp/Ns)*ylocal[j]);
    }
    
    FirstDerivativeMatrix(Dy,ylocal,N);
    SecondDerivativeMatrix(Dyy,Dy,N);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dy[i][j] *= ((2.0*Ns)/Lyp);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dyy[i][j] *= ((4.0*Ns*Ns)/(Lyp*Lyp));
    
    Deallocate(ylocal,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    for (int jj{}; jj < Ns; ++jj)
    {
        Deallocate(S[jj].ul,NX,N);
        Deallocate(S[jj].ur,NX,N);
        
        Deallocate(S[jj].zetal,NX);
        Deallocate(S[jj].zetar,NX);
        Deallocate(S[jj].b,NX);
        Deallocate(S[jj].c,NX);
        Deallocate(S[jj].rhs,NX);
    }
    
    delete [] S;
    
    S = nullptr;
    
    Deallocate(y,Ns*(N-1)+1);
    Deallocate(rhsPDE,Ns,NX,N);
    
    Deallocate(br,NX);
    Deallocate(rhsr,NX);
    Deallocate(leftBC,NX);
    Deallocate(rightBC,NX);
    
    Deallocate(Dy,N,N);
    Deallocate(Dyy,N,N);
    
    Deallocate(Operator,N-2,N-2);
    Deallocate(npivot,N-2);
}

//________________________________________________________________________________________________
// This program returns first derivative on 'k'th Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double FirstDerivative ( int k, 
                         double *h )
{
    return FirstDerivative(k,h,Dy,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void ConstructFundamentalSolutions ( FieldVariable U[], 
                                     const int i )
{
    for (int jj{}; jj < Ns; ++jj)
    {
        U[jj].phi[i][0]   = 0.0;
        U[jj].phi[i][N-1] = 0.0;
        
        U[jj].phi[i]++;
        
        for (int j{}; j < (N-2); ++j)
            U[jj].phi[i][j] = rhsPDE[jj][i][j+1];
        
        LUSolve(Operator,npivot,U[jj].phi[i],N-2);
        
        U[jj].phi[i]--;
    }
    
    S[0].ul[i][0]   = 1.0;
    S[0].ul[i][N-1] = 0.0;
    S[0].ur[i][0]   = 0.0;
    S[0].ur[i][N-1] = 1.0;
    
    S[0].ul[i]++;
    S[0].ur[i]++;
    
    for (int j{}; j < (N-2); ++j)
        S[0].ul[i][j] = -Dyy[j+1][0];
    
    LUSolve(Operator,npivot,S[0].ul[i],N-2);
    
    for (int j{}; j < (N-2); ++j)
        S[0].ur[i][j] = -Dyy[j+1][N-1];
    
    LUSolve(Operator,npivot,S[0].ur[i],N-2);
    
    S[0].ul[i]--;
    S[0].ur[i]--;
    
    for (int jj = 1; jj < Ns; ++jj)
    {
        for (int j{}; j < N; ++j)
        {
            S[jj].ul[i][j] = S[0].ul[i][j];
            S[jj].ur[i][j] = S[0].ur[i][j];
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( FieldVariable U[], 
                void (*ConstructOperator)(const int) )
{
    // Step 1: Take forward Fourier transform of the right hand side and BCs
    for (int jj{}; jj < Ns; ++jj)
        FourierTransformX(rhsPDE[jj],1,Nx,N);
    
    FourierTransformX(leftBC,1,Nx);
    FourierTransformX(rightBC,1,Nx);
    
    // Step 2: Construct fundamental solutions on every smallest sub-domains
    for (int iw{}; iw <= Nx/2; iw++)
    {
        // Construct 'Operator' matrix for wave numbers |iw| and compute its LU decomposition.
        ConstructOperator(iw);
        
        // Solve ODEs with this operator
        if (iw == 0)
        {
            ConstructFundamentalSolutions(U,0);
            ConstructFundamentalSolutions(U,1);
        }
        else if (iw == Nx/2)
        {
            ConstructFundamentalSolutions(U,Nx);
            ConstructFundamentalSolutions(U,Nx+1);
        }
        else
        {
            ConstructFundamentalSolutions(U,2*iw);
            ConstructFundamentalSolutions(U,2*iw+1);
            ConstructFundamentalSolutions(U,2*(Nx-iw));
            ConstructFundamentalSolutions(U,2*(Nx-iw)+1);
        }
    }
    
    // Step 3: TDMA forward sweep
    if (size == 1)
    {
        for (int i{}; i < NX; ++i)
        {
            S[0].b[i] = -alphal*FirstDerivative(0,S[0].ul[i]) + betal;
            S[0].c[i] = -alphal*FirstDerivative(0,S[0].ur[i]);
            S[0].rhs[i] = leftBC[i] + alphal*FirstDerivative(0,U[0].phi[i]);
        }
        
        for (int jj = 1; jj < Ns; ++jj)
        {
            for (int i{}; i < NX; ++i)
                leftBC[i] = FirstDerivative(N-1,S[jj-1].ul[i]);
            
            for (int i{}; i < NX; ++i)
            {
                S[jj].b[i] = FirstDerivative(N-1,S[jj-1].ur[i]) - FirstDerivative(0,S[jj].ul[i]) - (S[jj-1].c[i]/S[jj-1].b[i])*leftBC[i];
                S[jj].c[i] = -FirstDerivative(0,S[jj].ur[i]);
                S[jj].rhs[i] = FirstDerivative(0,U[jj].phi[i]) - FirstDerivative(N-1,U[jj-1].phi[i]) - (S[jj-1].rhs[i]/S[jj-1].b[i])*leftBC[i];
            }
        }
        
        for (int i{}; i < NX; ++i)
            leftBC[i] = alphar*FirstDerivative(N-1,S[Ns-1].ul[i]);
        
        for (int i{}; i < NX; ++i)
        {
            br[i] = alphar*FirstDerivative(N-1,S[Ns-1].ur[i]) + betar - (S[Ns-1].c[i]/S[Ns-1].b[i])*leftBC[i];
            rhsr[i] = rightBC[i] - alphar*FirstDerivative(N-1,U[Ns-1].phi[i]) - (S[Ns-1].rhs[i]/S[Ns-1].b[i])*leftBC[i];
        }
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == 0)
            {
                for (int i{}; i < NX; ++i)
                {
                    S[0].b[i] = -alphal*FirstDerivative(0,S[0].ul[i]) + betal;
                    S[0].c[i] = -alphal*FirstDerivative(0,S[0].ur[i]);
                    S[0].rhs[i] = leftBC[i] + alphal*FirstDerivative(0,U[0].phi[i]);
                }
            }
            else
            {
                MPI_Recv(&(S[0].zetal[0]), NX, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar[0]), NX, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int i{}; i < NX; ++i)
                {
                    S[0].b[i] =  -FirstDerivative(0,S[0].ul[i]) + S[0].zetal[i];
                    S[0].c[i] =  -FirstDerivative(0,S[0].ur[i]);
                    S[0].rhs[i] = FirstDerivative(0,U[0].phi[i]) + S[0].zetar[i];
                }
            }
            
            if (Rank == (size-1))
            {
                for (int i{}; i < NX; ++i)
                    leftBC[i] = alphar*FirstDerivative(N-1,S[0].ul[i]);
                
                for (int i{}; i < NX; ++i)
                {
                    br[i] = alphar*FirstDerivative(N-1,S[0].ur[i]) + betar - (S[0].c[i]/S[0].b[i])*leftBC[i];
                    rhsr[i] = rightBC[i] - alphar*FirstDerivative(N-1,U[0].phi[i]) - (S[0].rhs[i]/S[0].b[i])*leftBC[i];
                }
            }
            else
            {
                for (int i{}; i < NX; ++i)
                    leftBC[i] = FirstDerivative(N-1,S[0].ul[i]);
                
                for (int i{}; i < NX; ++i)
                {
                    S[0].zetal[i] = FirstDerivative(N-1,S[0].ur[i]) - (S[0].c[i]/S[0].b[i])*leftBC[i];
                    S[0].zetar[i] = -FirstDerivative(N-1,U[0].phi[i]) - (S[0].rhs[i]/S[0].b[i])*leftBC[i];
                }
                
                MPI_Send(&(S[0].zetal[0]), NX, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[0].zetar[0]), NX, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
        else
        {
            if (Rank == 0)
            {
                for (int i{}; i < NX; ++i)
                {
                    S[0].b[i] = -alphal*FirstDerivative(0,S[0].ul[i]) + betal;
                    S[0].c[i] = -alphal*FirstDerivative(0,S[0].ur[i]);
                    S[0].rhs[i] = leftBC[i] + alphal*FirstDerivative(0,U[0].phi[i]);
                }
            }
            else
            {
                MPI_Recv(&(S[0].zetal[0]), NX, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar[0]), NX, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int i{}; i < NX; ++i)
                {
                    S[0].b[i] =  -FirstDerivative(0,S[0].ul[i]) + S[0].zetal[i];
                    S[0].c[i] =  -FirstDerivative(0,S[0].ur[i]);
                    S[0].rhs[i] = FirstDerivative(0,U[0].phi[i]) + S[0].zetar[i];
                }
            }
            
            for (int jj = 1; jj < Ns; ++jj)
            {
                for (int i{}; i < NX; ++i)
                    leftBC[i] = FirstDerivative(N-1,S[jj-1].ul[i]);
                
                for (int i{}; i < NX; ++i)
                {
                    S[jj].b[i] = FirstDerivative(N-1,S[jj-1].ur[i]) - FirstDerivative(0,S[jj].ul[i]) - (S[jj-1].c[i]/S[jj-1].b[i])*leftBC[i];
                    S[jj].c[i] = -FirstDerivative(0,S[jj].ur[i]);
                    S[jj].rhs[i] = FirstDerivative(0,U[jj].phi[i]) - FirstDerivative(N-1,U[jj-1].phi[i]) - (S[jj-1].rhs[i]/S[jj-1].b[i])*leftBC[i];
                }
            }
            
            if (Rank == (size-1))
            {
                for (int i{}; i < NX; ++i)
                    leftBC[i] = alphar*FirstDerivative(N-1,S[Ns-1].ul[i]);
                
                for (int i{}; i < NX; ++i)
                {
                    br[i] = alphar*FirstDerivative(N-1,S[Ns-1].ur[i]) + betar - (S[Ns-1].c[i]/S[Ns-1].b[i])*leftBC[i];
                    rhsr[i] = rightBC[i] - alphar*FirstDerivative(N-1,U[Ns-1].phi[i]) - (S[Ns-1].rhs[i]/S[Ns-1].b[i])*leftBC[i];
                }
            }
            else
            {
                for (int i{}; i < NX; ++i)
                    leftBC[i] = FirstDerivative(N-1,S[Ns-1].ul[i]);
                
                for (int i{}; i < NX; ++i)
                {
                    S[Ns-1].zetal[i] = FirstDerivative(N-1,S[Ns-1].ur[i]) - (S[Ns-1].c[i]/S[Ns-1].b[i])*leftBC[i];
                    S[Ns-1].zetar[i] = -FirstDerivative(N-1,U[Ns-1].phi[i]) - (S[Ns-1].rhs[i]/S[Ns-1].b[i])*leftBC[i];
                }
                
                MPI_Send(&(S[Ns-1].zetal[0]), NX, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[Ns-1].zetar[0]), NX, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
    }
    
    // Step 4: TDMA backward sweep
    if (size == 1)
    {
        for (int jj = (Ns-1); jj >= 0; --jj)
        {
            for (int i{}; i < NX; ++i)
            {
                S[jj].zetar[i] = (jj == (Ns-1) ? rhsr[i]/br[i] : S[jj+1].zetal[i]);
                S[jj].zetal[i] = (S[jj].rhs[i] - S[jj].c[i] * S[jj].zetar[i])/S[jj].b[i];
            }
        }
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == (size-1))
                for (int i{}; i < NX; ++i)
                    S[0].zetar[i] = rhsr[i]/br[i];
            else
                MPI_Recv(&(S[0].zetar[0]), NX, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i{}; i < NX; ++i)
                S[0].zetal[i] = (S[0].rhs[i] - S[0].c[i] * S[0].zetar[i])/S[0].b[i];
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal[0]), NX, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
        else
        {
            if (Rank == (size-1))
                for (int i{}; i < NX; ++i)
                    S[Ns-1].zetar[i] = rhsr[i]/br[i];
            else
                MPI_Recv(&(S[Ns-1].zetar[0]), NX, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i{}; i < NX; ++i)
                S[Ns-1].zetal[i] = (S[Ns-1].rhs[i] - S[Ns-1].c[i] * S[Ns-1].zetar[i])/S[Ns-1].b[i];
            
            for (int jj = (Ns-2); jj >= 0; --jj)
            {
                for (int i{}; i < NX; ++i)
                {
                    S[jj].zetar[i] = S[jj+1].zetal[i];
                    S[jj].zetal[i] = (S[jj].rhs[i] - S[jj].c[i] * S[jj].zetar[i])/S[jj].b[i];
                }
            }
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal[0]), NX, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
    }
    
    // Step 5: Construct final solution
    for (int jj{}; jj < Ns; ++jj)
        for (int i{}; i < NX; ++i)
            for (int j{}; j < N; ++j)
                U[jj].phi[i][j] += S[jj].zetal[i]*S[jj].ul[i][j] + S[jj].zetar[i]*S[jj].ur[i][j];
    
    for (int jj{}; jj < Ns; ++jj)
        FourierTransformX(U[jj].phi,-1,Nx,N);
}

//_______________________________________________________________________________
// First derivative with respect to x
//_______________________________________________________________________________
void FirstDerivativeX ( FieldVariable Ux[], 
                        FieldVariable U[] )
{
    
}

//_______________________________________________________________________________
// Second derivative with respect to x
//_______________________________________________________________________________
void SecondDerivativeX ( FieldVariable Uxx[], 
                         FieldVariable U[] )
{
    
}

//_______________________________________________________________________________
// First derivative with respect to y
//_______________________________________________________________________________
void FirstDerivativeY ( FieldVariable Uy[], 
                        FieldVariable U[] )
{
    double ay;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int jj{}; jj < Ns; ++jj)
        {
            for (int j{}; j < N; ++j)
            {
                ay = 0.0;
                
                for (int k{}; k < N; ++k)
                    ay += Dy[j][k]*U[jj].phi[i][k];
                
                Uy[jj].phi[i][j] = ay;
            }
        }
    }
}

//_______________________________________________________________________________
// Second derivative with respect to y
//_______________________________________________________________________________
void SecondDerivativeY ( FieldVariable Uyy[], 
                         FieldVariable U[] )
{
    double ayy;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int jj{}; jj < Ns; ++jj)
        {
            for (int j{}; j < N; ++j)
            {
                ayy = 0.0;
                
                for (int k{}; k < N; ++k)
                    ayy += Dyy[j][k]*U[jj].phi[i][k];
                
                Uyy[jj].phi[i][j] = ayy;
            }
        }
    }
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( FieldVariable U[], 
                    double (*ExactFunction)(const double, const double) )
{
    double Error, LocalError = 0.0;
    
    for (int jj{}; jj < Ns; ++jj)
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < N; ++j)
                LocalError = Maximum(LocalError,Absolute(U[jj].phi[i][j]-ExactFunction((xl+i*dx),y[jj*(N-1)+j])));
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (size > 1)
        MPI_Allreduce(&LocalError, &Error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    else
        Error = LocalError;
    
    if (Rank == 0)
        std::cout << "Error = " << Error << std::endl << std::endl;
}

//_______________________________________________________________________________
// Write file in Tecplot format
//_______________________________________________________________________________
void WriteFile ( FieldVariable U[], 
                 double (*ExactFunction)(const double, const double) )
{
    char *s;
    
    Allocate(s,200);
    
    #ifdef TECPLOT
    sprintf(s,"Output/Solution-%04d.tec",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "TITLE = \"Helmholtz solver\"" << std::endl;
    FileWrite << "Variables = \"X\",\"Y\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nx << ", J = " << Ns*N << ", DATAPACKING=POINT" << std::endl;
    
    for (int jj{}; jj < Ns; ++jj)
        for (int j{}; j < N; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << (xl+i*dx) << "\t" << y[jj*(N-1)+j] << "\t" << U[jj].phi[i][j] << "\t" << ExactFunction((xl+i*dx),y[jj*(N-1)+j]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution-%04d.vtk",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "# vtk DataFile Version 3.1" << std::endl;
    FileWrite << "Helmholtz solver" << std::endl;
    FileWrite << "ASCII" << std::endl;
    FileWrite << "DATASET STRUCTURED_GRID" << std::endl;
    FileWrite << "DIMENSIONS " << Nx << " " << Ns*N << " " << 1 << std::endl; 
    FileWrite << "POINTS " << Nx*Ns*N << " FLOAT" << std::endl;
    
    for (int jj{}; jj < Ns; ++jj)
        for (int j{}; j < N; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << (xl+i*dx) << "\t" << y[jj*(N-1)+j] << "\t" << 0.0 << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nx*Ns*N << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int jj{}; jj < Ns; ++jj)
        for (int j{}; j < N; ++j)
            for (int i{}; i < Nx; ++i)
                FileWrite << U[jj].phi[i][j] << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::out);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        Output << "!NBLOCKS " << size << std::endl;
        
        Output.close();
    }
    
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::app);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        for (int i{}; i < size; ++i)
        {
            sprintf(s,"Solution-%04d.vtk",i);
            
            Output << s << std::endl;
        }
        
        Output.close();
    }
    #endif
    
    Deallocate(s,200);
}
