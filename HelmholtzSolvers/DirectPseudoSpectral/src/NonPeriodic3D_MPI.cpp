//_______________________________________________________________________________
//_______________________________________________________________________________
// MPI-parallel Program for 3D Helmholtz equation (Domain decomposition along z direction).
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 20.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

#include "NonPeriodic3D_MPI.h"

struct PDE
{
    double ***ul, ***ur;
    
    double **zetal, **zetar, **b, **c, **rhs;
};

PDE *S;

double **br, **rhsr;

int Rank, size;

int Nx, Ny, Nz, NY, N, Ns;

double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, dx, dy;

double alphal, betal, alphar, betar, Lambda;

double **Dz, **Dzz, **Operator, *z, ****rhsPDE, **leftBC, **rightBC, determinant;

int *npivot;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const int Ny, 
                     const int Nz, 
                     const int N, 
                     const double xl, 
                     const double xr, 
                     const double yl, 
                     const double yr, 
                     const double zl, 
                     const double zr )
{
    ::Nx = Nx;
    ::Ny = Ny;
    ::Nz = Nz;
    ::N  = N;
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    ::zl = zl;
    ::zr = zr;
    
    NY = 2*Ny;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    Lz = (zr-zl);
    dx = Lx/Nx;
    dy = Ly/Ny;
    
    Ns = ceil(((double)(Nz))/(size*N));
    
    S = new (std::nothrow) PDE [Ns];
    
    if (S == nullptr)
        ErrorMessage("Error: Memory can not be allocated!");
    
    for (int kk{}; kk < Ns; ++kk)
    {
        Allocate(S[kk].ul,Nx,NY,N);
        Allocate(S[kk].ur,Nx,NY,N);
        
        Allocate(S[kk].zetal,Nx,NY);
        Allocate(S[kk].zetar,Nx,NY);
        Allocate(S[kk].b,Nx,NY);
        Allocate(S[kk].c,Nx,NY);
        Allocate(S[kk].rhs,Nx,NY);
    }
    
    Allocate(z,Ns*(N-1)+1);
    Allocate(rhsPDE,Ns,2*Nx,NY,N);
    
    Allocate(br,Nx,NY);
    Allocate(rhsr,Nx,NY);
    Allocate(leftBC,2*Nx,NY);
    Allocate(rightBC,2*Nx,NY);
    
    Allocate(Dz,N,N);
    Allocate(Dzz,N,N);
    
    Allocate(Operator,N-2,N-2);
    Allocate(npivot,N-2);
    
    double Lzp = Lz/size;
    
    double *zlocal;
    
    Allocate(zlocal,N);
    
    GaussLobattoLegendrePoints(zlocal,N);
    
    for (int kk{}; kk < Ns; ++kk)
    {
        double zl, zr;
        
        zl = Rank*Lzp + kk*Lzp/Ns;
        zr = zl + Lzp/Ns;
        
        for (int k{}; k < N; ++k)
            z[kk*(N-1)+k] = 0.5*(zl + zr + (Lzp/Ns)*zlocal[k]);
    }
    
    FirstDerivativeMatrix(Dz,zlocal,N);
    SecondDerivativeMatrix(Dzz,Dz,N);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dz[i][j] *= ((2.0*Ns)/Lzp);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dzz[i][j] *= ((4.0*Ns*Ns)/(Lzp*Lzp));
    
    Deallocate(zlocal,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    for (int kk{}; kk < Ns; ++kk)
    {
        Deallocate(S[kk].ul,Nx,NY,N);
        Deallocate(S[kk].ur,Nx,NY,N);
        
        Deallocate(S[kk].zetal,Nx,NY);
        Deallocate(S[kk].zetar,Nx,NY);
        Deallocate(S[kk].b,Nx,NY);
        Deallocate(S[kk].c,Nx,NY);
        Deallocate(S[kk].rhs,Nx,NY);
    }
    
    delete [] S;
    
    S = nullptr;
    
    Deallocate(z,Ns*(N-1)+1);
    Deallocate(rhsPDE,Ns,2*Nx,NY,N);
    
    Deallocate(br,Nx,NY);
    Deallocate(rhsr,Nx,NY);
    Deallocate(leftBC,2*Nx,NY);
    Deallocate(rightBC,2*Nx,NY);
    
    Deallocate(Dz,N,N);
    Deallocate(Dzz,N,N);
    
    Deallocate(Operator,N-2,N-2);
    Deallocate(npivot,N-2);
}

//________________________________________________________________________________________________
// This program returns first derivative on 'k'th Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double FirstDerivative ( int k, 
                         double *h )
{
    return FirstDerivative(k,h,Dz,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void ConstructFundamentalSolutions ( FieldVariable U[], 
                                     const int i, 
                                     const int j )
{
    for (int kk{}; kk < Ns; ++kk)
    {
        U[kk].phi[i][j][0]   = 0.0;
        U[kk].phi[i][j][N-1] = 0.0;
        
        U[kk].phi[i][j]++;
        
        for (int k{}; k < (N-2); ++k)
            U[kk].phi[i][j][k] = rhsPDE[kk][i][j][k+1];
        
        LUSolve(Operator,npivot,U[kk].phi[i][j],N-2);
        
        U[kk].phi[i][j]--;
    }
    
    S[0].ul[i][j][0]   = 1.0;
    S[0].ul[i][j][N-1] = 0.0;
    S[0].ur[i][j][0]   = 0.0;
    S[0].ur[i][j][N-1] = 1.0;
    
    S[0].ul[i][j]++;
    S[0].ur[i][j]++;
    
    for (int k{}; k < (N-2); ++k)
        S[0].ul[i][j][k] = -Dzz[k+1][0];
    
    LUSolve(Operator,npivot,S[0].ul[i][j],N-2);
    
    for (int k{}; k < (N-2); ++k)
        S[0].ur[i][j][k] = -Dzz[k+1][N-1];
    
    LUSolve(Operator,npivot,S[0].ur[i][j],N-2);
    
    S[0].ul[i][j]--;
    S[0].ur[i][j]--;
    
    for (int kk = 1; kk < Ns; ++kk)
    {
        for (int k{}; k < N; ++k)
        {
            S[kk].ul[i][j][k] = S[0].ul[i][j][k];
            S[kk].ur[i][j][k] = S[0].ur[i][j][k];
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( FieldVariable U[], 
                void (*ConstructOperator)(const int, const int) )
{
    // Step 1: Take forward Fourier transform of the right hand side and BCs
    for (int kk{}; kk < Ns; ++kk)
        FourierTransformXY(rhsPDE[kk],1,Nx,Ny,N);
    
    FourierTransformXY(leftBC,1,Nx,Ny);
    FourierTransformXY(rightBC,1,Nx,Ny);
    
    // Step 2: Construct fundamental solutions on every smallest sub-domains
    for (int iw{}; iw <= Nx/2; iw++)
    {
        for (int jw{}; jw <= Ny/2; jw++)
        {
            // Construct 'Operator' matrix for wave numbers |iw| and |jw| and compute its LU decomposition.
            ConstructOperator(iw,jw);
            
            // Solve ODEs with this operator
            if (iw == 0)
            {
                if (jw == 0)
                {
                    ConstructFundamentalSolutions(U,0,0);
                    ConstructFundamentalSolutions(U,0,1);
                }
                else if (jw == Ny/2)
                {
                    ConstructFundamentalSolutions(U,0,Ny);
                    ConstructFundamentalSolutions(U,0,Ny+1);
                }
                else
                {
                    ConstructFundamentalSolutions(U,0,2*jw);
                    ConstructFundamentalSolutions(U,0,2*jw+1);
                    ConstructFundamentalSolutions(U,0,2*(Ny-jw));
                    ConstructFundamentalSolutions(U,0,2*(Ny-jw)+1);
                }
            }
            else if (iw == Nx/2)
            {
                if (jw == 0)
                {
                    ConstructFundamentalSolutions(U,Nx/2,0);
                    ConstructFundamentalSolutions(U,Nx/2,1);
                }
                else if (jw == Ny/2)
                {
                    ConstructFundamentalSolutions(U,Nx/2,Ny);
                    ConstructFundamentalSolutions(U,Nx/2,Ny+1);
                }
                else
                {
                    ConstructFundamentalSolutions(U,Nx/2,2*jw);
                    ConstructFundamentalSolutions(U,Nx/2,2*jw+1);
                    ConstructFundamentalSolutions(U,Nx/2,2*(Ny-jw));
                    ConstructFundamentalSolutions(U,Nx/2,2*(Ny-jw)+1);
                }
            }
            else
            {
                if (jw == 0)
                {
                    ConstructFundamentalSolutions(U,iw,0);
                    ConstructFundamentalSolutions(U,iw,1);
                    ConstructFundamentalSolutions(U,(Nx-iw),0);
                    ConstructFundamentalSolutions(U,(Nx-iw),1);
                }
                else if (jw == Ny/2)
                {
                    ConstructFundamentalSolutions(U,iw,Ny);
                    ConstructFundamentalSolutions(U,iw,Ny+1);
                    ConstructFundamentalSolutions(U,(Nx-iw),Ny);
                    ConstructFundamentalSolutions(U,(Nx-iw),Ny+1);
                }
                else
                {
                    ConstructFundamentalSolutions(U,iw,2*jw);
                    ConstructFundamentalSolutions(U,iw,2*jw+1);
                    ConstructFundamentalSolutions(U,iw,2*(Ny-jw));
                    ConstructFundamentalSolutions(U,iw,2*(Ny-jw)+1);
                    ConstructFundamentalSolutions(U,(Nx-iw),2*jw);
                    ConstructFundamentalSolutions(U,(Nx-iw),2*jw+1);
                    ConstructFundamentalSolutions(U,(Nx-iw),2*(Ny-jw));
                    ConstructFundamentalSolutions(U,(Nx-iw),2*(Ny-jw)+1);
                }
            }
        }
    }
    
    // Step 3: TDMA forward sweep
    if (size == 1)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < NY; ++j)
            {
                S[0].b[i][j] = -alphal*FirstDerivative(0,S[0].ul[i][j]) + betal;
                S[0].c[i][j] = -alphal*FirstDerivative(0,S[0].ur[i][j]);
                S[0].rhs[i][j] = leftBC[i][j] + alphal*FirstDerivative(0,U[0].phi[i][j]);
            }
        }
        
        for (int kk = 1; kk < Ns; ++kk)
        {
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < NY; ++j)
                    leftBC[i][j] = FirstDerivative(N-1,S[kk-1].ul[i][j]);
            
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < NY; ++j)
                {
                    S[kk].b[i][j] = FirstDerivative(N-1,S[kk-1].ur[i][j]) - FirstDerivative(0,S[kk].ul[i][j]) - (S[kk-1].c[i][j]/S[kk-1].b[i][j])*leftBC[i][j];
                    S[kk].c[i][j] = -FirstDerivative(0,S[kk].ur[i][j]);
                    S[kk].rhs[i][j] = FirstDerivative(0,U[kk].phi[i][j]) - FirstDerivative(N-1,U[kk-1].phi[i][j]) - (S[kk-1].rhs[i][j]/S[kk-1].b[i][j])*leftBC[i][j];
                }
            }
        }
        
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < NY; ++j)
                leftBC[i][j] = alphar*FirstDerivative(N-1,S[Ns-1].ul[i][j]);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < NY; ++j)
            {
                br[i][j] = alphar*FirstDerivative(N-1,S[Ns-1].ur[i][j]) + betar - (S[Ns-1].c[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
                rhsr[i][j] = rightBC[i][j] - alphar*FirstDerivative(N-1,U[Ns-1].phi[i][j]) - (S[Ns-1].rhs[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
            }
        }
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == 0)
            {
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[0].b[i][j] = -alphal*FirstDerivative(0,S[0].ul[i][j]) + betal;
                        S[0].c[i][j] = -alphal*FirstDerivative(0,S[0].ur[i][j]);
                        S[0].rhs[i][j] = leftBC[i][j] + alphal*FirstDerivative(0,U[0].phi[i][j]);
                    }
                }
            }
            else
            {
                MPI_Recv(&(S[0].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[0].b[i][j] =  -FirstDerivative(0,S[0].ul[i][j]) + S[0].zetal[i][j];
                        S[0].c[i][j] =  -FirstDerivative(0,S[0].ur[i][j]);
                        S[0].rhs[i][j] = FirstDerivative(0,U[0].phi[i][j]) + S[0].zetar[i][j];
                    }
                }
            }
            
            if (Rank == (size-1))
            {
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        leftBC[i][j] = alphar*FirstDerivative(N-1,S[0].ul[i][j]);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        br[i][j] = alphar*FirstDerivative(N-1,S[0].ur[i][j]) + betar - (S[0].c[i][j]/S[0].b[i][j])*leftBC[i][j];
                        rhsr[i][j] = rightBC[i][j] - alphar*FirstDerivative(N-1,U[0].phi[i][j]) - (S[0].rhs[i][j]/S[0].b[i][j])*leftBC[i][j];
                    }
                }
            }
            else
            {
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        leftBC[i][j] = FirstDerivative(N-1,S[0].ul[i][j]);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[0].zetal[i][j] = FirstDerivative(N-1,S[0].ur[i][j]) - (S[0].c[i][j]/S[0].b[i][j])*leftBC[i][j];
                        S[0].zetar[i][j] = -FirstDerivative(N-1,U[0].phi[i][j]) - (S[0].rhs[i][j]/S[0].b[i][j])*leftBC[i][j];
                    }
                }
                
                MPI_Send(&(S[0].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[0].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
        else
        {
            if (Rank == 0)
            {
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[0].b[i][j] = -alphal*FirstDerivative(0,S[0].ul[i][j]) + betal;
                        S[0].c[i][j] = -alphal*FirstDerivative(0,S[0].ur[i][j]);
                        S[0].rhs[i][j] = leftBC[i][j] + alphal*FirstDerivative(0,U[0].phi[i][j]);
                    }
                }
            }
            else
            {
                MPI_Recv(&(S[0].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[0].b[i][j] =  -FirstDerivative(0,S[0].ul[i][j]) + S[0].zetal[i][j];
                        S[0].c[i][j] =  -FirstDerivative(0,S[0].ur[i][j]);
                        S[0].rhs[i][j] = FirstDerivative(0,U[0].phi[i][j]) + S[0].zetar[i][j];
                    }
                }
            }
            
            for (int kk = 1; kk < Ns; ++kk)
            {
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        leftBC[i][j] = FirstDerivative(N-1,S[kk-1].ul[i][j]);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[kk].b[i][j] = FirstDerivative(N-1,S[kk-1].ur[i][j]) - FirstDerivative(0,S[kk].ul[i][j]) - (S[kk-1].c[i][j]/S[kk-1].b[i][j])*leftBC[i][j];
                        S[kk].c[i][j] = -FirstDerivative(0,S[kk].ur[i][j]);
                        S[kk].rhs[i][j] = FirstDerivative(0,U[kk].phi[i][j]) - FirstDerivative(N-1,U[kk-1].phi[i][j]) - (S[kk-1].rhs[i][j]/S[kk-1].b[i][j])*leftBC[i][j];
                    }
                }
            }
            
            if (Rank == (size-1))
            {
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        leftBC[i][j] = alphar*FirstDerivative(N-1,S[Ns-1].ul[i][j]);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        br[i][j] = alphar*FirstDerivative(N-1,S[Ns-1].ur[i][j]) + betar - (S[Ns-1].c[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
                        rhsr[i][j] = rightBC[i][j] - alphar*FirstDerivative(N-1,U[Ns-1].phi[i][j]) - (S[Ns-1].rhs[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
                    }
                }
            }
            else
            {
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        leftBC[i][j] = FirstDerivative(N-1,S[Ns-1].ul[i][j]);
                
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[Ns-1].zetal[i][j] = FirstDerivative(N-1,S[Ns-1].ur[i][j]) - (S[Ns-1].c[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
                        S[Ns-1].zetar[i][j] = -FirstDerivative(N-1,U[Ns-1].phi[i][j]) - (S[Ns-1].rhs[i][j]/S[Ns-1].b[i][j])*leftBC[i][j];
                    }
                }
                
                MPI_Send(&(S[Ns-1].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[Ns-1].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
    }
    
    // Step 4: TDMA backward sweep
    if (size == 1)
    {
        for (int kk = (Ns-1); kk >= 0; --kk)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < NY; ++j)
                {
                    S[kk].zetar[i][j] = (kk == (Ns-1) ? rhsr[i][j]/br[i][j] : S[kk+1].zetal[i][j]);
                    S[kk].zetal[i][j] = (S[kk].rhs[i][j] - S[kk].c[i][j] * S[kk].zetar[i][j])/S[kk].b[i][j];
                }
            }
        }
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == (size-1))
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        S[0].zetar[i][j] = rhsr[i][j]/br[i][j];
            else
                MPI_Recv(&(S[0].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < NY; ++j)
                    S[0].zetal[i][j] = (S[0].rhs[i][j] - S[0].c[i][j] * S[0].zetar[i][j])/S[0].b[i][j];
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
        else
        {
            if (Rank == (size-1))
                for (int i{}; i < Nx; ++i)
                    for (int j{}; j < NY; ++j)
                        S[Ns-1].zetar[i][j] = rhsr[i][j]/br[i][j];
            else
                MPI_Recv(&(S[Ns-1].zetar[0][0]), Nx*NY, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < NY; ++j)
                    S[Ns-1].zetal[i][j] = (S[Ns-1].rhs[i][j] - S[Ns-1].c[i][j] * S[Ns-1].zetar[i][j])/S[Ns-1].b[i][j];
            
            for (int kk = (Ns-2); kk >= 0; --kk)
            {
                for (int i{}; i < Nx; ++i)
                {
                    for (int j{}; j < NY; ++j)
                    {
                        S[kk].zetar[i][j] = S[kk+1].zetal[i][j];
                        S[kk].zetal[i][j] = (S[kk].rhs[i][j] - S[kk].c[i][j] * S[kk].zetar[i][j])/S[kk].b[i][j];
                    }
                }
            }
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal[0][0]), Nx*NY, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
    }
    
    // Step 5: Construct final solution
    for (int kk{}; kk < Ns; ++kk)
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < NY; ++j)
                for (int k{}; k < N; ++k)
                    U[kk].phi[i][j][k] += S[kk].zetal[i][j]*S[kk].ul[i][j][k] + S[kk].zetar[i][j]*S[kk].ur[i][j][k];
    
    for (int kk{}; kk < Ns; ++kk)
        FourierTransformXY(U[kk].phi,-1,Nx,Ny,N);
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
    
}

//_______________________________________________________________________________
// Second derivative with respect to y
//_______________________________________________________________________________
void SecondDerivativeY ( FieldVariable Uyy[], 
                         FieldVariable U[] )
{
    
}

//_______________________________________________________________________________
// First derivative with respect to z
//_______________________________________________________________________________
void FirstDerivativeZ ( FieldVariable Uz[], 
                        FieldVariable U[] )
{
    double az;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int kk{}; kk < Ns; ++kk)
            {
                for (int k{}; k < N; ++k)
                {
                    az = 0.0;
                    
                    for (int p{}; p < N; ++p)
                        az += Dz[k][p]*U[kk].phi[i][j][p];
                    
                    Uz[kk].phi[i][j][k] = az;
                }
            }
        }
    }
}

//_______________________________________________________________________________
// Second derivative with respect to z
//_______________________________________________________________________________
void SecondDerivativeZ ( FieldVariable Uzz[], 
                         FieldVariable U[] )
{
    double azz;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int kk{}; kk < Ns; ++kk)
            {
                for (int k{}; k < N; ++k)
                {
                    azz = 0.0;
                    
                    for (int p{}; p < N; ++p)
                        azz += Dzz[k][p]*U[kk].phi[i][j][p];
                    
                    Uzz[kk].phi[i][j][k] = azz;
                }
            }
        }
    }
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( FieldVariable U[], 
                    double (*ExactFunction)(const double, const double, const double) )
{
    double Error, LocalError = 0.0;
    
    for (int kk{}; kk < Ns; ++kk)
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                for (int k{}; k < N; ++k)
                    LocalError = Maximum(LocalError,Absolute(U[kk].phi[i][j][k]-ExactFunction((xl+i*dx),(yl+j*dy),z[kk*(N-1)+k])));
    
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
                 double (*ExactFunction)(const double, const double, const double) )
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
    FileWrite << "Variables = \"X\",\"Y\",\"Z\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nx << ", J = " << Ny << ", K = " << Ns*N << ", DATAPACKING=POINT" << std::endl;
    
    for (int kk{}; kk < Ns; ++kk)
        for (int k{}; k < N; ++k)
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    FileWrite << (xl+i*dx) << "\t" << (yl+j*dy) << "\t" << z[kk*(N-1)+k] << "\t" << U[kk].phi[i][j][k] << "\t" << ExactFunction((xl+i*dx),(yl+j*dy),z[kk*(N-1)+k]) << std::endl;
    
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
    FileWrite << "DIMENSIONS " << Nx << " " << Ny << " " << Ns*N << std::endl; 
    FileWrite << "POINTS " << Nx*Ny*Ns*N << " FLOAT" << std::endl;
    
    for (int kk{}; kk < Ns; ++kk)
        for (int k{}; k < N; ++k)
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    FileWrite << (xl+i*dx) << "\t" << (yl+j*dy) << "\t" << z[kk*(N-1)+k] << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nx*Ny*Ns*N << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int kk{}; kk < Ns; ++kk)
        for (int k{}; k < N; ++k)
            for (int j{}; j < Ny; ++j)
                for (int i{}; i < Nx; ++i)
                    FileWrite << U[kk].phi[i][j][k] << std::endl;
    
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
