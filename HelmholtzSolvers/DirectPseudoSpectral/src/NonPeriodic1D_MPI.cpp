//_______________________________________________________________________________
//_______________________________________________________________________________
// MPI-parallel Program for 1D Helmholtz equation.
// 
// Developed by: Dr. Pritam Giri
// Email: pritam.jumech@gmail.com
// Date : 15.03.2023
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

#include "NonPeriodic1D_MPI.h"

struct PDE
{
    double *ul, *ur;
    
    double zetal, zetar, b, c, rhs;
};

PDE *S;

double br, rhsr;

int Rank, size;

int Nx, N, Ns;

double xl, xr, Lx;

double alphal, betal, alphar, betar, Lambda;

double **Dx, **Dxx, **Operator, *x, **rhsPDE, leftBC, rightBC, determinant;

int *npivot;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( const int Nx, 
                     const int N, 
                     const double xl, 
                     const double xr )
{
    ::Nx = Nx;
    ::N  = N;
    ::xl = xl;
    ::xr = xr;
    
    Lx = (xr - xl);
    
    Ns = ceil(((double)(Nx))/(size*N));
    
    S = new (std::nothrow) PDE [Ns];
    
    if (S == nullptr)
        ErrorMessage("Error: Memory can not be allocated.");
    
    for (int ii{}; ii < Ns; ++ii)
    {
        Allocate(S[ii].ul,N);
        Allocate(S[ii].ur,N);
    }
    
    Allocate(x,(N-1)*Ns+1);
    Allocate(rhsPDE,Ns,N);
    
    Allocate(Dx,N,N);
    Allocate(Dxx,N,N);
    
    Allocate(Operator,N-2,N-2);
    Allocate(npivot,N-2);
    
    double Lxp = Lx/size;
    
    double *xlocal;
    
    Allocate(xlocal,N);
    
    GaussLobattoLegendrePoints(xlocal,N);
    
    for (int ii{}; ii < Ns; ++ii)
    {
        double xl, xr;
        
        xl = Rank*Lxp + ii*Lxp/Ns;
        xr = xl + Lxp/Ns;
        
        for (int i{}; i < N; ++i)
            x[ii*(N-1)+i] = 0.5*(xl + xr + (Lxp/Ns)*xlocal[i]);
    }
    
    FirstDerivativeMatrix(Dx,xlocal,N);
    SecondDerivativeMatrix(Dxx,Dx,N);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dx[i][j] *= ((2.0*Ns)/Lxp);
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dxx[i][j] *= ((4.0*Ns*Ns)/(Lxp*Lxp));
    
    Deallocate(xlocal,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    for (int ii{}; ii < Ns; ++ii)
    {
        Deallocate(S[ii].ul,N);
        Deallocate(S[ii].ur,N);
    }
    
    delete [] S;
    
    S = nullptr;
    
    Deallocate(x,(N-1)*Ns+1);
    Deallocate(rhsPDE,Ns,N);
    
    Deallocate(Dx,N,N);
    Deallocate(Dxx,N,N);
    
    Deallocate(Operator,N-2,N-2);
    Deallocate(npivot,N-2);
}

//________________________________________________________________________________________________
// This program returns first derivative on 'k'th Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double FirstDerivative ( int k, 
                         double *h )
{
    return FirstDerivative(k,h,Dx,N);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void ConstructFundamentalSolutions ( FieldVariable U[] )
{
    for (int jj{}; jj < Ns; ++jj)
    {
        U[jj].phi[0]   = 0.0;
        U[jj].phi[N-1] = 0.0;
        
        U[jj].phi++;
        
        for (int i{}; i < (N-2); ++i)
            U[jj].phi[i] = rhsPDE[jj][i+1];
        
        LUSolve(Operator,npivot,U[jj].phi,N-2);
        
        U[jj].phi--;
    }
    
    S[0].ul[0]   = 1.0;
    S[0].ul[N-1] = 0.0;
    S[0].ur[0]   = 0.0;
    S[0].ur[N-1] = 1.0;
    
    S[0].ul++;
    S[0].ur++;
    
    for (int i{}; i < (N-2); ++i)
        S[0].ul[i] = -Dxx[i+1][0];
    
    LUSolve(Operator,npivot,S[0].ul,N-2);
    
    for (int i{}; i < (N-2); ++i)
        S[0].ur[i] = -Dxx[i+1][N-1];
    
    LUSolve(Operator,npivot,S[0].ur,N-2);
    
    S[0].ul--;
    S[0].ur--;
    
    for (int jj = 1; jj < Ns; ++jj)
    {
        for (int i{}; i < N; ++i)
        {
            S[jj].ul[i] = S[0].ul[i];
            S[jj].ur[i] = S[0].ur[i];
        }
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SolvePDE ( FieldVariable U[], 
                void (*ConstructOperator)() )
{
    // Step 1: Construct fundamental solutions on every smallest sub-domains
    ConstructOperator();
    ConstructFundamentalSolutions(U);
    
    // Step 2: TDMA forward sweep
    if (size == 1)
    {
        S[0].b = -alphal*FirstDerivative(0,S[0].ul) + betal;
        S[0].c = -alphal*FirstDerivative(0,S[0].ur);
        S[0].rhs = leftBC + alphal*FirstDerivative(0,U[0].phi);
        
        for (int ii = 1; ii < Ns; ++ii)
        {
            leftBC = FirstDerivative(N-1,S[ii-1].ul);
            
            S[ii].b = FirstDerivative(N-1,S[ii-1].ur) - FirstDerivative(0,S[ii].ul) - (S[ii-1].c/S[ii-1].b)*leftBC;
            S[ii].c = -FirstDerivative(0,S[ii].ur);
            S[ii].rhs = FirstDerivative(0,U[ii].phi) - FirstDerivative(N-1,U[ii-1].phi) - (S[ii-1].rhs/S[ii-1].b)*leftBC;
        }
        
        leftBC = alphar*FirstDerivative(N-1,S[Ns-1].ul);
        
        br = alphar*FirstDerivative(N-1,S[Ns-1].ur) + betar - (S[Ns-1].c/S[Ns-1].b)*leftBC;
        rhsr = rightBC - alphar*FirstDerivative(N-1,U[Ns-1].phi) - (S[Ns-1].rhs/S[Ns-1].b)*leftBC;
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == 0)
            {
                S[0].b = -alphal*FirstDerivative(0,S[0].ul) + betal;
                S[0].c = -alphal*FirstDerivative(0,S[0].ur);
                S[0].rhs = leftBC + alphal*FirstDerivative(0,U[0].phi);
            }
            else
            {
                MPI_Recv(&(S[0].zetal), 1, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar), 1, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                S[0].b =  -FirstDerivative(0,S[0].ul) + S[0].zetal;
                S[0].c =  -FirstDerivative(0,S[0].ur);
                S[0].rhs = FirstDerivative(0,U[0].phi) + S[0].zetar;
            }
            
            if (Rank == (size-1))
            {
                leftBC = alphar*FirstDerivative(N-1,S[0].ul);
                
                br = alphar*FirstDerivative(N-1,S[0].ur) + betar - (S[0].c/S[0].b)*leftBC;
                rhsr = rightBC - alphar*FirstDerivative(N-1,U[0].phi) - (S[0].rhs/S[0].b)*leftBC;
            }
            else
            {
                leftBC = FirstDerivative(N-1,S[0].ul);
                
                S[0].zetal = FirstDerivative(N-1,S[0].ur) - (S[0].c/S[0].b)*leftBC;
                S[0].zetar = -FirstDerivative(N-1,U[0].phi) - (S[0].rhs/S[0].b)*leftBC;
                
                MPI_Send(&(S[0].zetal), 1, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[0].zetar), 1, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
        else
        {
            if (Rank == 0)
            {
                S[0].b = -alphal*FirstDerivative(0,S[0].ul) + betal;
                S[0].c = -alphal*FirstDerivative(0,S[0].ur);
                S[0].rhs = leftBC + alphal*FirstDerivative(0,U[0].phi);
            }
            else
            {
                MPI_Recv(&(S[0].zetal), 1, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&(S[0].zetar), 1, MPI_DOUBLE, (Rank-1), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                S[0].b =  -FirstDerivative(0,S[0].ul) + S[0].zetal;
                S[0].c =  -FirstDerivative(0,S[0].ur);
                S[0].rhs = FirstDerivative(0,U[0].phi) + S[0].zetar;
            }
            
            for (int ii = 1; ii < Ns; ++ii)
            {
                leftBC = FirstDerivative(N-1,S[ii-1].ul);
                
                S[ii].b = FirstDerivative(N-1,S[ii-1].ur) - FirstDerivative(0,S[ii].ul) - (S[ii-1].c/S[ii-1].b)*leftBC;
                S[ii].c = -FirstDerivative(0,S[ii].ur);
                S[ii].rhs = FirstDerivative(0,U[ii].phi) - FirstDerivative(N-1,U[ii-1].phi) - (S[ii-1].rhs/S[ii-1].b)*leftBC;
            }
            
            if (Rank == (size-1))
            {
                leftBC = alphar*FirstDerivative(N-1,S[Ns-1].ul);
                
                br = alphar*FirstDerivative(N-1,S[Ns-1].ur) + betar - (S[Ns-1].c/S[Ns-1].b)*leftBC;
                rhsr = rightBC - alphar*FirstDerivative(N-1,U[Ns-1].phi) - (S[Ns-1].rhs/S[Ns-1].b)*leftBC;
            }
            else
            {
                leftBC = FirstDerivative(N-1,S[Ns-1].ul);
                
                S[Ns-1].zetal = FirstDerivative(N-1,S[Ns-1].ur) - (S[Ns-1].c/S[Ns-1].b)*leftBC;
                S[Ns-1].zetar = -FirstDerivative(N-1,U[Ns-1].phi) - (S[Ns-1].rhs/S[Ns-1].b)*leftBC;
                
                MPI_Send(&(S[Ns-1].zetal), 1, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD);
                MPI_Send(&(S[Ns-1].zetar), 1, MPI_DOUBLE, (Rank+1), 1, MPI_COMM_WORLD);
            }
        }
    }
    
    // Step 3: TDMA backward sweep
    if (size == 1)
    {
        for (int ii = (Ns-1); ii >= 0; --ii)
        {
            S[ii].zetar = (ii == (Ns-1) ? rhsr/br : S[ii+1].zetal);
            S[ii].zetal = (S[ii].rhs - S[ii].c * S[ii].zetar)/S[ii].b;
        }
    }
    else
    {
        if (Ns == 1)
        {
            if (Rank == (size-1))
                S[0].zetar = rhsr/br;
            else
                MPI_Recv(&(S[0].zetar), 1, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            S[0].zetal = (S[0].rhs - S[0].c * S[0].zetar)/S[0].b;
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal), 1, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
        else
        {
            if (Rank == (size-1))
                S[Ns-1].zetar = rhsr/br;
            else
                MPI_Recv(&(S[Ns-1].zetar), 1, MPI_DOUBLE, (Rank+1), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            S[Ns-1].zetal = (S[Ns-1].rhs - S[Ns-1].c * S[Ns-1].zetar)/S[Ns-1].b;
            
            for (int ii = (Ns-2); ii >= 0; --ii)
            {
                S[ii].zetar = S[ii+1].zetal;
                S[ii].zetal = (S[ii].rhs - S[ii].c * S[ii].zetar)/S[ii].b;
            }
            
            if (Rank != 0)
                MPI_Send(&(S[0].zetal), 1, MPI_DOUBLE, (Rank-1), 0, MPI_COMM_WORLD);
        }
    }
    
    // Step 4: Construct final solution
    for (int ii{}; ii < Ns; ++ii)
        for (int i{}; i < N; ++i)
            U[ii].phi[i] += S[ii].zetal*S[ii].ul[i] + S[ii].zetar*S[ii].ur[i];
}

//_______________________________________________________________________________
// First derivative with respect to x
//_______________________________________________________________________________
void FirstDerivativeX ( FieldVariable Ux[], 
                        FieldVariable U[] )
{
    double ax;
    
    for (int ii{}; ii < Ns; ++ii)
    {
        for (int i{}; i < N; ++i)
        {
            ax = 0.0;
            
            for (int j{}; j < N; ++j)
                ax += Dx[i][j]*U[ii].phi[j];
            
            Ux[ii].phi[i] = ax;
        }
    }
}

//_______________________________________________________________________________
// Second derivative with respect to x
//_______________________________________________________________________________
void SecondDerivativeX ( FieldVariable Uxx[], 
                         FieldVariable U[] )
{
    double axx;
    
    for (int ii{}; ii < Ns; ++ii)
    {
        for (int i{}; i < N; ++i)
        {
            axx = 0.0;
            
            for (int j{}; j < N; ++j)
                axx += Dxx[i][j]*U[ii].phi[j];
            
            Uxx[ii].phi[i] = axx;
        }
    }
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( FieldVariable U[], 
                    double (*ExactFunction)(const double) )
{
    double Error, LocalError = 0.0;
    
    for (int ii{}; ii < Ns; ++ii)
        for (int i{}; i < N; ++i)
            LocalError = Maximum(LocalError,Absolute(U[ii].phi[i]-ExactFunction(x[ii*(N-1)+i])));
    
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
                 double (*ExactFunction)(const double) )
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
    
    for (int ii{}; ii < Ns; ++ii)
        for (int i{}; i < N; ++i)
            FileWrite << x[ii*(N-1)+i] << "\t" << U[ii].phi[i] << "\t" << ExactFunction(x[ii*(N-1)+i]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution-%04d.curve",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    if (Rank == 0)
        FileWrite << "#Computed" << std::endl;
    
    for (int ii{}; ii < Ns; ++ii)
        for (int i{}; i < N; ++i)
            FileWrite << x[ii*(N-1)+i] << "\t" << U[ii].phi[i] << std::endl;
    
    FileWrite.close();
    #endif
    
    Deallocate(s,200);
}
