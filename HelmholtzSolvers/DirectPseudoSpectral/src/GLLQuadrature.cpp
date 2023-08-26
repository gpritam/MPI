#include "GLLQuadrature.h"

//________________________________________________________________________________________________
// This program returns the value of Legendre polynomial of order 'n' at a given point 'x', 
// where |x| <= 1 and n >= 0. 
//________________________________________________________________________________________________
double Legendre ( double x, 
                  const int n )
{
    // Check the validity of input variables.
    if (Absolute(x) > (1.0 + EPSILON))
        ErrorMessage("Absolute value of the argument of the Legendre polynomial cannot be more than 1!");
    
    if (x > 1.0)
        x = 1.0;
    if (x < -1.0)
        x = -1.0;
    
    if (n < 0)
        ErrorMessage("Order of Legendre polynomial cannot be less than zero!");
    
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return x;
    else
    {
        double P0 = 1.0, P1 = x, Pn;
        
        for (int i{2}; i <= n; ++i)
        {
            Pn = ((2.0*i-1.0)*x*P1 - (i-1.0)*P0)/i;
            
            P0 = P1;
            P1 = Pn;
        }
        
        return Pn;
    }
}

//________________________________________________________________________________________________
// This program returns the value of Legendre polynomial of any degree <= 'n' at points 'x', 
// where |magnitude of any element of x| <= 1 and n >= 0. 
// 
// x               : N
// LegendreValues : N x N : (number of node points) x (maximum degree - 1)
//________________________________________________________________________________________________
void GetLegendreValues ( double **&LegendreValues, 
                         double *&x, 
                         const int N )
{
    // Check the validity of input variables.
    for (int i{}; i < N; ++i)
        if (Absolute(x[i]) > (1.0 + EPSILON))
            ErrorMessage("Absolute value of the argument of the Legendre polynomial cannot be more than 1!");
    
    // Loop over all nodes
    for (int i{}; i < N; ++i)
    {
        if (x[i] > 1.0)
            x[i] = 1.0;
        
        if (x[i] < -1.0)
            x[i] = -1.0;
        
        for (int n{}; n < N; ++n)
            LegendreValues[i][n] = (n == 0 ? 1.0 : (n == 1 ? x[i] : ((2.0*n-1.0)*x[i]*LegendreValues[i][n-1] - (n-1.0)*LegendreValues[i][n-2])/n));
    }
}

//________________________________________________________________________________________________
// This program returns 'N' Gauss-Lobatto-Legendre points, where N >= 2.
// 
// References : 1. 'High-order Gauss–Lobatto formulae' by Walter Gautschi
//              2. 'Gaussian Quadrature and the Eigenvalue Problem' by John A. Gubner
//________________________________________________________________________________________________
void GaussLobattoLegendrePoints ( double *&x, 
                                  const int N )
{
    // Check the validity of input variables.
    if (N < 2)
        ErrorMessage("Number of points has to be more than or equal to 2.");
    
    if (N == 2)
    {
        x[0] = -1.0;
        x[1] = 1.0;
    }
    else
    {
        double *As, **Q;
        
        Allocate(As,N);
        Allocate(Q,N,N);
        
        for (int i{}; i < (N-1); ++i)
            x[i] = 0.0;
        
        for (int i{}; i < (N-2); ++i)
            As[i] = sqrt((i+1.0)*(i+1.0)/(4.0*(i+1.0)*(i+1.0)-1.0));
        
        double sign = (N % 2 == 0 ? 1.0 : -1.0);
        double a = -sign*(N-1.0)/(2.0*N-3.0);
        double b = sign;
        double c = (N-1.0)/(2.0*N-3.0);
        double d = 1.0;
        
        double e = sign*(N-1.0)/(2.0*N-3.0);
        double f = (N-1.0)/(2.0*N-3.0);
        
        double alpha = (d*e-b*f)/(a*d-b*c);
        double beta  = (a*f-e*c)/(a*d-b*c);
        
        x[N-1] = sqrt(alpha);
        As[N-2] = sqrt(beta);
        
        FrancisMethodTridiagonalSymmetric(x,As,Q,N);
        
        int M = (N % 2 == 0 ? N/2 : (N-1)/2);
        
        for (int i{}; i < M; ++i)
            Swap(x[i],x[N-i-1]);
        
        Deallocate(As,N);
        Deallocate(Q,N,N);
    }
}

//________________________________________________________________________________________________
// This program returns 'N' weights for integration on Gauss-Lobatto-Legendre points, where N >= 2. 
// We can compute weights more accurately this way.
// Reference : 'High-order Gauss–Lobatto formulae' by Walter Gautschi
//________________________________________________________________________________________________
void GaussLobattoLegendreWeights ( double *&w, 
                                   double *x, 
                                   const int N )
{
    int points = (N % 2 == 0 ? N/2 : (N+1)/2);
    
    double numerator = 2.0/((N-1.0)*N), Pn;
    
    for (int i{}; i < points; ++i)
    {
        Pn = Legendre(x[i],N-1);
        
        w[i] = numerator/(Pn*Pn);
        
        w[N-i-1] = w[i];
    }
}

//________________________________________________________________________________________________
// This program returns 'N' spectral coefficients from functional values on Gauss-Lobatto-Legendre 
// points, where N >= 2. 
//________________________________________________________________________________________________
void PhysicalToSpectral ( double *&H, 
                          double *h, 
                          double *x, 
                          double *w, 
                          const int N )
{
    double Gammak;
    
    for (int k{}; k < N; ++k)
    {
        Gammak = (k == (N-1) ? 2.0/(N-1.0) : 2.0/(2.0*k+1.0));
        
        H[k] = 0.0;
        
        for (int i{}; i < N; ++i)
            H[k] += h[i]*Legendre(x[i],k)*w[i];
        
        H[k] /= Gammak;
    }
}

//________________________________________________________________________________________________
// This program returns 'N' functional values from spectral coefficients on Gauss-Lobatto-Legendre 
// points, where N >= 2. 
//________________________________________________________________________________________________
void SpectralToPhysical ( double *&h, 
                          double *H, 
                          double *x, 
                          const int N )
{
    for (int i{}; i < N; ++i)
        h[i] = 0.0;
    
    for (int i{}; i < N; ++i)
        for (int k{}; k < N; ++k)
            h[i] += H[k]*Legendre(x[i],k);
}

//________________________________________________________________________________________________
// This program returns 'NxN' derivative matrix (Dx) on Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
void FirstDerivativeMatrix ( double **&Dx, 
                             double *x, 
                             const int N )
{
    Dx[0][0] = -0.25*N*(N-1);
    
    Dx[N-1][N-1] = 0.25*N*(N-1);
    
    for (int i{1}; i < (N-1); ++i)
        Dx[i][i] = 0.0;
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            if (i != j)
                Dx[i][j] = Legendre(x[i],N-1)/(Legendre(x[j],N-1)*(x[i]-x[j]));
}

//________________________________________________________________________________________________
// This program calculates derivatives on 'N' Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
void FirstDerivative ( double *&hp, 
                       double *h, 
                       double **Dx, 
                       const int N )
{
    for (int i{}; i < N; ++i)
        hp[i] = 0.0;
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            hp[i] += Dx[i][j]*h[j];
}

//________________________________________________________________________________________________
// This program returns first derivative on 'k'th Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double FirstDerivative ( int k, 
                         double *h, 
                         double **Dx, 
                         const int N )
{
    double sum = 0.0;
    
    for (int i{}; i < N; ++i)
        sum += Dx[k][i]*h[i];
    
    return sum;
}

//________________________________________________________________________________________________
// This program returns 'NxN' second derivative matrix (Dxx) on Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
void SecondDerivativeMatrix ( double **&Dxx, 
                              double **Dx, 
                              const int N )
{
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            Dxx[i][j] = 0.0;
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            for (int k{}; k < N; ++k)
                Dxx[i][j] += Dx[i][k]*Dx[k][j];
}

//________________________________________________________________________________________________
// This program calculates second derivatives on 'N' Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
void SecondDerivative ( double *&hpp, 
                        double *h, 
                        double **Dxx, 
                        const int N )
{
    for (int i{}; i < N; ++i)
        hpp[i] = 0.0;
    
    for (int i{}; i < N; ++i)
        for (int j{}; j < N; ++j)
            hpp[i] += Dxx[i][j]*h[j];
}

//________________________________________________________________________________________________
// This program returns second derivative on 'k'th Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double SecondDerivative ( int k, 
                          double *h, 
                          double **Dxx, 
                          const int N )
{
    double sum = 0.0;
    
    for (int i{}; i < N; ++i)
        sum += Dxx[k][i]*h[i];
    
    return sum;
}

//________________________________________________________________________________________________
// This program returns integration on [-1,1] from 'N' functional values and weights on 
// Gauss-Lobatto-Legendre points, where N >= 2. 
//________________________________________________________________________________________________
double Integration ( double *h, 
                     double *w, 
                     const int N )
{
    double sum = 0.0;
    
    for (int i{}; i < N; ++i)
        sum += h[i]*w[i];
    
    return sum;
}

//________________________________________________________________________________________________
// (P_i,P_j)
//________________________________________________________________________________________________
double FirstInnerProduct ( const int i, 
                           const int j )
{
    return (i == j ? 1.0/(i+0.5) : 0.0);
}

//________________________________________________________________________________________________
// (d^2 P_i/d zeta^2,P_j)
//________________________________________________________________________________________________
double SecondInnerProduct ( const int i, 
                            const int j )
{
    if ( (j > (i-2)) || (i < 2) )
        return 0.0;
    else
    {
        double sum = 0.0;
        
        if (i % 2 == 0)
            for (int k{}, kmax{(i-2)/2}; k <= kmax; ++k)
                sum += FirstInnerProduct(2*k,j)*(4.0*k+1.0)*(i-2.0*k)*(i+2.0*k+1.0)*0.5;
        else
            for (int k{}, kmax{(i-3)/2}; k <= kmax; ++k)
                sum += FirstInnerProduct(2*k+1,j)*(4.0*k+3.0)*(i-2.0*k-1.0)*(i+2.0*k+2.0)*0.5;
        
        return sum;
    }
}
