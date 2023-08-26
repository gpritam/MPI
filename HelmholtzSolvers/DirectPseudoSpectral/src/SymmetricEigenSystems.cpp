#include "SymmetricEigenSystems.h"

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void CheckSymmetry ( double **A, 
                     const int N )
{
    for (int i{}; i < (N-1); ++i)
    {
        for (int j{i+1}; j < N; ++j)
        {
            if (Absolute(A[i][j]-A[j][i]) > 1.0E-8)
            {
                std::cout << "Difference = " << Absolute(A[i][j]-A[j][i]) << std::endl;
                
                ErrorMessage("This matrix is not a symmetric matrix!");
            }
        }
    }
}

//_______________________________________________________________________________
// 
//_______________________________________________________________________________
void Rotate ( double **&a, 
              const double s, 
              const double tau, 
              const int i, 
              const int j, 
              const int k, 
              const int l )
{
    double g = a[i][j];
    double h = a[k][l];
    
    a[i][j] = g - s*(h+g*tau);
    a[k][l] = h + s*(g-h*tau);
}

//_______________________________________________________________________________
// Jacobi rotation method for symmetric matrix
// 'A' is NxN symmetric matrix. 'Q' contains all N eigenvectors. At the end of 
// the execution, 'Lambda' will contain N eigenvalues.
// 
// A Q = Q Lambda
// 
// Reference : Numerical recipes : The art of scientific computing, 3rd edition, pp. 576
//_______________________________________________________________________________
void JacobiMethodSymmetric ( double **&A,
                             double **&Q, 
                             double *&Lambda, 
                             const int N )
{
    CheckSymmetry(A,N);
    
    int i, j, ip, iq, nrot{};
    
    double tresh, theta, tau, t, sm, s, h, g, c;
    
    double *b, *z;
    
    Allocate(b,N);
    Allocate(z,N);
    
    for (ip = 0; ip < N; ++ip)
    {
        for (iq = 0; iq < N; ++iq)
            Q[ip][iq] = 0.0;
        
        Q[ip][ip] = 1.0;
    }
    
    for (ip = 0; ip < N; ++ip)
    {
        b[ip] = Lambda[ip] = A[ip][ip];
        z[ip] = 0.0;
    }
    
    for (i = 1; i <= 50; ++i)
    {
        sm = 0.0;
        
        for (ip = 0; ip < N-1; ++ip)
        {
            for (iq = ip+1; iq < N; ++iq)
                sm += Absolute(A[ip][iq]);
        }
        
        if (sm == 0.0)
        {
            int k;
            
            for (int i{}; i < N-1; ++i)
            {
                double p = Lambda[k = i];
                
                for (int j{i}; j < N; ++j)
                    if (Lambda[j] >= p)
                        p = Lambda[k = j];
                
                if (k != i)
                {
                    Lambda[k] = Lambda[i];
                    Lambda[i] = p;
                    
                    for (int j{}; j < N; ++j)
                        Swap(Q[j][i],Q[j][k]);
                }
            }
            
            return;
        }
        
        if (i < 4)
            tresh = 0.2*sm/(N*N);
        else
            tresh = 0.0;
        
        for (ip = 0; ip < N-1; ++ip)
        {
            for (iq = ip+1; iq < N; ++iq)
            {
                g = 100.0*Absolute(A[ip][iq]);
                
                if (i > 4 && g <= EPSILON*Absolute(Lambda[ip]) && g <= EPSILON*Absolute(Lambda[iq]))
                    A[ip][iq]=0.0;
                else if (Absolute(A[ip][iq]) > tresh)
                {
                    h = Lambda[iq] - Lambda[ip];
                    
                    if (g <= EPSILON*Absolute(h))
                        t = A[ip][iq]/h;
                    else
                    {
                        theta = 0.5*h/A[ip][iq];
                        t = 1.0/(Absolute(theta)+sqrt(1.0+theta*theta));
                        
                        if (theta < 0.0)
                            t = -t;
                    }
                    
                    c = 1.0/sqrt(1+t*t);
                    s = t*c;
                    tau = s/(1.0+c);
                    h = t*A[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    Lambda[ip] -= h;
                    Lambda[iq] += h;
                    A[ip][iq] = 0.0;
                    
                    for (j = 0; j < ip; ++j)
                        Rotate(A,s,tau,j,ip,j,iq);
                    
                    for (j = ip+1; j < iq; ++j)
                        Rotate(A,s,tau,ip,j,j,iq);
                    
                    for (j = iq+1; j < N; ++j)
                        Rotate(A,s,tau,ip,j,iq,j);
                    
                    for (j = 0; j < N; ++j)
                        Rotate(Q,s,tau,j,ip,j,iq);
                    
                    nrot++;
                }
            }
        }
        
        for (ip = 0; ip < N; ++ip)
        {
            b[ip] += z[ip];
            Lambda[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    ErrorMessage("Too many iterations in routine jacobi");
    
    Deallocate(b,N);
    Deallocate(z,N);
}

//_______________________________________________________________________________
// Householder transformation to make a full symmetric matrix to symmetric 
// tridiagonal matrix. 'A' is a real, symmetric matrix of dimension NxN. On output, 
// 'A' is replaced by the orthogonal matrix 'Q' effecting the transformation. 
// 'Ad' returns the diagonal elements of the tridiagonal matrix, and 
// 'As' the subdiagonal elements. 
// 'Ad' and 'As' are vectors of length N. As[N-1] is arbitrary.
// 
// T = QAQ, where 'T' is the tridiagonal matrix having diagonal elements Ad and subdiagonal elements As
// 'A' is the input matrix and 'Q' is the output orthogonal matrix that is returned in A.
// 
// Reference : Numerical recipes : The art of scientific computing, 3rd edition, pp. 582
//_______________________________________________________________________________
void HouseholderTransformationSymmetric ( double **&A, 
                                          double *&Ad, 
                                          double *&As, 
                                          const int N )
{
    CheckSymmetry(A,N);
    
    double scale, hh, h, g, f;
    
    for (int i{N-1}; i > 0; i--)
    {
        h = scale = 0.0;
        
        if (i > 1)
        {
            for (int k{}; k < i; ++k)
                scale += Absolute(A[i][k]);
            
            if (scale == 0.0)
                As[i] = A[i][i-1];
            else
            {
                for (int k{}; k < i; ++k)
                {
                    A[i][k] /= scale;
                    h += A[i][k]*A[i][k];
                }
                
                f = A[i][i-1];
                
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                
                As[i] = scale*g;
                
                h -= f*g;
                
                A[i][i-1] = f-g;
                
                f = 0.0;
                
                for (int j{}; j < i; ++j)
                {
                    A[j][i] = A[i][j]/h;
                    
                    g = 0.0;
                    
                    for (int k{}; k < j+1; ++k)
                        g += A[j][k]*A[i][k];
                    
                    for (int k{j+1}; k < i; ++k)
                        g += A[k][j]*A[i][k];
                    
                    As[j] = g/h;
                    f += As[j]*A[i][j];
                }
                
                hh = 0.5*f/h;
                
                for (int j{}; j < i; ++j)
                {
                    f = A[i][j];
                    
                    As[j] = g = As[j]-hh*f;
                    
                    for (int k{}; k < j+1; ++k)
                        A[j][k] -= (f*As[k]+g*A[i][k]);
                }
            }
        }
        else
            As[i] = A[i][i-1];
        
        Ad[i] = h;
    }
    
    Ad[0] = 0.0;
    
    for (int i{}; i < N; ++i)
    {
        if (Ad[i] != 0.0)
        {
            for (int j{}; j < i; ++j)
            {
                g = 0.0;
                
                for (int k{}; k < i; ++k)
                    g += A[i][k]*A[k][j];
                
                for (int k{}; k < i; ++k)
                    A[k][j] -= g*A[k][i];
            }
        }
        
        Ad[i] = A[i][i];
        
        A[i][i] = 1.0;
        
        for (int j{}; j < i; ++j)
            A[j][i]=A[i][j]=0.0;
    }
    
    for (int i{1}; i < N; ++i)
        As[i-1] = As[i];
}

//_______________________________________________________________________________
// QR algorithm with implicit shifts, to determine the eigenvalues and 
// eigenvectors of a real, symmetric, tridiagonal matrix. Ad[N] contains the 
// diagonal elements of the tridiagonal matrix. On output, it returns the 
// eigenvalues. The vector As[N] inputs the subdiagonal elements of the 
// tridiagonal matrix (As[N-1] is arbitrary). The kth column of Q[N][N] 
// returns the normalized eigenvector corresponding to Ad[k].
// 
// AQ = Q Lambda
// 
// Reference : Numerical recipes : The art of scientific computing, 3rd edition, pp. 588
//_______________________________________________________________________________
void FrancisMethodTridiagonalSymmetric ( double *&Ad,
                                         double *As,
                                         double **&Q,
                                         const int N, 
                                         const bool InitiallyTridiagonal )
{
    // Initialize Q = I.
    if (InitiallyTridiagonal)
    {
        for (int i{}; i < N; ++i)
            for (int j{}; j < N; ++j)
                Q[i][j] = (i == j ? 1.0 : 0.0);
    }
    
    int i, l, m, k;
    
    double s, r, p, g, f, d, c, b;
    
    for (l = 0; l < N; ++l)
    {
        int Iteration = 0;
        
        do
        {
            // Look for a single small subdiagonal element to split the matrix.
            for (m = l; m < (N-1); ++m)
            {
                d = Absolute(Ad[m])+Absolute(Ad[m+1]);
                
                if ((double)(Absolute(As[m]) + d) == d)
                    break;
            }
            
            if (m != l)
            {
                if (Iteration++ == 30)
                    std::cout << "Too many iterations!" << std::endl;
                
                // Form shift.
                g = 0.5*(Ad[l+1]-Ad[l])/As[l];
                r = SquareRootSquaredSum(g,1.0);
                g = Ad[m]-Ad[l]+As[l]/(g+Sign(r,g));
                s = c = 1.0;
                p = 0.0;
                
                // A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
                for (i = (m-1); i >= l; i--)
                {
                    f = s*As[i];
                    b = c*As[i];
                    As[i+1] = (r = SquareRootSquaredSum(f,g));
                    
                    if (r == 0.0)
                    {
                        Ad[i+1] -= p;
                        As[m] = 0.0;
                        break;
                    }
                    
                    s = (f / r);
                    c = (g / r);
                    g = (Ad[i+1] - p);
                    r = (Ad[i] - g)*s + 2.0*c*b;
                    Ad[i+1] = g + (p = s*r);
                    g = (c*r - b);
                    
                     // Form eigenvectors.                    
                    for (k = 0; k < N; ++k)
                    {
                        f = Q[k][i+1];
                        Q[k][i+1] = (s*Q[k][i] + c*f);
                        Q[k][i] = (c*Q[k][i] - s*f);
                    }
                }
                
                if ((r == 0.0) && (i >= l))
                    continue;
                
                Ad[l] -= p;
                As[l] = g;
                As[m] = 0.0;
            }
        } while (m != l);
    }
    
    // Sort eigenvalues
    for (i = 0; i < (N-1); ++i)
    {
        p = Ad[k = i];
        
        for (int j{i}; j < N; ++j)
            if (Ad[j] >= p)
                p = Ad[k = j];
        
        if (k != i)
        {
            Ad[k] = Ad[i];
            Ad[i] = p;
            
            for (int j{}; j < N; ++j)
                Swap(Q[j][i],Q[j][k]);
        }
    }
}

//_______________________________________________________________________________
// AQ = Q Lambda
// 
// Eigenvectors are returned in 'A' and eigenvalues are returned in 'Lambda'.
// A      : N x N
// Lambda : N
//_______________________________________________________________________________
void EigenSystemSymmetricMatrix ( double **&A, 
                                  double *&Lambda, 
                                  const int N )
{
    double *As;
    
    Allocate(As,N);
    
    HouseholderTransformationSymmetric(A,Lambda,As,N);
    FrancisMethodTridiagonalSymmetric(Lambda,As,A,N,false);
    
    Deallocate(As,N);
}
