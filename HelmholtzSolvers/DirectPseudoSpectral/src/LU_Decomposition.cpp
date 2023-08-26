#include "LU_Decomposition.h"

//_______________________________________________________________________________
// Solution of a linear system Ax = b using LU decomposition.
// LU decomposition with pivoting.
//_______________________________________________________________________________
void LUDecomposition ( double**& A,
                       int*& npivot,
                       double& determinant,
                       const int N )
{
    int index;
    double *s, *c;
    
    Allocate(s,N);
    Allocate(c,N);
    
    determinant = 1.0;
    
    for (int p{}; p < N; ++p)
    {
        s[p] = Absolute(A[p][0]);
        
        for (int q{1}; q < N; ++q)
        {
            if (Absolute(A[p][q]) > s[p]) 
                s[p] = Absolute(A[p][q]);
        }
    }
    
    for (int k{}; k < (N-1); ++k)
    {
        c[k] = Absolute(A[k][k])/s[k];
        
        index = k;
        
        for (int i{k+1}; i < N; ++i)
        {
            if ((Absolute(A[i][k])/s[i]) > c[k])
            {
                index = i;
                
                c[k] = (Absolute(A[i][k])/s[i]);
            }
        }
        
        npivot[k] = index;
        
        if (Absolute(c[k]) <= EPSILON)
            ErrorMessage("ERROR: Matrix is singular!");
        
        if (index != k)
        {
            determinant = -determinant;
            
            for (int j{k}; j < N; ++j)
                Swap(A[k][j],A[index][j]);
        }
        
        for (int l{k+1}; l < N; ++l)
        {
            A[l][k] = A[l][k]/A[k][k];
            
            for (int n{k+1}; n < N; ++n)
                A[l][n] -= A[l][k]*A[k][n];
        }
        
        determinant = determinant*A[k][k];
    }
    
    determinant = determinant*A[N-1][N-1];
}

//_______________________________________________________________________________
// LU Solver.
// Right hand side 'b' is overwritten by the solution 'x'
//_______________________________________________________________________________
void LUSolve ( double** A,
               int* npivot,
               double*& b,
               const int N )
{
    double temp;
    int index;
    
    for (int k{}; k < (N-1); ++k)
    {
        index = npivot[k];
        
        if (index != k)
            Swap(b[index],b[k]);
        
        for (int i{k+1}; i < N; ++i) 
            b[i] -= A[i][k]*b[k];
    }
    
    b[N-1] /= A[N-1][N-1];
    
    for (int p{N-2}; p >= 0; p--)
    {
        temp = 0;
        
        for (int q{p+1}; q < N; ++q)
            temp += A[p][q]*b[q];
        
        b[p] = (b[p]-temp)/A[p][p];
    }
}
