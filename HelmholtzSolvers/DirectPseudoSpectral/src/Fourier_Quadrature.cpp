#include "Fourier_Quadrature.h"

#ifndef UseFFTWRoutines
//________________________________________________________________________________________________
// Fourier quadratures using in-house routines
//________________________________________________________________________________________________
bool FFTPlanned = false;

int FFT_N = -1;

double *FFT_Temporary, *FFT_Temporary0;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreateFFTPlan ( const int N )
{
    if (!FFTPlanned)
    {
        FFT_N = N;
        FFTPlanned = true;
        
        // Allocate 'FFT_Temporary'
        Allocate(FFT_Temporary,2*N);
        
        // Allocate 'FFT_Temporary0'
        Allocate(FFT_Temporary0,2*N);
    }
    else
    {
        if (N > FFT_N)
        {
            // Deallocate memory first
            Deallocate(FFT_Temporary,2*FFT_N);
            Deallocate(FFT_Temporary0,2*FFT_N);
            
            // Reallocate 'FFT_Temporary'
            Allocate(FFT_Temporary,2*N);
            
            // Reallocate 'FFT_Temporary0'
            Allocate(FFT_Temporary0,2*N);
            
            FFT_N = N;
        }
        
        //ErrorMessage("One FFT plan already exists!");
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyFFTPlan ( const int N )
{
    if (FFTPlanned)
    {
        FFTPlanned = false;
        
        // Deallocate 'FFT_Temporary'
        Deallocate(FFT_Temporary,2*N);
        
        // Deallocate 'FFT_Temporary0'
        Deallocate(FFT_Temporary0,2*N);
        
        FFT_N = -1;
    }
    else
        ErrorMessage("At first, create a FFT plan!");
}

//____________________________________________________________________________________
// This function returns the the Discrete Fourier Transform in 1D.
// 'isign' can be 1 or -1.
//  A : 2*Nx 
//____________________________________________________________________________________
void DFT (  double *&A,
            const int isign,
            const int Nx )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
        RealComplexX(A,Nx,1);
    
    for (int k{}; k < 2*Nx; ++k)
        FFT_Temporary[k] = 0.0;
    
    for (int k{}; k < Nx; ++k)
    {
        for (int j{}; j < Nx; ++j)
        {
            double C = cos(FFT_2PI*j*k/Nx);
            double S = sin(FFT_2PI*j*k/Nx);
            
            FFT_Temporary[2*k]   += A[2*j]  *C + isign*A[2*j+1]*S;
            FFT_Temporary[2*k+1] += A[2*j+1]*C - isign*A[2*j]  *S;
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFT_Temporary[i]/(double)Nx;
        
        RealComplexX(A,Nx,0);
    }
    else
    {
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFT_Temporary[i];
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
//  A : 2*Nx
//____________________________________________________________________________________
void RealComplexX ( double *&A,
                    const int Nx,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int i{}; i < Nx; ++i)
        {
            FFT_Temporary[i] = A[2*i];
            FFT_Temporary[Nx+i] = A[2*i+1];
        }
        
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFT_Temporary[i];
    }
    else
    {
        for (int i{}; i < Nx; ++i)
        {
            FFT_Temporary[2*i] = A[i];
            FFT_Temporary[2*i+1] = A[Nx+i];
        }
        
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFT_Temporary[i];
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
// A : 2*Nx x Ny
//____________________________________________________________________________________
void RealComplexX ( double **&A,
                    const int Nx,
                    const int Ny,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                FFT_Temporary[i] = A[2*i][j];
                FFT_Temporary[Nx+i] = A[2*i+1][j];
            }
            
            for (int i{}; i < 2*Nx; ++i)
                A[i][j] = FFT_Temporary[i];
        }
    }
    else
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                FFT_Temporary[2*i] = A[i][j];
                FFT_Temporary[2*i+1] = A[Nx+i][j];
            }
            
            for (int i{}; i < 2*Nx; ++i)
                A[i][j] = FFT_Temporary[i];
        }
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
// A : Nx x 2*Ny
//____________________________________________________________________________________
void RealComplexY ( double **&A,
                    const int Nx,
                    const int Ny,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Ny > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                FFT_Temporary[j] = A[i][2*j];
                FFT_Temporary[Ny+j] = A[i][2*j+1];
            }
            
            for (int j{}; j < 2*Ny; ++j)
                A[i][j] = FFT_Temporary[j];
        }
    }
    else
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                FFT_Temporary[2*j] = A[i][j];
                FFT_Temporary[2*j+1] = A[i][Ny+j];
            }
            
            for (int j{}; j < 2*Ny; ++j)
                A[i][j] = FFT_Temporary[j];
        }
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
// A : 2*Nx x Ny x Nz
//____________________________________________________________________________________

void RealComplexX ( double ***&A,
                    const int Nx,
                    const int Ny,
                    const int Nz,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int i{}; i < Nx; ++i)
                {
                    FFT_Temporary[i] = A[2*i][j][k];
                    FFT_Temporary[Nx+i] = A[2*i+1][j][k];
                }
                
                for (int i{}; i < 2*Nx; ++i)
                    A[i][j][k] = FFT_Temporary[i];
            }
        }
    }
    else
    {
        for (int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int i{}; i < Nx; ++i)
                {
                    FFT_Temporary[2*i] = A[i][j][k];
                    FFT_Temporary[2*i+1] = A[Nx+i][j][k];
                }
                
                for (int i{}; i < 2*Nx; ++i)
                    A[i][j][k] = FFT_Temporary[i];
            }
        }
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
// A : Nx x 2*Ny x Nz
//____________________________________________________________________________________
void RealComplexY ( double ***&A,
                    const int Nx,
                    const int Ny,
                    const int Nz,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Ny > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < Ny; ++j)
                {
                    FFT_Temporary[j] = A[i][2*j][k];
                    FFT_Temporary[Ny+j] = A[i][2*j+1][k];
                }
                
                for (int j{}; j < 2*Ny; ++j)
                    A[i][j][k] = FFT_Temporary[j];
            }
        }
    }
    else
    {
        for (int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < Ny; ++j)
                {
                    FFT_Temporary[2*j] = A[i][j][k];
                    FFT_Temporary[2*j+1] = A[i][Ny+j][k];
                }
                
                for (int j{}; j < 2*Ny; ++j)
                    A[i][j][k] = FFT_Temporary[j];
            }
        }
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
// A : Nx x Ny x 2*Nz
//____________________________________________________________________________________
void RealComplexZ ( double ***&A,
                    const int Nx,
                    const int Ny,
                    const int Nz,
                    const int arrangement )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nz > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    FFT_Temporary[k] = A[i][j][2*k];
                    FFT_Temporary[Nz+k] = A[i][j][2*k+1];
                }
                
                for (int k{}; k < 2*Nz; ++k)
                    A[i][j][k] = FFT_Temporary[k];
            }
        }
    }
    else
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    FFT_Temporary[2*k] = A[i][j][k];
                    FFT_Temporary[2*k+1] = A[i][j][Nz+k];
                }
                
                for (int k{}; k < 2*Nz; ++k)
                    A[i][j][k] = FFT_Temporary[k];
            }
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1   : Fourier transform
// 'isign' = -1   : Inverse Fourier transform
// 'real' = true  : Sequence of real numbers
// 'real' = false : Sequence of complex numbers
// A              : 2*Nx
//____________________________________________________________________________________
void FourierTransformX ( double *&A,
                         const int isign,
                         const int Nx,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                A[Nx+i] = 0.0;
        
        RealComplexX(A,Nx,1);
    }
    
    int mmax = 2, m, j = 1, istep, nn = Nx;
    double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
    
    // reverse-binary reindexing
    for (int i = 1; i < 2*Nx; i += 2)
    {
        if (j > i)
        {
            Swap(A[j-1],A[i-1]);
            Swap(A[j],A[i]);
        }
            m = nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m /= 2;
        }
        
        j += m;
    }
    
    // Danielson-Lanczos section
    while ( 2*Nx > mmax )
    {
        istep = mmax * 2;
        theta = isign * (-FFT_2PI/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        
        for (int m = 1; m < mmax; m += 2)
        {
            for (int i = m; i <= 2*Nx; i += istep)
            {
                j = i + mmax;
                tempr = wr*A[j-1] - wi*A[j];
                tempi = wr * A[j] + wi*A[j-1];
 
                A[j-1] = A[i-1] - tempr;
                A[j] = A[i] - tempi;
                A[i-1] += tempr;
                A[i] += tempi;
            }
            
            wtemp=wr;
            wr += (wr*wpr - wi*wpi);
            wi += (wi*wpr + wtemp*wpi);
        }
        
        mmax = istep;
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Nx; ++i)
            A[i] /= Nx;
        
        RealComplexX(A,Nx,0);
    }
}

//____________________________________________________________________________________
// 'isign' =  1   : Fourier transform
// 'isign' = -1   : Inverse Fourier transform
// 'real' = true  : Sequence of real numbers
// 'real' = false : Sequence of complex numbers
// A              : 2*Nx x Ny
//____________________________________________________________________________________
void FourierTransformX ( double **&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    A[Nx+i][j] = 0.0;
        
        RealComplexX(A,Nx,Ny,1);
    }
    
    for (int j{}; j < Ny; ++j)
    {
        int mmax = 2, m, jj = 1, istep, nn = Nx;
        double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
         
        // reverse-binary reindexing
        for (int ii = 1; ii < 2*Nx; ii += 2)
        {
            if (jj > ii)
            {
                Swap(A[jj-1][j],A[ii-1][j]);
                Swap(A[jj][j],A[ii][j]);
            }
            
            m = nn;
            while (m >= 2 && jj > m)
            {
                jj -= m;
                m /= 2;
            }
            
            jj += m;
        }
        
        // Danielson-Lanczos section
        while (2*Nx > mmax)
        {
            istep = mmax * 2;
            theta = isign * (-FFT_2PI/mmax);
            wtemp = sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            
            for (int m = 1; m < mmax; m += 2)
            {
                for (int ii = m; ii <= 2*Nx; ii += istep)
                {
                    jj = ii + mmax;
                    tempr = wr*A[jj-1][j] - wi*A[jj][j];
                    tempi = wr*A[jj][j] + wi*A[jj-1][j];
                     
                    A[jj-1][j] = A[ii-1][j]-tempr;
                    A[jj][j] = A[ii][j]-tempi;
                    A[ii-1][j] += tempr;
                    A[ii][j] += tempi;
                }
                
                wtemp=wr;
                wr += (wr*wpr - wi*wpi);
                wi += (wi*wpr + wtemp*wpi);
            }
            
            mmax = istep;
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Nx; ++i)
            for (int j{}; j < Ny; ++j)
                A[i][j] /= Nx;
        
        RealComplexX(A,Nx,Ny,0);
    }
}

//____________________________________________________________________________________
// 'isign' =  1   : Fourier transform
// 'isign' = -1   : Inverse Fourier transform
// 'real' = true  : Sequence of real numbers
// 'real' = false : Sequence of complex numbers
// A              : Nx x 2*Ny
//____________________________________________________________________________________
void FourierTransformY ( double **&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Ny > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    A[i][Ny+j] = 0.0;
        
        RealComplexY(A,Nx,Ny,1);
    }
    
    for (int i{}; i < Nx; ++i)
    {
        int mmax = 2, m, jj = 1, istep, nn = Ny;
        double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
        
        // reverse-binary reindexing
        for (int ii = 1; ii < 2*Ny; ii += 2)
        {
            if (jj > ii)
            {
                Swap(A[i][jj-1],A[i][ii-1]);
                Swap(A[i][jj],A[i][ii]);
            }
            
            m = nn;
            while (m >= 2 && jj > m)
            {
                jj -= m;
                m /= 2;
            }
            
            jj += m;
        }
        // Danielson-Lanczos section
        while (2*Ny > mmax)
        {
            istep = mmax * 2;
            theta = isign * (-FFT_2PI/mmax);
            wtemp = sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            
            for (int m = 1; m < mmax; m += 2)
            {
                for (int ii = m; ii <= 2*Ny; ii += istep)
                {
                    jj = ii + mmax;
                    tempr = wr*A[i][jj-1] - wi*A[i][jj];
                    tempi = wr*A[i][jj] + wi*A[i][jj-1];
                    
                    A[i][jj-1] = A[i][ii-1] - tempr;
                    A[i][jj] = A[i][ii] - tempi;
                    A[i][ii-1] += tempr;
                    A[i][ii] += tempi;
                }
                
                wtemp=wr;
                wr += (wr*wpr - wi*wpi);
                wi += (wi*wpr + wtemp*wpi);
            }
            
            mmax = istep;
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Ny; ++i)
            for (int j{}; j < Nx; ++j)
                A[j][i] /= Ny;
        
        RealComplexY(A,Nx,Ny,0);
    }
}

//____________________________________________________________________________________
// 'isign' =  1   : Fourier transform
// 'isign' = -1   : Inverse Fourier transform
// 'real' = true  : Sequence of real numbers
// 'real' = false : Sequence of complex numbers
// A              : 2*Nx x Ny x Nz
//____________________________________________________________________________________
void FourierTransformX ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Nx > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[Nx+i][j][k] = 0.0;
        
        RealComplexX(A,Nx,Ny,Nz,1);
    }
    
    for (int j{}; j < Ny; ++j)
    {
        for (int k{}; k < Nz; ++k)
        {
            int mmax = 2, m, jj = 1, istep, nn = Nx;
            double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
            
            // reverse-binary reindexing
            for (int ii = 1; ii < 2*Nx; ii += 2)
            {
                if (jj > ii)
                {
                    Swap(A[jj-1][j][k],A[ii-1][j][k]);
                    Swap(A[jj][j][k],A[ii][j][k]);
                }
                
                m = nn;
                while (m >= 2 && jj > m)
                {
                    jj -= m;
                    m /= 2;
                }
                
                jj += m;
            }
            
            // Danielson-Lanczos section
            while (2*Nx > mmax)
            {
                istep = mmax * 2;
                theta = isign * (-FFT_2PI/mmax);
                wtemp = sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi = sin(theta);
                wr = 1.0;
                wi = 0.0;
                
                for (int m = 1; m < mmax; m += 2)
                {
                    for (int ii = m; ii <= 2*Nx; ii += istep)
                    {
                        jj = ii + mmax;
                        tempr = wr*A[jj-1][j][k] - wi*A[jj][j][k];
                        tempi = wr*A[jj][j][k] + wi*A[jj-1][j][k];
                         
                        A[jj-1][j][k] = A[ii-1][j][k]-tempr;
                        A[jj][j][k] = A[ii][j][k]-tempi;
                        A[ii-1][j][k] += tempr;
                        A[ii][j][k] += tempi;
                    }
                    
                    wtemp=wr;
                    wr += (wr*wpr - wi*wpi);
                    wi += (wi*wpr + wtemp*wpi);
                }
                
                mmax = istep;
            }
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Nx; ++i)
            for (int j{}; j < Ny; ++j)
                for (int k{}; k < Nz; ++k)
                    A[i][j][k] /= Nx;
        
        RealComplexX(A,Nx,Ny,Nz,0);
    }
}

//____________________________________________________________________________________
// 'isign' =  1  : Fourier transform
// 'isign' = -1  : Inverse Fourier transform
// 'real' = true : Sequence of real numbers
// 'real' = false: Sequence of complex numbers
// A             : Nx x 2*Ny x Nz
//____________________________________________________________________________________
void FourierTransformY ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Ny > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[i][Ny+j][k] = 0.0;
        
        RealComplexY(A,Nx,Ny,Nz,1);
    }
    
    for (int i{}; i < Nx; ++i)
    {
        for (int k{}; k < Nz; ++k)
        {
            int mmax = 2, m, jj = 1, istep, nn = Ny;
            double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
            
            // reverse-binary reindexing
            for (int ii = 1; ii < 2*Ny; ii += 2)
            {
                if (jj > ii)
                {
                    Swap(A[i][jj-1][k],A[i][ii-1][k]);
                    Swap(A[i][jj][k],A[i][ii][k]);
                }
                
                m = nn;
                while (m >= 2 && jj > m)
                {
                    jj -= m;
                    m /= 2;
                }
                
                jj += m;
            }
            
            // Danielson-Lanczos section
            while (2*Ny > mmax)
            {
                istep = mmax * 2;
                theta = isign * (-FFT_2PI/mmax);
                wtemp = sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi = sin(theta);
                wr = 1.0;
                wi = 0.0;
                
                for (int m = 1; m < mmax; m += 2)
                {
                    for (int ii = m; ii <= 2*Ny; ii += istep)
                    {
                        jj = ii + mmax;
                        tempr = wr*A[i][jj-1][k] - wi*A[i][jj][k];
                        tempi = wr*A[i][jj][k] + wi*A[i][jj-1][k];
                        
                        A[i][jj-1][k] = A[i][ii-1][k] - tempr;
                        A[i][jj][k] = A[i][ii][k] - tempi;
                        A[i][ii-1][k] += tempr;
                        A[i][ii][k] += tempi;
                    }
                    
                    wtemp=wr;
                    wr += (wr*wpr - wi*wpi);
                    wi += (wi*wpr + wtemp*wpi);
                }
                
                mmax = istep;
            }
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < 2*Ny; ++j)
                for (int k{}; k < Nz; ++k)
                    A[i][j][k] /= Ny;
        
        RealComplexY(A,Nx,Ny,Nz,0);
    }
}

//____________________________________________________________________________________
// 'isign' =  1   : Fourier transform
// 'isign' = -1   : Inverse Fourier transform
// 'real' = true  : Sequence of real numbers
// 'real' = false : Sequence of complex numbers
// A              : Nx x Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformZ ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nz)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Nz > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[i][j][Nz+k] = 0.0;
        
        RealComplexZ(A,Nx,Ny,Nz,1);
    }
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
        {
            int mmax = 2, m, jj = 1, istep, nn = Nz;
            double wtemp, wr, wi, wpr, wpi, tempr, tempi, theta;
            
            // reverse-binary reindexing
            for (int ii = 1; ii < 2*Nz; ii += 2)
            {
                if (jj > ii)
                {
                    Swap(A[i][j][jj-1],A[i][j][ii-1]);
                    Swap(A[i][j][jj],A[i][j][ii]);
                }
                
                m = nn;
                while (m >= 2 && jj > m)
                {
                    jj -= m;
                    m /= 2;
                }
                
                jj += m;
            }
            
            // Danielson-Lanczos section
            while (2*Nz > mmax)
            {
                istep = mmax * 2;
                theta = isign * (-FFT_2PI/mmax);
                wtemp = sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi = sin(theta);
                wr = 1.0;
                wi = 0.0;
                
                for (int m = 1; m < mmax; m += 2)
                {
                    for (int ii = m; ii <= 2*Nz; ii += istep)
                    {
                        jj = ii + mmax;
                        tempr = wr*A[i][j][jj-1] - wi*A[i][j][jj];
                        tempi = wr*A[i][j][jj] + wi*A[i][j][jj-1];
                        
                        A[i][j][jj-1] = A[i][j][ii-1] - tempr;
                        A[i][j][jj] = A[i][j][ii] - tempi;
                        A[i][j][ii-1] += tempr;
                        A[i][j][ii] += tempi;
                    }
                    
                    wtemp=wr;
                    wr += (wr*wpr - wi*wpi);
                    wi += (wi*wpr + wtemp*wpi);
                }
                
                mmax = istep;
            }
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
            for (int j{}; j < Ny; ++j)
                for (int k{}; k < 2*Nz; ++k)
                    A[i][j][k] /= Nz;
        
        RealComplexZ(A,Nx,Ny,Nz,0);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny
//____________________________________________________________________________________
void FourierTransformXY ( double **&A,
                          const int isign,
                          const int Nx,
                          const int Ny )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Maximum(Nx,Ny) > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                A[i][j] = A[2*i][j];
                A[i][Ny+j] = A[2*i+1][j];
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,false);
    }
    else
    {
        FourierTransformY(A,-1,Nx,Ny,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                A[2*i][j] = A[i][j];
                A[2*i+1][j] = A[i][Ny+j];
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny x Nz
//____________________________________________________________________________________
void FourierTransformXY ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Maximum(Nx,Ny) > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][Ny+j][k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformY(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][Ny+j][k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformXZ ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    exponent = log10((double)Nz)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Maximum(Nx,Nz) > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][j][Nz+k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : Nx x 2*Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformYZ ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    exponent = log10((double)Nz)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("N has to be 2^n for some integer n.");
    
    if (Maximum(Ny,Nz) > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformY(A,1,Nx,Ny,Nz);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[i][2*j][k];
                    A[i][j][Nz+k] = A[i][2*j+1][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int j{Ny-1}; j >= 0; j--)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][2*j][k] = A[i][j][k];
                    A[i][2*j+1][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformY(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformXYZ ( double ***&A,
                           const int isign,
                           const int Nx,
                           const int Ny,
                           const int Nz )
{
    if (FFT_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    double exponent = log10((double)Nx)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("Nx has to be 2^n for some integer n.");
    
    exponent = log10((double)Ny)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("Ny has to be 2^n for some integer n.");
    
    exponent = log10((double)Nz)/log10(2.0);
    
    if (floor(exponent) != ceil(exponent))
        ErrorMessage("Nz has to be 2^n for some integer n.");
    
    if (Maximum(Nx,Ny,Nz) > FFT_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][Ny+j][k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,Nz,false);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[i][2*j][k];
                    A[i][j][Nz+k] = A[i][2*j+1][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int j{Ny-1}; j >= 0; j--)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][2*j][k] = A[i][j][k];
                    A[i][2*j+1][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformY(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][Ny+j][k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// First derivative along first direction.
// Fx : Nx
// F  : Nx
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double *&Fx,
                               double *F,
                               const int Nx,
                               const double Lx )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
        FFT_Temporary0[i] = F[i];
    
    FourierTransformX (FFT_Temporary0,1,Nx);
    
    // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k).
    for (int i{}; i < Nx; ++i)
    {
        K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
        
        scalar = FFT_Temporary0[2*i];
        FFT_Temporary0[2*i] = -FFT_2PI*K*FFT_Temporary0[2*i+1];
        FFT_Temporary0[2*i+1] = FFT_2PI*K*scalar;
    }
    
    FourierTransformX (FFT_Temporary0,-1,Nx);
    
    for (int i{}; i < Nx; ++i) 
        Fx[i] = FFT_Temporary0[i]/Lx;
}

//____________________________________________________________________________________
// First derivative along first direction.
// Fx : Nx x Ny
// F  : Nx x Ny
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double **&Fx,
                               double **F,
                               const int Nx,
                               const int Ny,
                               const double Lx )
{
    double scalar, K;
    
    for (int j{}; j < Ny; ++j)
    {
        for (int i{}; i < Nx; ++i)
            FFT_Temporary0[i] = F[i][j];
        
        FourierTransformX (FFT_Temporary0,1,Nx);
        
        // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)
        for (int i{}; i < Nx; ++i)
        {
            K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
            
            scalar = FFT_Temporary0[2*i];
            FFT_Temporary0[2*i] = -FFT_2PI*K*FFT_Temporary0[2*i+1];
            FFT_Temporary0[2*i+1] = FFT_2PI*K*scalar;
        }
        
        FourierTransformX (FFT_Temporary0,-1,Nx);
        
        for (int i{}; i < Nx; ++i)
            Fx[i][j] = FFT_Temporary0[i]/Lx;
    }
}

//____________________________________________________________________________________
// First derivative along second direction.
// Fy : Nx x Ny
// F  : Nx x Ny
//____________________________________________________________________________________
void FirstDerivativeFourierY ( double **&Fy,
                               double **F,
                               const int Nx,
                               const int Ny,
                               const double Ly )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
            FFT_Temporary0[j] = F[i][j];
        
        FourierTransformX (FFT_Temporary0,1,Ny);
        
        // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)
        for (int j{}; j < Ny; ++j)
        {
            K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
            
            scalar = FFT_Temporary0[2*j];
            FFT_Temporary0[2*j] = -FFT_2PI*K*FFT_Temporary0[2*j+1];
            FFT_Temporary0[2*j+1] = FFT_2PI*K*scalar;
        }
        
        FourierTransformX (FFT_Temporary0,-1,Ny);
        
        for (int j{}; j < Ny; ++j) 
            Fy[i][j] = FFT_Temporary0[j]/Ly;
    }
}

//____________________________________________________________________________________
// First derivative along first direction.
// Fx : Nx x Ny x Nz
// F  : Nx x Ny x Nz
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double ***&Fx,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Lx )
{
    double scalar, K;
    
    for (int j{}; j < Ny; ++j)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
                FFT_Temporary0[i] = F[i][j][k];
            
            FourierTransformX (FFT_Temporary0,1,Nx);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)
            for (int i{}; i < Nx; ++i)
            {
                K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
                
                scalar = FFT_Temporary0[2*i];
                FFT_Temporary0[2*i] = -FFT_2PI*K*FFT_Temporary0[2*i+1];
                FFT_Temporary0[2*i+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFT_Temporary0,-1,Nx);
            
            for (int i{}; i < Nx; ++i) 
                Fx[i][j][k] = FFT_Temporary0[i]/Lx;
        }
    }
}

//____________________________________________________________________________________
// First derivative along second direction.
// Fy : Nx x Ny x Nz
// F  : Nx x Ny x Nz
//____________________________________________________________________________________
void FirstDerivativeFourierY ( double ***&Fy,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Ly )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
                FFT_Temporary0[j] = F[i][j][k];
            
            FourierTransformX(FFT_Temporary0,1,Ny);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)
            for (int j{}; j < Ny; ++j)
            {
                K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
                
                scalar = FFT_Temporary0[2*j];
                FFT_Temporary0[2*j] = -FFT_2PI*K*FFT_Temporary0[2*j+1];
                FFT_Temporary0[2*j+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFT_Temporary0,-1,Ny);
            
            for (int j{}; j < Ny; ++j) 
                Fy[i][j][k] = FFT_Temporary0[j]/Ly;
        }
    }
}

//____________________________________________________________________________________
// First derivative along third direction.
// Fz : Nx x Ny x Nz
// F  : Nx x Ny x Nz
//____________________________________________________________________________________
void FirstDerivativeFourierZ ( double ***&Fz,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Lz )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
                FFT_Temporary0[k] = F[i][j][k];
            
            FourierTransformX (FFT_Temporary0,1,Nz);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)
            for (int k{}; k < Nz; ++k)
            {
                K = (k < Nz/2 ? double(k) : double(k-Nz));
                
                scalar = FFT_Temporary0[2*k];
                FFT_Temporary0[2*k] = -FFT_2PI*K*FFT_Temporary0[2*k+1];
                FFT_Temporary0[2*k+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFT_Temporary0,-1,Nz);
            
            for (int k{}; k < Nz; ++k) 
                Fz[i][j][k] = FFT_Temporary0[k]/Lz;
        }
    }
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx
// F   : Nx
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double *&Fxx,
                                double *F,
                                const int Nx,
                                const double Lx )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
        FFT_Temporary0[i] = F[i];
    
    FourierTransformX (FFT_Temporary0,1,Nx);
    
    // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2
    for (int i{}; i < Nx; ++i)
    {
        K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
        
        FFT_Temporary0[2*i] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i];
        FFT_Temporary0[2*i+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i+1];
    }
    
    FourierTransformX (FFT_Temporary0,-1,Nx);
    
    for (int i{}; i < Nx; ++i) 
        Fxx[i] = FFT_Temporary0[i]/(Lx*Lx);
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx x Ny
// F   : Nx x Ny
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double **&Fxx,
                                double **F,
                                const int Nx,
                                const int Ny,
                                const double Lx )
{
    double K;
    
    for (int j{}; j < Ny; ++j)
    {
        for (int i{}; i < Nx; ++i)
            FFT_Temporary0[i] = F[i][j];
        
        FourierTransformX (FFT_Temporary0,1,Nx);
        
        // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2
        for (int i{}; i < Nx; ++i)
        {
            K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
            
            FFT_Temporary0[2*i] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i];
            FFT_Temporary0[2*i+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i+1];
        }
        
        FourierTransformX (FFT_Temporary0,-1,Nx);
        
        for (int i{}; i < Nx; ++i) 
            Fxx[i][j] = FFT_Temporary0[i]/(Lx*Lx);
    }
}

//____________________________________________________________________________________
// Second derivative along second direction.
// Fyy : Nx x Ny
// F   : Nx x Ny
//____________________________________________________________________________________
void SecondDerivativeFourierY ( double **&Fyy,
                                double **F,
                                const int Nx,
                                const int Ny,
                                const double Ly )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
            FFT_Temporary0[j] = F[i][j];
        
        FourierTransformX (FFT_Temporary0,1,Ny);
        
        // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2
        for (int j{}; j < Ny; ++j)
        {
            K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
            
            FFT_Temporary0[2*j] = -FFT_4PIPI*K*K*FFT_Temporary0[2*j];
            FFT_Temporary0[2*j+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*j+1];
        }
        
        FourierTransformX (FFT_Temporary0,-1,Ny);
        
        for (int j{}; j < Ny; ++j) 
            Fyy[i][j] = FFT_Temporary0[j]/(Ly*Ly);
    }
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double ***&Fxx,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Lx )
{
    double K;
    
    for (int j{}; j < Ny; ++j)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
                FFT_Temporary0[i] = F[i][j][k];
            
            FourierTransformX (FFT_Temporary0,1,Nx);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2
            for (int i{}; i < Nx; ++i)
            {
                K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
                
                FFT_Temporary0[2*i] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i];
                FFT_Temporary0[2*i+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*i+1];
            }
            
            FourierTransformX (FFT_Temporary0,-1,Nx);
            
            for (int i{}; i < Nx; ++i) 
                Fxx[i][j][k] = FFT_Temporary0[i]/(Lx*Lx);
        }
    }
}

//____________________________________________________________________________________
// Second derivative along second direction.
// Fyy : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierY ( double ***&Fyy,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Ly )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
                FFT_Temporary0[j] = F[i][j][k];
            
            FourierTransformX (FFT_Temporary0,1,Ny);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2
            for (int j{}; j < Ny; ++j)
            {
                K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
                
                FFT_Temporary0[2*j] = -FFT_4PIPI*K*K*FFT_Temporary0[2*j];
                FFT_Temporary0[2*j+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*j+1];
            }
            
            FourierTransformX (FFT_Temporary0,-1,Ny);
            
            for (int j{}; j < Ny; ++j) 
                Fyy[i][j][k] = FFT_Temporary0[j]/(Ly*Ly);
        }
    }
}

//____________________________________________________________________________________
// Second derivative along third direction.
// Fzz : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierZ ( double ***&Fzz,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Lz )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
                FFT_Temporary0[k] = F[i][j][k];
            
            FourierTransformX (FFT_Temporary0,1,Nz);
            
            // Rearrangement of FFT_Temporary0 and multiplication by (2*PI*i*k)^2.
            for (int k{}; k < Nz; ++k)
            {
                K = (k < Nz/2 ? (double)(k) : (double)(k - Nz));
                
                FFT_Temporary0[2*k] = -FFT_4PIPI*K*K*FFT_Temporary0[2*k];
                FFT_Temporary0[2*k+1] = -FFT_4PIPI*K*K*FFT_Temporary0[2*k+1];
            }
            
            FourierTransformX (FFT_Temporary0,-1,Nz);
            
            for (int k{}; k < Nz; ++k) 
                Fzz[i][j][k] = FFT_Temporary0[k]/(Lz*Lz);
        }
    }
}
#else
//________________________________________________________________________________________________
// Fourier quadratures using FFTW routine
//________________________________________________________________________________________________
#include <fftw3.h>

bool FFTWPlanned = false;

int FFTW_N = -1;

fftw_complex *FFTW_in, *FFTW_out;

double *FFTW_Temporary;

fftw_plan plan_backward, plan_forward;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreateFFTPlan ( const int N )
{
    if (!FFTWPlanned)
    {
        FFTW_N = N;
        FFTWPlanned = true;
        
        // Allocate 'FFTW_in'
        if ( ( FFTW_in = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)) ) == NULL )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        for (int i{}; i < N; ++i)
        {
            FFTW_in[i][0] = 0.0;
            FFTW_in[i][1] = 0.0;
        }
        
        // Allocate 'FFTW_out'
        if ( ( FFTW_out = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)) ) == NULL )
            ErrorMessage("Error: Memory cannot be allocated!");
        
        for (int i{}; i < N; ++i)
        {
            FFTW_out[i][0] = 0.0;
            FFTW_out[i][1] = 0.0;
        }
        
        // Allocate 'FFTW_Temporary'
        Allocate(FFTW_Temporary,2*N);
        
        // Create FFTW plans
        plan_forward  = fftw_plan_dft_1d(N,FFTW_in,FFTW_out,FFTW_FORWARD,FFTW_ESTIMATE);
        plan_backward = fftw_plan_dft_1d(N,FFTW_in,FFTW_out,FFTW_BACKWARD,FFTW_ESTIMATE);
    }
    else
    {
        if (N > FFTW_N)
        {
            // Deallocate all previously allocated memory first
            fftw_free(FFTW_in);
            fftw_free(FFTW_out);
            
            // Deallocate 'FFTW_Temporary'
            Deallocate(FFTW_Temporary,2*FFTW_N);
            
            // Reallocate 'FFTW_in'
            if ( ( FFTW_in = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)) ) == NULL )
                ErrorMessage("Error: Memory cannot be allocated!");
            
            for (int i{}; i < N; ++i)
            {
                FFTW_in[i][0] = 0.0;
                FFTW_in[i][1] = 0.0;
            }
            
            // Reallocate 'FFTW_out'
            if ( ( FFTW_out = (fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)) ) == NULL )
                ErrorMessage("Error: Memory cannot be allocated!");
            
            for (int i{}; i < N; ++i)
            {
                FFTW_out[i][0] = 0.0;
                FFTW_out[i][1] = 0.0;
            }
            
            // Reallocate 'FFTW_Temporary'
            Allocate(FFTW_Temporary,2*N);
            
            FFTW_N = N;
        }
        
        //ErrorMessage("One FFT plan already exists!");
    }
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyFFTPlan ( const int N )
{
    if (FFTWPlanned)
    {
        FFTWPlanned = false;
        
        fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);
        
        fftw_free(FFTW_in);
        fftw_free(FFTW_out);
        
        // Deallocate 'FFTW_Temporary'
        Deallocate(FFTW_Temporary,2*N);
        
        FFTW_N = -1;
    }
    else
        ErrorMessage("At first, create a FFT plan!");
}

//____________________________________________________________________________________
// This function returns the the Discrete Fourier Transform in 1D.
// 'isign' can be 1 or -1.
// A : 2*Nx 
//____________________________________________________________________________________
void DFT (  double *&A,
            const int isign,
            const int Nx )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
        RealComplexX(A,Nx,1);
    
    for (int k{}; k < 2*Nx; ++k)
        FFTW_Temporary[k] = 0.0;
    
    for (int k{}; k < Nx; ++k)
    {
        for (int j{}; j < Nx; ++j)
        {
            double C = cos(FFT_2PI*j*k/Nx);
            double S = sin(FFT_2PI*j*k/Nx);
            
            FFTW_Temporary[2*k]   += A[2*j]  *C + isign*A[2*j+1]*S;
            FFTW_Temporary[2*k+1] += A[2*j+1]*C - isign*A[2*j]  *S;
        }
    }
    
    if (isign == -1)
    {
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFTW_Temporary[i]/(double)Nx;
        
        RealComplexX(A,Nx,0);
    }
    else
    {
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFTW_Temporary[i];
    }
}

//____________________________________________________________________________________
// 'arrangement' = 0 means [real real ... complex complex ... ]
// 'arrangement' = 1 means [real complex real complex ... ]
//  A : 2*Nx
//____________________________________________________________________________________
void RealComplexX ( double *&A,
                    const int Nx,
                    const int arrangement )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (arrangement == 0)
    {
        for (int i{}; i < Nx; ++i)
        {
            FFTW_Temporary[i] = A[2*i];
            FFTW_Temporary[Nx+i] = A[2*i+1];
        }
        
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFTW_Temporary[i];
    }
    else
    {
        for (int i{}; i < Nx; ++i)
        {
            FFTW_Temporary[2*i] = A[i];
            FFTW_Temporary[2*i+1] = A[Nx+i];
        }
        
        for (int i{}; i < 2*Nx; ++i)
            A[i] = FFTW_Temporary[i];
    }
}

//____________________________________________________________________________________
// 'isign' =  1 : Fourier transform
// 'isign' = -1 : Inverse Fourier transform
//  A           : 2*Nx
//____________________________________________________________________________________
void FourierTransformX ( double *&A,
                         const int isign,
                         const int Nx,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    for (int k{}; k < 2*Nx; ++k)
        FFTW_Temporary[k] = 0.0;
    
    // Forward transform
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                A[Nx+i] = 0.0;
        
        for (int i{}; i < Nx; ++i)
        {
            FFTW_in[i][0] = A[i];
            FFTW_in[i][1] = A[Nx+i];
        }
        
        fftw_execute(plan_forward);
        
        for (int i{}; i < Nx; ++i)
        {
            A[2*i] = FFTW_out[i][0];
            A[2*i+1] = FFTW_out[i][1];
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
        {
            FFTW_in[i][0] = A[2*i];
            FFTW_in[i][1] = A[2*i+1];
        }
        
        fftw_execute(plan_backward);
        
        for (int i{}; i < Nx; ++i)
        {
            A[i] = FFTW_out[i][0]/Nx;
            A[Nx+i] = FFTW_out[i][1]/Nx;
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1 : Fourier transform
// 'isign' = -1 : Inverse Fourier transform
//  A           : 2*Nx x Ny
//____________________________________________________________________________________
void FourierTransformX ( double **&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    for (int k{}; k < 2*Nx; ++k)
        FFTW_Temporary[k] = 0.0;
    
    // Forward transform
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    A[Nx+i][j] = 0.0;
        
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                FFTW_in[i][0] = A[i][j];
                FFTW_in[i][1] = A[Nx+i][j];
            }
            
            fftw_execute(plan_forward);
            
            for (int i{}; i < Nx; ++i)
            {
                A[2*i][j] = FFTW_out[i][0];
                A[2*i+1][j] = FFTW_out[i][1];
            }
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                FFTW_in[i][0] = A[2*i][j];
                FFTW_in[i][1] = A[2*i+1][j];
            }
            
            fftw_execute(plan_backward);
            
            for (int i{}; i < Nx; ++i)
            {
                A[i][j] = FFTW_out[i][0]/Nx;
                A[Nx+i][j] = FFTW_out[i][1]/Nx;
            }
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1 : Fourier transform
// 'isign' = -1 : Inverse Fourier transform
//  A           : Nx x 2*Ny
//____________________________________________________________________________________
void FourierTransformY ( double **&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Ny > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    // Forward transform
    if (isign == 1)
    {
        if (real)
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    A[i][Ny+j] = 0.0;
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                FFTW_in[j][0] = A[i][j];
                FFTW_in[j][1] = A[i][Ny+j];
            }
            
            fftw_execute(plan_forward);
            
            for (int j{}; j < Ny; ++j)
            {
                A[i][2*j] = FFTW_out[j][0];
                A[i][2*j+1] = FFTW_out[j][1];
            }
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                FFTW_in[j][0] = A[i][2*j];
                FFTW_in[j][1] = A[i][2*j+1];
            }
            
            fftw_execute(plan_backward);
            
            for (int j{}; j < Ny; ++j)
            {
                A[i][j] = FFTW_out[j][0]/Ny;
                A[i][Ny+j] = FFTW_out[j][1]/Ny;
            }
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1 : Fourier transform
// 'isign' = -1 : Inverse Fourier transform
//  A           : 2*Nx x Ny x Nz
//____________________________________________________________________________________
void FourierTransformX ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nx > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    // Forward transform
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[Nx+i][j][k] = 0.0;
        
        for (int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int i{}; i < Nx; ++i)
                {
                    FFTW_in[i][0] = A[i][j][k];
                    FFTW_in[i][1] = A[Nx+i][j][k];
                }
                
                fftw_execute(plan_forward);
                
                for (int i{}; i < Nx; ++i)
                {
                    A[2*i][j][k] = FFTW_out[i][0];
                    A[2*i+1][j][k] = FFTW_out[i][1];
                }
            }
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int i{}; i < Nx; ++i)
                {
                    FFTW_in[i][0] = A[2*i][j][k];
                    FFTW_in[i][1] = A[2*i+1][j][k];
                }
                
                fftw_execute(plan_backward);
                
                for (int i{}; i < Nx; ++i)
                {
                    A[i][j][k] = FFTW_out[i][0]/Nx;
                    A[Nx+i][j][k] = FFTW_out[i][1]/Nx;
                }
            }
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1  : Fourier transform
// 'isign' = -1  : Inverse Fourier transform
//  A            : Nx x 2*Ny x Nz
//____________________________________________________________________________________
void FourierTransformY ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Ny > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    // Forward transform
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[i][Ny+j][k] = 0.0;
        
        for (int i{}; i < Nx; ++i)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    FFTW_in[j][0] = A[i][j][k];
                    FFTW_in[j][1] = A[i][Ny+j][k];
                }
                
                fftw_execute(plan_forward);
                
                for (int j{}; j < Ny; ++j)
                {
                    A[i][2*j][k] = FFTW_out[j][0];
                    A[i][2*j+1][k] = FFTW_out[j][1];
                }
            }
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    FFTW_in[j][0] = A[i][2*j][k];
                    FFTW_in[j][1] = A[i][2*j+1][k];
                }
                
                fftw_execute(plan_backward);
                
                for (int j{}; j < Ny; ++j)
                {
                    A[i][j][k] = FFTW_out[j][0]/Ny;
                    A[i][Ny+j][k] = FFTW_out[j][1]/Ny;
                }
            }
        }
    }
}

//____________________________________________________________________________________
// 'isign' =  1 : Fourier transform
// 'isign' = -1 : Inverse Fourier transform
// A            : Nx x Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformZ ( double ***&A,
                         const int isign,
                         const int Nx,
                         const int Ny,
                         const int Nz,
                         const bool real = true )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if (Nz > FFTW_N)
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    // Forward transform
    if (isign == 1)
    {
        if ( real )
            for (int i{}; i < Nx; ++i)
                for (int j{}; j < Ny; ++j)
                    for (int k{}; k < Nz; ++k)
                        A[i][j][Nz+k] = 0.0;
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    FFTW_in[k][0] = A[i][j][k];
                    FFTW_in[k][1] = A[i][j][Nz+k];
                }
                
                fftw_execute(plan_forward);
                
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][2*k] = FFTW_out[k][0];
                    A[i][j][2*k+1] = FFTW_out[k][1];
                }
            }
        }
    }
    
    // Inverse transform
    if (isign == -1)
    {
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    FFTW_in[k][0] = A[i][j][2*k];
                    FFTW_in[k][1] = A[i][j][2*k+1];
                }
                
                fftw_execute(plan_backward);
                
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = FFTW_out[k][0]/Nz;
                    A[i][j][Nz+k] = FFTW_out[k][1]/Nz;
                }
            }
        }
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny
//____________________________________________________________________________________
void FourierTransformXY ( double **&A,
                          const int isign,
                          const int Nx,
                          const int Ny )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if ( (Nx > FFTW_N) || (Ny > FFTW_N) )
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                A[i][j] = A[2*i][j];
                A[i][Ny+j] = A[2*i+1][j];
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,false);
    }
    else
    {
        FourierTransformY(A,-1,Nx,Ny,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                A[2*i][j] = A[i][j];
                A[2*i+1][j] = A[i][Ny+j];
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny x Nz
//____________________________________________________________________________________
void FourierTransformXY ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if ( (Nx > FFTW_N) || (Ny > FFTW_N) )
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][Ny+j][k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformY(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][Ny+j][k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformXZ ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if ( (Nx > FFTW_N) || (Nz > FFTW_N) )
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][j][Nz+k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : Nx x 2*Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformYZ ( double ***&A,
                          const int isign,
                          const int Nx,
                          const int Ny,
                          const int Nz )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if ( (Nz > FFTW_N) || (Ny > FFTW_N) )
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformY(A,1,Nx,Ny,Nz);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[i][2*j][k];
                    A[i][j][Nz+k] = A[i][2*j+1][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int j{Ny-1}; j >= 0; j--)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][2*j][k] = A[i][j][k];
                    A[i][2*j+1][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformY(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// Multidirectional Fourier Transform.
// A : 2*Nx x 2*Ny x 2*Nz
//____________________________________________________________________________________
void FourierTransformXYZ ( double ***&A,
                           const int isign,
                           const int Nx,
                           const int Ny,
                           const int Nz )
{
    if (FFTW_N == -1)
        ErrorMessage("Create a FFT plan first!");
    
    if ( (Nx > FFTW_N) || (Ny > FFTW_N) || (Nz > FFTW_N) )
        ErrorMessage("Computation of DFT is not possible with the present FFT plan!");
    
    if (isign == 1)
    {
        FourierTransformX(A,1,Nx,Ny,Nz);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int k{}; k < Nz; ++k)
            {
                for (int j{}; j < Ny; ++j)
                {
                    A[i][j][k] = A[2*i][j][k];
                    A[i][Ny+j][k] = A[2*i+1][j][k];
                }
            }
        }
        
        FourierTransformY(A,1,Nx,Ny,Nz,false);
        
        for (int i{}; i < Nx; ++i)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][j][k] = A[i][2*j][k];
                    A[i][j][Nz+k] = A[i][2*j+1][k];
                }
            }
        }
        
        FourierTransformZ(A,1,Nx,Ny,Nz,false);
    }
    else
    {
        FourierTransformZ(A,-1,Nx,Ny,Nz,false);
        
        for (int j{Ny-1}; j >= 0; j--)
        {
            for (int i{}; i < Nx; ++i)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[i][2*j][k] = A[i][j][k];
                    A[i][2*j+1][k] = A[i][j][Nz+k];
                }
            }
        }
        
        FourierTransformY(A,-1,Nx,Ny,Nz,false);
        
        for (int i{Nx-1}; i >= 0; i--)
        {
            for (int j{}; j < Ny; ++j)
            {
                for (int k{}; k < Nz; ++k)
                {
                    A[2*i][j][k] = A[i][j][k];
                    A[2*i+1][j][k] = A[i][Ny+j][k];
                }
            }
        }
        
        FourierTransformX(A,-1,Nx,Ny,Nz);
    }
}

//____________________________________________________________________________________
// First derivative along the first direction.
// Fx : Nx
// F  : Nx
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double *&Fx,
                               double *F,
                               const int Nx,
                               const double Lx )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
        FFTW_Temporary[i] = F[i];
    
    FourierTransformX (FFTW_Temporary,1,Nx);
    
    // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
    for (int i{}; i < Nx; ++i)
    {
        K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
        
        scalar = FFTW_Temporary[2*i];
        FFTW_Temporary[2*i] = -FFT_2PI*K*FFTW_Temporary[2*i+1];
        FFTW_Temporary[2*i+1] = FFT_2PI*K*scalar;
    }
    
    FourierTransformX (FFTW_Temporary,-1,Nx);
    
    for (int i{}; i < Nx; ++i) 
        Fx[i] = FFTW_Temporary[i]/Lx;
}

//____________________________________________________________________________________
// First derivative along first direction.
// Fx : Nx x Ny
// F  : Nx x Ny
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double **&Fx,
                               double **F,
                               const int Nx,
                               const int Ny,
                               const double Lx )
{
    double scalar, K;
    
    for (int j{}; j < Ny; ++j)
    {
        for (int i{}; i < Nx; ++i)
            FFTW_Temporary[i] = F[i][j];
        
        FourierTransformX (FFTW_Temporary,1,Nx);
        
        // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
        for (int i{}; i < Nx; ++i)
        {
            K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
            
            scalar = FFTW_Temporary[2*i];
            FFTW_Temporary[2*i] = -FFT_2PI*K*FFTW_Temporary[2*i+1];
            FFTW_Temporary[2*i+1] = FFT_2PI*K*scalar;
        }
        
        FourierTransformX (FFTW_Temporary,-1,Nx);
        
        for (int i{}; i < Nx; ++i)
            Fx[i][j] = FFTW_Temporary[i]/Lx;
    }
}

//____________________________________________________________________________________
// First derivative along second direction.
// Fy : Nx x Ny
// F  : Nx x Ny
//____________________________________________________________________________________
void FirstDerivativeFourierY ( double **&Fy,
                               double **F,
                               const int Nx,
                               const int Ny,
                               const double Ly )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
            FFTW_Temporary[j] = F[i][j];
        
        FourierTransformX (FFTW_Temporary,1,Ny);
        
        // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
        for (int j{}; j < Ny; ++j)
        {
            K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
            
            scalar = FFTW_Temporary[2*j];
            FFTW_Temporary[2*j] = -FFT_2PI*K*FFTW_Temporary[2*j+1];
            FFTW_Temporary[2*j+1] = FFT_2PI*K*scalar;
        }
        
        FourierTransformX (FFTW_Temporary,-1,Ny);
        
        for (int j{}; j < Ny; ++j) 
            Fy[i][j] = FFTW_Temporary[j]/Ly;
    }
}

//____________________________________________________________________________________
// First derivative along first direction.
// Fx : Nx x Ny x Nz
// F  : Nx x Ny x Nz
//____________________________________________________________________________________
void FirstDerivativeFourierX ( double ***&Fx,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Lx )
{
    double scalar, K;
    
    for (int j{}; j < Ny; ++j)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
                FFTW_Temporary[i] = F[i][j][k];
            
            FourierTransformX (FFTW_Temporary,1,Nx);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
            for (int i{}; i < Nx; ++i)
            {
                K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
                
                scalar = FFTW_Temporary[2*i];
                FFTW_Temporary[2*i] = -FFT_2PI*K*FFTW_Temporary[2*i+1];
                FFTW_Temporary[2*i+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFTW_Temporary,-1,Nx);
            
            for (int i{}; i < Nx; ++i) 
                Fx[i][j][k] = FFTW_Temporary[i]/Lx;
        }
    }
}

//____________________________________________________________________________________
// First derivative along second direction.
// Fy : Nx x Ny x Nz
// F  : Nx x Ny x Nz
//____________________________________________________________________________________
void FirstDerivativeFourierY ( double ***&Fy,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Ly )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
                FFTW_Temporary[j] = F[i][j][k];
            
            FourierTransformX(FFTW_Temporary,1,Ny);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
            for (int j{}; j < Ny; ++j)
            {
                K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
                
                scalar = FFTW_Temporary[2*j];
                FFTW_Temporary[2*j] = -FFT_2PI*K*FFTW_Temporary[2*j+1];
                FFTW_Temporary[2*j+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFTW_Temporary,-1,Ny);
            
            for (int j{}; j < Ny; ++j) 
                Fy[i][j][k] = FFTW_Temporary[j]/Ly;
        }
    }
}

//____________________________________________________________________________________
// First derivative along third direction.
// Fz         : Nx x Ny x Nz
// F          : Nx x Ny x Nz
// Temporary  : 2*Nz
//____________________________________________________________________________________
void FirstDerivativeFourierZ ( double ***&Fz,
                               double ***F,
                               const int Nx,
                               const int Ny,
                               const int Nz,
                               const double Lz )
{
    double scalar, K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
                FFTW_Temporary[k] = F[i][j][k];
            
            FourierTransformX (FFTW_Temporary,1,Nz);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k).
            for (int k{}; k < Nz; ++k)
            {
                K = (k < Nz/2 ? (double)(k) : (double)(k - Nz));
                
                scalar = FFTW_Temporary[2*k];
                FFTW_Temporary[2*k] = -FFT_2PI*K*FFTW_Temporary[2*k+1];
                FFTW_Temporary[2*k+1] = FFT_2PI*K*scalar;
            }
            
            FourierTransformX (FFTW_Temporary,-1,Nz);
            
            for (int k{}; k < Nz; ++k) 
                Fz[i][j][k] = FFTW_Temporary[k]/Lz;
        }
    }
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx
// F   : Nx
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double *&Fxx,
                                double *F,
                                const int Nx,
                                const double Lx )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
        FFTW_Temporary[i] = F[i];
    
    FourierTransformX (FFTW_Temporary,1,Nx);
    
    // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k)^2.
    for (int i{}; i < Nx; ++i)
    {
        K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
        
        FFTW_Temporary[2*i] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i];
        FFTW_Temporary[2*i+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i+1];
    }
    
    FourierTransformX (FFTW_Temporary,-1,Nx);
    
    for (int i{}; i < Nx; ++i) 
        Fxx[i] = FFTW_Temporary[i]/(Lx*Lx);
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx x Ny
// F   : Nx x Ny
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double **&Fxx,
                                double **F,
                                const int Nx,
                                const int Ny,
                                const double Lx )
{
    double K;
    
    for (int j{}; j < Ny; ++j)
    {
        for (int i{}; i < Nx; ++i)
            FFTW_Temporary[i] = F[i][j];
        
        FourierTransformX (FFTW_Temporary,1,Nx);
        
        // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k)^2.
        for (int i{}; i < Nx; ++i)
        {
            K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
            
            FFTW_Temporary[2*i] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i];
            FFTW_Temporary[2*i+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i+1];
        }
        
        FourierTransformX (FFTW_Temporary,-1,Nx);
        
        for (int i{}; i < Nx; ++i) 
            Fxx[i][j] = FFTW_Temporary[i]/(Lx*Lx);
    }
}

//____________________________________________________________________________________
// Second derivative along second direction.
// Fyy : Nx x Ny
// F   : Nx x Ny
//____________________________________________________________________________________
void SecondDerivativeFourierY ( double **&Fyy,
                                double **F,
                                const int Nx,
                                const int Ny,
                                const double Ly )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for (int j{}; j < Ny; ++j)
            FFTW_Temporary[j] = F[i][j];
        
        FourierTransformX (FFTW_Temporary,1,Ny);
        
        // Rearrangement of FFTW_Temporary and multiplication by (2*PI*i*k)^2.
        for (int j{}; j < Ny; ++j)
        {
            K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
            
            FFTW_Temporary[2*j] = -FFT_4PIPI*K*K*FFTW_Temporary[2*j];
            FFTW_Temporary[2*j+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*j+1];
        }
        
        FourierTransformX (FFTW_Temporary,-1,Ny);
        
        for (int j{}; j < Ny; ++j) 
            Fyy[i][j] = FFTW_Temporary[j]/(Ly*Ly);
    }
}

//____________________________________________________________________________________
// Second derivative along first direction.
// Fxx : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierX ( double ***&Fxx,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Lx )
{
    double K;
    
    for (int j{}; j < Ny; ++j)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int i{}; i < Nx; ++i)
                FFTW_Temporary[i] = F[i][j][k];
            
            FourierTransformX (FFTW_Temporary,1,Nx);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*i*PI*k)^2.
            for (int i{}; i < Nx; ++i)
            {
                K = (i < Nx/2 ? (double)(i) : (double)(i - Nx));
                
                FFTW_Temporary[2*i] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i];
                FFTW_Temporary[2*i+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*i+1];
            }
            
            FourierTransformX (FFTW_Temporary,-1,Nx);
            
            for (int i{}; i < Nx; ++i) 
                Fxx[i][j][k] = FFTW_Temporary[i]/(Lx*Lx);
        }
    }
}

//____________________________________________________________________________________
// Second derivative along second direction.
// Fyy : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierY ( double ***&Fyy,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Ly )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int k{}; k < Nz; ++k)
        {
            for (int j{}; j < Ny; ++j)
                FFTW_Temporary[j] = F[i][j][k];
            
            FourierTransformX (FFTW_Temporary,1,Ny);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*i*PI*k)^2.
            for (int j{}; j < Ny; ++j)
            {
                K = (j < Ny/2 ? (double)(j) : (double)(j - Ny));
                
                FFTW_Temporary[2*j] = -FFT_4PIPI*K*K*FFTW_Temporary[2*j];
                FFTW_Temporary[2*j+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*j+1];
            }
            
            FourierTransformX (FFTW_Temporary,-1,Ny);
            
            for (int j{}; j < Ny; ++j) 
                Fyy[i][j][k] = FFTW_Temporary[j]/(Ly*Ly);
        }
    }
}

//____________________________________________________________________________________
// Second derivative along third direction.
// Fzz : Nx x Ny x Nz
// F   : Nx x Ny x Nz
//____________________________________________________________________________________
void SecondDerivativeFourierZ ( double ***&Fzz,
                                double ***F,
                                const int Nx,
                                const int Ny,
                                const int Nz,
                                const double Lz )
{
    double K;
    
    for (int i{}; i < Nx; ++i)
    {
        for(int j{}; j < Ny; ++j)
        {
            for (int k{}; k < Nz; ++k)
                FFTW_Temporary[k] = F[i][j][k];
            
            FourierTransformX (FFTW_Temporary,1,Nz);
            
            // Rearrangement of FFTW_Temporary and multiplication by (2*i*PI*k)^2.
            for (int k{}; k < Nz; ++k)
            {
                K = (k < Nz/2 ? (double)(k) : (double)(k - Nz));
                
                FFTW_Temporary[2*k] = -FFT_4PIPI*K*K*FFTW_Temporary[2*k];
                FFTW_Temporary[2*k+1] = -FFT_4PIPI*K*K*FFTW_Temporary[2*k+1];
            }
            
            FourierTransformX (FFTW_Temporary,-1,Nz);
            
            for (int k{}; k < Nz; ++k) 
                Fzz[i][j][k] = FFTW_Temporary[k]/(Lz*Lz);
        }
    }
}
#endif
