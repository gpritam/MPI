#pragma once

#ifndef Fourier_Quadrature_H
#define Fourier_Quadrature_H

#include "General.h"

//_______________________________________________________________________________
//_______________________________________________________________________________
// C++ Library for the FFT, first and second drivative computation of periodic data.
// 
// Developed by: Dr. Pritam Giri
// Date : 23.01.2019
// Bangalore
//_______________________________________________________________________________
//_______________________________________________________________________________

//#define UseFFTWRoutines

//____________________________________________________________________________________
// 1. Input array for FFT (isign = 1) or, Output array for inverse of FFT (isign = -1)
// 
//         [ R(h(0)), R(h(1)), ... , R(h(N-1)), C(h(0)), C(h(1)), ..., C(h(N-1)) ]
// 
// 2. Output array for FFT (isign = 1) or, Input array for inverse of FFT (isign = -1)
// 
//         [ R(H(0)), C(H(0)), ..., R(H(N/2-1)), C(H(N/2-1)), R(H(-N/2)), C(H(-N/2)), ... , R(H(-1)), C(H(-1)) ]
// 
// 3. Arrangement of the wave numbers on the output array of this FFT routine
// 
//         0 1 2 3 ...... (N/2-2) (N/2-1) (-N/2) -(N/2-1) -(N/2-2) ....... -3 -2 -1
//____________________________________________________________________________________

//____________________________________________________________________________________
// Function prototypes
//____________________________________________________________________________________

const double FFT_2PI = 2.0*PI;
const double FFT_4PIPI = 4.0*PI*PI;

void DFT(double *&A, const int isign, const int Nx);

void CreateFFTPlan(const int N);
void DestroyFFTPlan(const int N);

void RealComplexX(double *&A,   const int Nx, const int arrangement);
void RealComplexX(double **&A,  const int Nx, const int Ny, const int arrangement);
void RealComplexY(double **&A,  const int Nx, const int Ny, const int arrangement);
void RealComplexX(double ***&A, const int Nx, const int Ny, const int Nz, const int arrangement);
void RealComplexY(double ***&A, const int Nx, const int Ny, const int Nz, const int arrangement);
void RealComplexZ(double ***&A, const int Nx, const int Ny, const int Nz, const int arrangement);

void FourierTransformX(double *&A,     const int isign, const int Nx, const bool real = true);
void FourierTransformX(double **&A,    const int isign, const int Nx, const int Ny, const bool real = true);
void FourierTransformY(double **&A,    const int isign, const int Nx, const int Ny, const bool real = true);
void FourierTransformX(double ***&A,   const int isign, const int Nx, const int Ny, const int Nz, const bool real = true);
void FourierTransformY(double ***&A,   const int isign, const int Nx, const int Ny, const int Nz, const bool real = true);
void FourierTransformZ(double ***&A,   const int isign, const int Nx, const int Ny, const int Nz, const bool real = true);
void FourierTransformXY(double **&A,   const int isign, const int Nx, const int Ny);
void FourierTransformXY(double ***&A,  const int isign, const int Nx, const int Ny, const int Nz);
void FourierTransformXZ(double ***&A,  const int isign, const int Nx, const int Ny, const int Nz);
void FourierTransformYZ(double ***&A,  const int isign, const int Nx, const int Ny, const int Nz);
void FourierTransformXYZ(double ***&A, const int isign, const int Nx, const int Ny, const int Nz);

void FirstDerivativeFourierX(double *&Fx,   double *F,   const int Nx, const double Lx);
void FirstDerivativeFourierX(double **&Fx,  double **F,  const int Nx, const int Ny, const double Lx);
void FirstDerivativeFourierY(double **&Fy,  double **F,  const int Nx, const int Ny, const double Ly);
void FirstDerivativeFourierX(double ***&Fx, double ***F, const int Nx, const int Ny, const int Nz, const double Lx);
void FirstDerivativeFourierY(double ***&Fy, double ***F, const int Nx, const int Ny, const int Nz, const double Ly);
void FirstDerivativeFourierZ(double ***&Fz, double ***F, const int Nx, const int Ny, const int Nz, const double Lz);

void SecondDerivativeFourierX(double *&Fxx,   double *F,   const int Nx, const double Lx);
void SecondDerivativeFourierX(double **&Fxx,  double **F,  const int Nx, const int Ny, const double Lx);
void SecondDerivativeFourierY(double **&Fyy,  double **F,  const int Nx, const int Ny, const double Ly);
void SecondDerivativeFourierX(double ***&Fxx, double ***F, const int Nx, const int Ny, const int Nz, const double Lx);
void SecondDerivativeFourierY(double ***&Fyy, double ***F, const int Nx, const int Ny, const int Nz, const double Ly);
void SecondDerivativeFourierZ(double ***&Fzz, double ***F, const int Nx, const int Ny, const int Nz, const double Lz);
#endif
