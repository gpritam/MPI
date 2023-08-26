#pragma once

#ifndef Periodic3D_MPI_H
#define Periodic3D_MPI_H

#include "mpi.h"
#include "General.h"
#include "GLLQuadrature.h"
#include "LU_Decomposition.h"
#include "Fourier_Quadrature.h"
#include <bits/stdc++.h>

#define TECPLOT
//#define VISIT

struct FieldVariable
{
	double ***phi;
};

extern int Rank, size;

extern int Nx, Ny, Nz, NY, N, Ns;

extern double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, dx, dy, Lambda;

extern double **Dz, **Dzz, **Operator, *z, ****rhsPDE, determinant, **alphaPDE;

extern int *npivot;

void CreatePDEPlan(const int Nx, const int Ny, const int Nz, const int N, const double xl, const double xr, const double yl, const double yr, const double zl, const double zr);
void DestroyPDEPlan();

void ConstructFundamentalSolutions(FieldVariable U[], const int i, const int j);
void SolvePDE(FieldVariable U[], void (*ConstructOperator)(const int, const int));

void FirstDerivativeX(FieldVariable Ux[], FieldVariable U[]);
void SecondDerivativeX(FieldVariable Uxx[], FieldVariable U[]);
void FirstDerivativeY(FieldVariable Uy[], FieldVariable U[]);
void SecondDerivativeY(FieldVariable Uyy[], FieldVariable U[]);
void FirstDerivativeZ(FieldVariable Uz[], FieldVariable U[]);
void SecondDerivativeZ(FieldVariable Uzz[], FieldVariable U[]);

void ComputeError(FieldVariable U[], double (*ExactFunction)(const double, const double, const double));
void WriteFile(FieldVariable U[], double (*ExactFunction)(const double, const double, const double));
#endif
