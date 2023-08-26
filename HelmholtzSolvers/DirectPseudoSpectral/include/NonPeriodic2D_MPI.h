#pragma once

#ifndef NonPeriodic2D_MPI_H
#define NonPeriodic2D_MPI_H

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
	double **phi;
};

extern int Rank, size;

extern int Nx, Ny, NX, N, Ns;

extern double xl, xr, yl, yr, Lx, Ly, dx;

extern double alphal, betal, alphar, betar, Lambda;

extern double **Dy, **Dyy, **Operator, *y, ***rhsPDE, *leftBC, *rightBC, determinant;

extern int *npivot;

void CreatePDEPlan(const int Nx, const int Ny, const int N, const double xl, const double xr, const double yl, const double yr);
void DestroyPDEPlan();

void ConstructFundamentalSolutions(FieldVariable U[], const int i);
void SolvePDE(FieldVariable U[], void (*ConstructOperator)(const int) );

void FirstDerivativeX(FieldVariable Ux[], FieldVariable U[]);
void SecondDerivativeX(FieldVariable Uxx[], FieldVariable U[]);
void FirstDerivativeY(FieldVariable Uy[], FieldVariable U[]);
void SecondDerivativeY(FieldVariable Uyy[], FieldVariable U[]);

void ComputeError(FieldVariable U[], double (*ExactFunction)(const double, const double));
void WriteFile(FieldVariable U[], double (*ExactFunction)(const double, const double));
#endif
