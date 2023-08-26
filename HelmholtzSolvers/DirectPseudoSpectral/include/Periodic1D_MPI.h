#pragma once

#ifndef Periodic1D_MPI_H
#define Periodic1D_MPI_H

#include "mpi.h"
#include "General.h"
#include "GLLQuadrature.h"
#include "LU_Decomposition.h"
#include <bits/stdc++.h>

#define TECPLOT
//#define VISIT

struct FieldVariable
{
	double *phi;
};

extern int Rank, size;

extern int Nx, N, Ns;

extern double xl, xr, Lambda, Lx;

extern double **Dx, **Dxx, **Operator, *x, **rhsPDE, determinant, alphaPDE;

extern int *npivot;

void CreatePDEPlan(const int Nx, const int N, const double xl, const double xr);
void DestroyPDEPlan();

void ConstructFundamentalSolutions(FieldVariable U[]);
void SolvePDE(FieldVariable U[], void (*ConstructOperator)());

void FirstDerivativeX(FieldVariable Ux[], FieldVariable U[]);
void SecondDerivativeX(FieldVariable Uxx[], FieldVariable U[]);

void ComputeError(FieldVariable U[], double (*ExactFunction)(const double));
void WriteFile(FieldVariable U[], double (*ExactFunction)(const double));
#endif
