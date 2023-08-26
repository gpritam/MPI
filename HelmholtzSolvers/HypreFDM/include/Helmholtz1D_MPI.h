#pragma once

#ifndef Helmholtz1D_MPI_H
#define Helmholtz1D_MPI_H

#include "mpi.h"
#include "General.h"
#include "HypreSolvers.h"

#define TECPLOT
//#define VISIT

extern int Rank, size;

extern int Nx, offset;

extern double xl, xr, Lx, *x, dx;

extern double alphal, betal, alphar, betar, Lambda;

extern double leftBC, rightBC, *rhsPDE;

extern bool PeriodicLeftRight;

void CreatePDEPlan(int Nx, double xl, double xr, bool PeriodicLeftRight = false);
void DestroyPDEPlan();
void SolvePDE(double *&U);
double FirstDerivativeX(double *U, const int i);
double SecondDerivativeX(double *U, const int i);
void ExchangeData(double *&U);
void ComputeError(double *U, double (*ExactFunction)(const double));
void WriteFile(double *U, double (*ExactFunction)(const double));
#endif
