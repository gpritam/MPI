#pragma once

#ifndef Helmholtz2D_MPI_H
#define Helmholtz2D_MPI_H

#include "mpi.h"
#include "General.h"
#include "HypreSolvers.h"

#define TECPLOT
//#define VISIT

extern int Rank, size;

extern int Nx, Ny, Npx, Npy, offset;

extern double xl, xr, yl, yr, Lx, Ly, *x, *y, dx, dy;

extern double alphal, betal, alphar, betar, Lambda;
extern double alphab, betab, alphat, betat;

extern double *leftBC, *rightBC, *bottomBC, *topBC, **rhsPDE;

extern bool PeriodicBottomTop, PeriodicLeftRight;

void CreatePDEPlan(int Nx, int Ny, double xl, double xr, double yl, double yr, bool PeriodicBottomTop = false, bool PeriodicLeftRight = false);
void DestroyPDEPlan();
void SolvePDE(double **&U);
double FirstDerivativeX(double **U, const int i, const int j);
double FirstDerivativeY(double **U, const int i, const int j);
double SecondDerivativeX(double **U, const int i, const int j);
double SecondDerivativeY(double **U, const int i, const int j);
void ExchangeData(double **&U);
void ExchangeData(int **&U);
void ComputeError(double **U, double (*ExactFunction)(const double, const double));
void WriteFile(double **U, double (*ExactFunction)(const double, const double));
#endif
