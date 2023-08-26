#pragma once

#ifndef Helmholtz3D_MPI_H
#define Helmholtz3D_MPI_H

#include "mpi.h"
#include "General.h"
#include "HypreSolvers.h"

#define TECPLOT
//#define VISIT

extern int Rank, size;

extern int Nx, Ny, Nz, Npx, Npy, Npz, offset;

extern double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, *x, *y, *z, dx, dy, dz;

extern double alphal, betal, alphar, betar, Lambda;
extern double alphab, betab, alphat, betat;
extern double alphaback, betaback, alphaf, betaf;

extern double **leftBC, **rightBC, **bottomBC, **topBC, **backBC, **frontBC, ***rhsPDE;

extern bool PeriodicBackFront, PeriodicBottomTop, PeriodicLeftRight;

void CreatePDEPlan(int Nx, int Ny, int Nz, double xl, double xr, double yl, double yr, double zl, double zr, bool PeriodicBackFront = false, bool PeriodicBottomTop = false, bool PeriodicLeftRight = false);
void DestroyPDEPlan();
void SolvePDE(double ***&U);
double FirstDerivativeX(double ***U, const int i, const int j, const int k);
double FirstDerivativeY(double ***U, const int i, const int j, const int k);
double FirstDerivativeZ(double ***U, const int i, const int j, const int k);
double SecondDerivativeX(double ***U, const int i, const int j, const int k);
double SecondDerivativeY(double ***U, const int i, const int j, const int k);
double SecondDerivativeZ(double ***U, const int i, const int j, const int k);
void ExchangeData(double ***&U);
void ExchangeData(int ***&U);
void ComputeError(double ***U, double (*ExactFunction)(const double, const double, const double));
void WriteFile(double ***U, double (*ExactFunction)(const double, const double, const double));
#endif
