#pragma once

#ifndef SymmetricEigenSystems_H
#define SymmetricEigenSystems_H

#include "General.h"

//_______________________________________________________________________________
// This contains codes related to computation of eigensystems of symmetric matrices.
// Developed by: Pritam Giri
// Date : 10.10.21
// IIT Kanpur
//_______________________________________________________________________________

void CheckSymmetry(double **A, const int N);
void Rotate ( double **&a, const double s, const double tau, const int i, const int j, const int k, const int l);
void JacobiMethodSymmetric(double **&A, double **&Q, double *&Lambda, const int N);
void HouseholderTransformationSymmetric(double **&A, double *&Ad, double *&As, const int N);
void FrancisMethodTridiagonalSymmetric(double *&Ad, double *As, double **&Q, const int N, const bool InitiallyTridiagonal = true);
void EigenSystemSymmetricMatrix(double **&A, double *&Lambda, const int N);
#endif
