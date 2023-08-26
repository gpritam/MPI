#pragma once

#ifndef GLLQuadrature_H
#define GLLQuadrature_H

#include "General.h"
#include "SymmetricEigenSystems.h"

//________________________________________________________________________________________________
//________________________________________________________________________________________________
// This contains codes related to Gauss-Lobatto-Legendre quadrature.
// Developed by: Pritam Giri
// Date : 10.10.21
// IIT Kanpur
//________________________________________________________________________________________________
//________________________________________________________________________________________________

double Legendre(double x, const int n);
void GetLegendreValues(double **&LegendreValues, double *&x, const int N);
void GaussLobattoLegendrePoints(double *&x, const int N);
void GaussLobattoLegendreWeights(double *&w, double *x, const int N);
void PhysicalToSpectral(double *&H, double *h, double *x, double *w, const int N);
void SpectralToPhysical(double *&h, double *H, double *x, const int N);
void FirstDerivativeMatrix(double **&Dx, double *x, const int N);
void FirstDerivative(double *&hp, double *h, double **Dx, const int N);
double FirstDerivative(int k, double *h, double **Dx, const int N);
void SecondDerivativeMatrix(double **&Dxx, double **Dx, const int N);
void SecondDerivative(double *&hpp, double *h, double **Dxx, const int N);
double SecondDerivative(int k, double *h, double **Dxx, const int N);
double Integration(double *h, double *w, const int N);
double FirstInnerProduct(const int i, const int j);
double SecondInnerProduct(const int i, const int j);
#endif
