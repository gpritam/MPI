#pragma once

#ifndef LU_Decomposition_H
#define LU_Decomposition_H

#include "General.h"

void LUDecomposition(double **&A, int *&npivot, double &determinant, const int N);
void LUSolve(double **A, int *npivot, double *&b, const int N);
#endif
