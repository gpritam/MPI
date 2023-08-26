#pragma once

#ifndef HypreSolvers_H
#define HypreSolvers_H

#include <bits/stdc++.h>

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "General.h"

int hypre_FlexGMRESModifyPCAMGExample ( void *precond_data, 
                                        int iterations, 
                                        double rel_residual_norm );

void HypreSolve ( HYPRE_ParCSRMatrix AHypre, 
                  HYPRE_ParVector rhsHypre, 
                  HYPRE_ParVector solutionHypre, 
                  int &iterations, 
                  double &residual, 
                  const int solver_id = 0, 
                  const double tolerance = 1E-6, 
                  const int MaximumCycle = 20, 
                  const bool PrintSolutionInformation = false );

void PlotHypreMatrix ( HYPRE_IJMatrix AHypre, 
                       const int N );
#endif
