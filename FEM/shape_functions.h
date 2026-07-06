#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

#include "general.h"

//______________________________________________________________________________
// Shape functions and reference-space gradients for the Gmsh element families
// 
//   1. P1, P2, P3 line elements
//   2. P1, P2, P3 quadrilateral elements
//   3. P1, P2, P3 hexahedral elements
//   4. P1, P2, P3 triangular elements
//   5. P1, P2, P3 tetrahedral elements
// 
// Node orderings match Gmsh MSH 2.2 native ordering.
// 
// Reference domains:
// 
//   1. Line         : [-1, 1]
//   2. Quadrilateral: [-1, 1] x [-1, 1]
//   3. Hexahedron   : [-1, 1] x [-1, 1] x [-1, 1]
//   4. Triangle     : {(zeta, eta) >= 0, (zeta + eta) <= 1}
//   5. Tetrahedron  : {(zeta, eta, xi) >= 0, (zeta + eta + xi) <= 1}
//   
// Conventions (say, number of nodes = n)
// 
//   1D boundary elements (line):
//       zeta         = zeta in [-1, 1] is the input
//       Psi[n]       = Psi[i] is the value of i-th nodal basis (Psi_i) at the input coordinate
//       dPsidzeta[n] = dPsidzeta[i] is the value of dPsi_i/dzeta at the input coordinate
// 
//   2D elements (triangle, quadrilateral):
//       zeta[2]          = (zeta, eta) in [-1, 1] x [-1, 1] is the input
//       Psi[n]           = Psi[i] is the value of i-th nodal basis (Psi_i) at the input coordinate
//       dPsidzeta[n][2]  = dPsidzeta[i][:] is [dPsi_i/dzeta, dPsi_i/deta] at the input coordinate
//
//   3D elements (tetrahedron, hexahedron):
//       zeta[3]          = (zeta, eta, xi) [-1, 1] x [-1, 1] x [-1, 1] is the input
//       Psi[n]           = Psi[i] is the value of i-th nodal basis (Psi_i) at the input coordinate
//       dPsidzeta[n][3]  = dPsidzeta[i][:] is [dPsi_i/dzeta, dPsi_i/deta, dPsi_i/dxi] at the input coordinate
//______________________________________________________________________________

//------------------------------------------------------------------------------
// 1D Lagrange polynomials at order 1, 2, 3 on [-1, 1].
//
// Order 1: nodes at zeta = -1, 1
// Order 2: nodes at zeta = -1, 0, 1
// Order 3: nodes at zeta = -1, -1/3, 1/3, 1
//
// The local node index is the same as the position index in this list.
//------------------------------------------------------------------------------
inline void lagrange1d_p1 ( const double zeta,
                            double       L[2],
                            double       dLdzeta[2] )
{
    const double t0 = 0.5*zeta;

    L[0]  = 0.5 - t0;
    L[1]  = 0.5 + t0;

    dLdzeta[0] = -0.5;
    dLdzeta[1] =  0.5;
}

inline void lagrange1d_p2 ( const double zeta,
                            double       L[3],
                            double       dLdzeta[3] )
{
    const double t0 = zeta*zeta;
    const double t1 = 0.5*t0;
    const double t2 = 0.5*zeta;

    L[0]  = t1 - t2;
    L[1]  = 1.0 - t0;
    L[2]  = t1 + t2;

    dLdzeta[0] = zeta - 0.5;
    dLdzeta[1] = -2.0 * zeta;
    dLdzeta[2] = zeta + 0.5;
}

inline void lagrange1d_p3 ( const double zeta,
                            double       L[4],
                            double       dLdzeta[4] )
{
    const double c0 = 0.0625;

    const double t0 = zeta*zeta;
    const double t1 = zeta*t0;
    const double t2 = 18.0*zeta;
    const double t3 = 27.0*zeta;
    const double t4 = 9.0*t0;
    const double t5 = 27.0*t0;
    const double t6 = 81.0*t0;
    const double t7 = 9.0*t1;
    const double t8 = 27.0*t1;

    L[0]  = c0 * (-t7 +  t4 + zeta -  1.0);
    L[1]  = c0 * ( t8 -  t4 - t3   +  9.0);
    L[2]  = c0 * (-t8 -  t4 + t3   +  9.0);
    L[3]  = c0 * ( t7 +  t4 - zeta -  1.0);

    dLdzeta[0] = c0 * (-t5 + t2 +  1.0);
    dLdzeta[1] = c0 * ( t6 - t2 - 27.0);
    dLdzeta[2] = c0 * (-t6 - t2 + 27.0);
    dLdzeta[3] = c0 * ( t5 + t2 -  1.0);
}

//------------------------------------------------------------------------------
// Line2 (Gmsh type 1)
//------------------------------------------------------------------------------
inline void shape_line2 ( const double zeta,
                          double       Psi[2],
                          double       dPsidzeta[2] )
{
    double L[2], dLdzeta[2];

    static const int node_sequence[2] = {0, 1};

    lagrange1d_p1(zeta, L, dLdzeta);

    for (int i = 0; i < 2; ++i)
    {
        Psi[i]       = L[node_sequence[i]];
        dPsidzeta[i] = dLdzeta[node_sequence[i]];
    }
}

//------------------------------------------------------------------------------
// Line3 (Gmsh type 8)
//------------------------------------------------------------------------------
inline void shape_line3 ( const double zeta,
                          double       Psi[3],
                          double       dPsidzeta[3] )
{
    double L[3], dLdzeta[3];

    static const int node_sequence[3] = {0, 2, 1};

    lagrange1d_p2(zeta, L, dLdzeta);

    for (int i = 0; i < 3; ++i)
    {
        Psi[i]       = L[node_sequence[i]];
        dPsidzeta[i] = dLdzeta[node_sequence[i]];
    }
}

//------------------------------------------------------------------------------
// Line4 (Gmsh type 26)
//------------------------------------------------------------------------------
inline void shape_line4 ( const double zeta,
                          double       Psi[4],
                          double       dPsidzeta[4] )
{
    double L[4], dLdzeta[4];

    static const int node_sequence[4] = {0, 3, 1, 2};

    lagrange1d_p3(zeta, L, dLdzeta);

    for (int i = 0; i < 4; ++i)
    {
        Psi[i]       = L[node_sequence[i]];
        dPsidzeta[i] = dLdzeta[node_sequence[i]];
    }
}

//------------------------------------------------------------------------------
// Quad4 (Gmsh type 3)
//------------------------------------------------------------------------------
inline void shape_quad4 ( const double zeta[2],
                          double       Psi[4],
                          double       dPsidzeta[4][2] )
{
    double Lzeta[2], dLzetadzeta[2];
    double Leta[2], dLetadeta[2];
    
    lagrange1d_p1(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p1(zeta[1], Leta, dLetadeta);

    static const int node_sequence[4][2] = {
        {0,0}, {1,0}, {1,1}, {0,1}
    };

    for (int node = 0; node < 4; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        
        Psi[node]          = Lzeta[i] * Leta[j];
        dPsidzeta[node][0] = dLzetadzeta[i] * Leta[j];
        dPsidzeta[node][1] = Lzeta[i]  * dLetadeta[j];
    }
}

//------------------------------------------------------------------------------
// Quad9 (Gmsh type 10)
//------------------------------------------------------------------------------
inline void shape_quad9 ( const double zeta[2],
                          double       Psi[9],
                          double       dPsidzeta[9][2] )
{
    double Lzeta[3], dLzetadzeta[3];
    double Leta[3], dLetadeta[3];
    
    lagrange1d_p2(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p2(zeta[1], Leta, dLetadeta);

    static const int node_sequence[9][2] = {
        {0,0}, {2,0}, {2,2}, 
        {0,2}, {1,0}, {2,1}, 
        {1,2}, {0,1}, {1,1}
    };

    for (int node = 0; node < 9; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        
        Psi[node]          = Lzeta[i] * Leta[j];
        dPsidzeta[node][0] = dLzetadzeta[i] * Leta[j];
        dPsidzeta[node][1] = Lzeta[i]  * dLetadeta[j];
    }
}

//------------------------------------------------------------------------------
// Quad16 (Gmsh type 36)
//------------------------------------------------------------------------------
inline void shape_quad16 ( const double zeta[2],
                           double       Psi[16],
                           double       dPsidzeta[16][2] )
{
    double Lzeta[4], dLzetadzeta[4];
    double Leta[4], dLetadeta[4];
    
    lagrange1d_p3(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p3(zeta[1], Leta, dLetadeta);

    static const int node_sequence[16][2] = {
        {0,0}, {3,0}, {3,3}, {0,3},
        {1,0}, {2,0}, {3,1}, {3,2},
        {2,3}, {1,3}, {0,2}, {0,1},
        {1,1}, {2,1}, {2,2}, {1,2}
    };

    for (int node = 0; node < 16; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        
        Psi[node]          = Lzeta[i] * Leta[j];
        dPsidzeta[node][0] = dLzetadzeta[i] * Leta[j];
        dPsidzeta[node][1] = Lzeta[i]  * dLetadeta[j];
    }
}

//------------------------------------------------------------------------------
// Hex8 (Gmsh type 5)
//------------------------------------------------------------------------------
inline void shape_hex8 ( const double zeta[3],
                         double       Psi[8],
                         double       dPsidzeta[8][3] )
{
    double Lzeta[2], dLzetadzeta[2];
    double Leta[2], dLetadeta[2];
    double Lxi[2], dLxidxi[2];
    
    lagrange1d_p1(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p1(zeta[1], Leta, dLetadeta);
    lagrange1d_p1(zeta[2], Lxi, dLxidxi);

    static const int node_sequence[8][3] = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };

    for (int node = 0; node < 8; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        const int k = node_sequence[node][2];
        
        Psi[node]          = Lzeta[i] * Leta[j] * Lxi[k];
        dPsidzeta[node][0] = dLzetadzeta[i]* Leta[j] * Lxi[k];
        dPsidzeta[node][1] = Lzeta[i] * dLetadeta[j]* Lxi[k];
        dPsidzeta[node][2] = Lzeta[i] * Leta[j] * dLxidxi[k];
    }
}

//------------------------------------------------------------------------------
// Hex27 (Gmsh type 12)
//------------------------------------------------------------------------------
inline void shape_hex27 ( const double zeta[3],
                          double       Psi[27],
                          double       dPsidzeta[27][3] )
{
    double Lzeta[3], dLzetadzeta[3];
    double Leta[3], dLetadeta[3];
    double Lxi[3], dLxidxi[3];
    
    lagrange1d_p2(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p2(zeta[1], Leta, dLetadeta);
    lagrange1d_p2(zeta[2], Lxi, dLxidxi);

    static const int node_sequence[27][3] = {
        {0,0,0}, {2,0,0}, {2,2,0}, 
        {0,2,0}, {0,0,2}, {2,0,2}, 
        {2,2,2}, {0,2,2}, {1,0,0}, 
        {0,1,0}, {0,0,1}, {2,1,0}, 
        {2,0,1}, {1,2,0}, {2,2,1}, 
        {0,2,1}, {1,0,2}, {0,1,2},
        {2,1,2}, {1,2,2}, {1,1,0}, 
        {1,0,1}, {0,1,1}, {2,1,1}, 
        {1,2,1}, {1,1,2}, {1,1,1}
    };

    for (int node = 0; node < 27; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        const int k = node_sequence[node][2];
        
        Psi[node]          = Lzeta[i] * Leta[j] * Lxi[k];
        dPsidzeta[node][0] = dLzetadzeta[i]* Leta[j] * Lxi[k];
        dPsidzeta[node][1] = Lzeta[i] * dLetadeta[j]* Lxi[k];
        dPsidzeta[node][2] = Lzeta[i] * Leta[j] * dLxidxi[k];
    }
}

//------------------------------------------------------------------------------
// Hex64 (Gmsh type 92)
//------------------------------------------------------------------------------
inline void shape_hex64 ( const double zeta[3],
                          double       Psi[64],
                          double       dPsidzeta[64][3] )
{
    double Lzeta[4], dLzetadzeta[4];
    double Leta[4], dLetadeta[4];
    double Lxi[4], dLxidxi[4];
    
    lagrange1d_p3(zeta[0], Lzeta, dLzetadzeta);
    lagrange1d_p3(zeta[1], Leta, dLetadeta);
    lagrange1d_p3(zeta[2], Lxi, dLxidxi);

    static const int node_sequence[64][3] = {
        {0,0,0}, {3,0,0}, {3,3,0}, {0,3,0},
        {0,0,3}, {3,0,3}, {3,3,3}, {0,3,3},
        {1,0,0}, {2,0,0}, {0,1,0}, {0,2,0},
        {0,0,1}, {0,0,2}, {3,1,0}, {3,2,0}, 
        {3,0,1}, {3,0,2}, {2,3,0}, {1,3,0}, 
        {3,3,1}, {3,3,2}, {0,3,1}, {0,3,2},
        {1,0,3}, {2,0,3}, {0,1,3}, {0,2,3},
        {3,1,3}, {3,2,3}, {2,3,3}, {1,3,3},
        {1,1,0}, {1,2,0}, {2,2,0}, {2,1,0},
        {1,0,1}, {2,0,1}, {2,0,2}, {1,0,2},
        {0,1,1}, {0,1,2}, {0,2,2}, {0,2,1},
        {3,1,1}, {3,2,1}, {3,2,2}, {3,1,2},
        {2,3,1}, {1,3,1}, {1,3,2}, {2,3,2},
        {1,1,3}, {2,1,3}, {2,2,3}, {1,2,3},
        {1,1,1}, {2,1,1}, {2,2,1}, {1,2,1},
        {1,1,2}, {2,1,2}, {2,2,2}, {1,2,2}
    };

    for (int node = 0; node < 64; node++)
    {
        const int i = node_sequence[node][0];
        const int j = node_sequence[node][1];
        const int k = node_sequence[node][2];
        
        Psi[node]          = Lzeta[i] * Leta[j] * Lxi[k];
        dPsidzeta[node][0] = dLzetadzeta[i]* Leta[j] * Lxi[k];
        dPsidzeta[node][1] = Lzeta[i] * dLetadeta[j]* Lxi[k];
        dPsidzeta[node][2] = Lzeta[i] * Leta[j] * dLxidxi[k];
    }
}

//------------------------------------------------------------------------------
// Tri3 (Gmsh type 2)
//------------------------------------------------------------------------------
inline void shape_tri3 ( const double zeta[2],
                         double       Psi[3],
                         double       dPsidzeta[3][2] )
{
    Psi[0] = 1.0 - zeta[0] - zeta[1];
    Psi[1] = zeta[0];
    Psi[2] = zeta[1];
    
    dPsidzeta[0][0] = -1.0;
    dPsidzeta[1][0] =  1.0;
    dPsidzeta[2][0] =  0.0;
    
    dPsidzeta[0][1] = -1.0;
    dPsidzeta[1][1] =  0.0;
    dPsidzeta[2][1] =  1.0;
}

//------------------------------------------------------------------------------
// Tri6 (Gmsh type 9)
//------------------------------------------------------------------------------
inline void shape_tri6 ( const double zeta[2],
                         double       Psi[6],
                         double       dPsidzeta[6][2] )
{
    // barycentric coordinate
    const double L0 = 1.0 - zeta[0] - zeta[1];
    const double L1 = zeta[0];
    const double L2 = zeta[1];

    const double t0 = 4.0*L0;
    const double t1 = 4.0*L1;
    const double t2 = 4.0*L2;
    
    Psi[0] = L0 * (2.0*L0 - 1.0);
    Psi[1] = L1 * (2.0*L1 - 1.0);
    Psi[2] = L2 * (2.0*L2 - 1.0);
    Psi[3] = t0 * L1;
    Psi[4] = t1 * L2;
    Psi[5] = t2 * L0;
    
    double dPsidL0[6];
    double dPsidL1[6];
    double dPsidL2[6];
    
    dPsidL0[0] = t0 - 1.0;
    dPsidL0[1] = 0.0;
    dPsidL0[2] = 0.0;
    dPsidL0[3] = t1;
    dPsidL0[4] = 0.0;
    dPsidL0[5] = t2;
    
    dPsidL1[0] = 0.0;
    dPsidL1[1] = t1 - 1.0;
    dPsidL1[2] = 0.0;
    dPsidL1[3] = t0;
    dPsidL1[4] = t2;
    dPsidL1[5] = 0.0;
    
    dPsidL2[0] = 0.0;
    dPsidL2[1] = 0.0;
    dPsidL2[2] = t2 - 1.0;
    dPsidL2[3] = 0.0;
    dPsidL2[4] = t1;
    dPsidL2[5] = t0;

    for (int i = 0; i < 6; i++)
    {
        // dPsi/dzeta = dPsi/dL1 - dPsi/dL0
        dPsidzeta[i][0] = dPsidL1[i] - dPsidL0[i];
        
        // dPsi/deta  = dPsi/dL2 - dPsi/dL0
        dPsidzeta[i][1] = dPsidL2[i] - dPsidL0[i];
    }
}

//------------------------------------------------------------------------------
// Tri10 (Gmsh type 21)
//------------------------------------------------------------------------------
inline void shape_tri10 ( const double zeta[2],
                          double       Psi[10],
                          double       dPsidzeta[10][2] )
{
    // barycentric coordinate
    const double L0 = 1.0 - zeta[0] - zeta[1];
    const double L1 = zeta[0];
    const double L2 = zeta[1];

    const double t0 = 3.0*L0;
    const double t1 = 3.0*L1;
    const double t2 = 3.0*L2;
    
    Psi[0] = 0.5 * L0 * (t0 - 1.0) * (t0 - 2.0);
    Psi[1] = 0.5 * L1 * (t1 - 1.0) * (t1 - 2.0);
    Psi[2] = 0.5 * L2 * (t2 - 1.0) * (t2 - 2.0);
    Psi[3] = 4.5 * L0 * L1 * (t0 - 1.0);
    Psi[4] = 4.5 * L0 * L1 * (t1 - 1.0);
    Psi[5] = 4.5 * L1 * L2 * (t1 - 1.0);
    Psi[6] = 4.5 * L1 * L2 * (t2 - 1.0);
    Psi[7] = 4.5 * L2 * L0 * (t2 - 1.0);
    Psi[8] = 4.5 * L2 * L0 * (t0 - 1.0);
    Psi[9] = t0 * t1 * t2;

    double dPsidL0[10];
    double dPsidL1[10];
    double dPsidL2[10];

    dPsidL0[0] = 0.5 * (3.0*t0*t0 - 6.0*t0 + 2.0);
    dPsidL0[1] = 0.0;
    dPsidL0[2] = 0.0;
    dPsidL0[3] = 9.0 * L1 * (t0 - 0.5);
    dPsidL0[4] = 4.5 * L1 * (t1 - 1.0);
    dPsidL0[5] = 0.0;
    dPsidL0[6] = 0.0;
    dPsidL0[7] = 4.5 * L2 * (t2 - 1.0);
    dPsidL0[8] = 9.0 * L2 * (t0 - 0.5);
    dPsidL0[9] = 3.0 * t1 * t2;
    
    dPsidL1[0] = 0.0;
    dPsidL1[1] = 0.5 * (3.0*t1*t1 - 6.0*t1 + 2.0);
    dPsidL1[2] = 0.0;
    dPsidL1[3] = 4.5 * L0 * (t0 - 1.0);
    dPsidL1[4] = 9.0 * L0 * (t1 - 0.5);
    dPsidL1[5] = 9.0 * L2 * (t1 - 0.5);
    dPsidL1[6] = 4.5 * L2 * (t2 - 1.0);
    dPsidL1[7] = 0.0;
    dPsidL1[8] = 0.0;
    dPsidL1[9] = 3.0 * t0 * t2;
    
    dPsidL2[0] = 0.0;
    dPsidL2[1] = 0.0;
    dPsidL2[2] = 0.5 * (3.0*t2*t2 - 6.0*t2 + 2.0);
    dPsidL2[3] = 0.0;
    dPsidL2[4] = 0.0;
    dPsidL2[5] = 4.5 * L1 * (t1 - 1.0);
    dPsidL2[6] = 9.0 * L1 * (t2 - 0.5);
    dPsidL2[7] = 9.0 * L0 * (t2 - 0.5);
    dPsidL2[8] = 4.5 * L0 * (t0 - 1.0);
    dPsidL2[9] = 3.0 * t0 * t1;

    for (int i = 0; i < 10; i++)
    {
        // dPsi/dzeta = dPsi/dL1 - dPsi/dL0
        dPsidzeta[i][0] = dPsidL1[i] - dPsidL0[i];
        
        // dPsi/deta  = dPsi/dL2 - dPsi/dL0
        dPsidzeta[i][1] = dPsidL2[i] - dPsidL0[i];
    }
}

//------------------------------------------------------------------------------
// Tet4 (Gmsh type 4)
//------------------------------------------------------------------------------
inline void shape_tet4 ( const double zeta[3],
                         double       Psi[4],
                         double       dPsidzeta[4][3] )
{
    Psi[0] = 1.0 - zeta[0] - zeta[1] - zeta[2];
    Psi[1] = zeta[0];
    Psi[2] = zeta[1];
    Psi[3] = zeta[2];

    dPsidzeta[0][0] = -1.0;
    dPsidzeta[1][0] =  1.0;
    dPsidzeta[2][0] =  0.0;
    dPsidzeta[3][0] =  0.0;
    
    dPsidzeta[0][1] = -1.0;
    dPsidzeta[1][1] =  0.0;
    dPsidzeta[2][1] =  1.0;
    dPsidzeta[3][1] =  0.0;
    
    dPsidzeta[0][2] = -1.0;
    dPsidzeta[1][2] =  0.0;
    dPsidzeta[2][2] =  0.0;
    dPsidzeta[3][2] =  1.0;
}

//------------------------------------------------------------------------------
// Tet10 (Gmsh type 11)
//------------------------------------------------------------------------------
inline void shape_tet10 ( const double zeta[3],
                          double       Psi[10],
                          double       dPsidzeta[10][3] )
{
    // barycentric coordinate
    const double L0 = 1.0 - zeta[0] - zeta[1] - zeta[2];
    const double L1 = zeta[0];
    const double L2 = zeta[1];
    const double L3 = zeta[2];

    const double t0 = 4.0*L0;
    const double t1 = 4.0*L1;
    const double t2 = 4.0*L2;
    const double t3 = 4.0*L3;
    
    Psi[0] = L0 * (2.0*L0 - 1.0);
    Psi[1] = L1 * (2.0*L1 - 1.0);
    Psi[2] = L2 * (2.0*L2 - 1.0);
    Psi[3] = L3 * (2.0*L3 - 1.0);
    Psi[4] = t0 * L1;
    Psi[5] = t1 * L2;
    Psi[6] = t2 * L0;
    Psi[7] = t3 * L0;
    Psi[8] = t3 * L2;
    Psi[9] = t3 * L1;

    double dPsidL0[10];
    double dPsidL1[10];
    double dPsidL2[10];
    double dPsidL3[10];
    
    dPsidL0[0] = t0 - 1.0;
    dPsidL0[1] = 0.0;
    dPsidL0[2] = 0.0;
    dPsidL0[3] = 0.0;
    dPsidL0[4] = t1;
    dPsidL0[5] = 0.0;
    dPsidL0[6] = t2;
    dPsidL0[7] = t3;
    dPsidL0[8] = 0.0;
    dPsidL0[9] = 0.0;
    
    dPsidL1[0] = 0.0;
    dPsidL1[1] = t1 - 1.0;
    dPsidL1[2] = 0.0;
    dPsidL1[3] = 0.0;
    dPsidL1[4] = t0;
    dPsidL1[5] = t2;
    dPsidL1[6] = 0.0;
    dPsidL1[7] = 0.0;
    dPsidL1[8] = 0.0;
    dPsidL1[9] = t3;
    
    dPsidL2[0] = 0.0;
    dPsidL2[1] = 0.0;
    dPsidL2[2] = t2 - 1.0;
    dPsidL2[3] = 0.0;
    dPsidL2[4] = 0.0;
    dPsidL2[5] = t1;
    dPsidL2[6] = t0;
    dPsidL2[7] = 0.0;
    dPsidL2[8] = t3;
    dPsidL2[9] = 0.0;
    
    dPsidL3[0] = 0.0;
    dPsidL3[1] = 0.0;
    dPsidL3[2] = 0.0;
    dPsidL3[3] = t3 - 1.0;
    dPsidL3[4] = 0.0;
    dPsidL3[5] = 0.0;
    dPsidL3[6] = 0.0;
    dPsidL3[7] = t0;
    dPsidL3[8] = t2;
    dPsidL3[9] = t1;
    
    for (int i = 0; i < 10; i++)
    {
        // dPsi/dzeta = dPsi/dL1 - dPsi/dL0
        dPsidzeta[i][0] = dPsidL1[i] - dPsidL0[i];
        
        // dPsi/deta = dPsi/dL2 - dPsi/dL0
        dPsidzeta[i][1] = dPsidL2[i] - dPsidL0[i];
        
        // dPsi/dxi = dPsi/dL3 - dPsi/dL0
        dPsidzeta[i][2] = dPsidL3[i] - dPsidL0[i];
    }
}

//------------------------------------------------------------------------------
// Tet20 (Gmsh type 29)
//------------------------------------------------------------------------------
inline void shape_tet20 ( const double zeta[3],
                          double       Psi[20],
                          double       dPsidzeta[20][3] )
{
    // barycentric coordinate
    const double L0 = 1.0 - zeta[0] - zeta[1] - zeta[2];
    const double L1 = zeta[0];
    const double L2 = zeta[1];
    const double L3 = zeta[2];

    const double t0 = 3.0*L0;
    const double t1 = 3.0*L1;
    const double t2 = 3.0*L2;
    const double t3 = 3.0*L3;
    
    Psi[0] = 0.5 * L0 * (t0 - 1.0) * (t0 - 2.0);
    Psi[1] = 0.5 * L1 * (t1 - 1.0) * (t1 - 2.0);
    Psi[2] = 0.5 * L2 * (t2 - 1.0) * (t2 - 2.0);
    Psi[3] = 0.5 * L3 * (t3 - 1.0) * (t3 - 2.0);
    Psi[4] = 4.5 * L0 * L1 * (t0 - 1.0);
    Psi[5] = 4.5 * L0 * L1 * (t1 - 1.0);
    Psi[6] = 4.5 * L1 * L2 * (t1 - 1.0);
    Psi[7] = 4.5 * L1 * L2 * (t2 - 1.0);
    Psi[8] = 4.5 * L2 * L0 * (t2 - 1.0);
    Psi[9] = 4.5 * L2 * L0 * (t0 - 1.0);
    Psi[10] = 4.5 * L3 * L0 * (t3 - 1.0);
    Psi[11] = 4.5 * L3 * L0 * (t0 - 1.0);
    Psi[12] = 4.5 * L3 * L2 * (t3 - 1.0);
    Psi[13] = 4.5 * L3 * L2 * (t2 - 1.0);
    Psi[14] = 4.5 * L3 * L1 * (t3 - 1.0);
    Psi[15] = 4.5 * L3 * L1 * (t1 - 1.0);
    Psi[16] = t0 * t1 * t2;
    Psi[17] = t0 * t1 * t3;
    Psi[18] = t0 * t2 * t3;
    Psi[19] = t1 * t2 * t3;

    double dPsidL0[20];
    double dPsidL1[20];
    double dPsidL2[20];
    double dPsidL3[20];

    dPsidL0[0] = 0.5 * (27.0*L0*L0 - 18.0*L0 + 2.0);
    dPsidL0[1] = 0.0;
    dPsidL0[2] = 0.0;
    dPsidL0[3] = 0.0;
    dPsidL0[4] = 4.5 * L1 * (6.0*L0 - 1.0);
    dPsidL0[5] = 4.5 * L1 * (3.0*L1 - 1.0);
    dPsidL0[6] = 0.0;
    dPsidL0[7] = 0.0;
    dPsidL0[8] = 4.5 * L2 * (3.0*L2 - 1.0);
    dPsidL0[9] = 4.5 * L2 * (6.0*L0 - 1.0);
    dPsidL0[10] = 4.5 * L3 * (3.0*L3 - 1.0);
    dPsidL0[11] = 4.5 * L3 * (6.0*L0 - 1.0);
    dPsidL0[12] = 0.0;
    dPsidL0[13] = 0.0;
    dPsidL0[14] = 0.0;
    dPsidL0[15] = 0.0;
    dPsidL0[16] = 27.0 * L1 * L2;
    dPsidL0[17] = 27.0 * L1 * L3;
    dPsidL0[18] = 27.0 * L2 * L3;
    dPsidL0[19] = 0.0;
    
    dPsidL1[0] = 0.0;
    dPsidL1[1] = 0.5 * (27.0*L1*L1 - 18.0*L1 + 2.0);
    dPsidL1[2] = 0.0;
    dPsidL1[3] = 0.0;
    dPsidL1[4] = 4.5 * L0 * (3.0*L0 - 1.0);
    dPsidL1[5] = 4.5 * L0 * (6.0*L1 - 1.0);
    dPsidL1[6] = 4.5 * L2 * (6.0*L1 - 1.0);
    dPsidL1[7] = 4.5 * L2 * (3.0*L2 - 1.0);
    dPsidL1[8] = 0.0;
    dPsidL1[9] = 0.0;
    dPsidL1[10] = 0.0;
    dPsidL1[11] = 0.0;
    dPsidL1[12] = 0.0;
    dPsidL1[13] = 0.0;
    dPsidL1[14] = 4.5 * L3 * (3.0*L3 - 1.0);
    dPsidL1[15] = 4.5 * L3 * (6.0*L1 - 1.0);
    dPsidL1[16] = 27.0 * L0 * L2;
    dPsidL1[17] = 27.0 * L0 * L3;
    dPsidL1[18] = 0.0;
    dPsidL1[19] = 27.0 * L2 * L3;
    
    dPsidL2[0] = 0.0;
    dPsidL2[1] = 0.0;
    dPsidL2[2] = 0.5 * (27.0*L2*L2 - 18.0*L2 + 2.0);
    dPsidL2[3] = 0.0;
    dPsidL2[4] = 0.0;
    dPsidL2[5] = 0.0;
    dPsidL2[6] = 4.5 * L1 * (3.0*L1 - 1.0);
    dPsidL2[7] = 4.5 * L1 * (6.0*L2 - 1.0);
    dPsidL2[8] = 4.5 * L0 * (6.0*L2 - 1.0);
    dPsidL2[9] = 4.5 * L0 * (3.0*L0 - 1.0);
    dPsidL2[10] = 0.0;
    dPsidL2[11] = 0.0;
    dPsidL2[12] = 4.5 * L3 * (3.0*L3 - 1.0);
    dPsidL2[13] = 4.5 * L3 * (6.0*L2 - 1.0);
    dPsidL2[14] = 0.0;
    dPsidL2[15] = 0.0;
    dPsidL2[16] = 27.0 * L0 * L1;
    dPsidL2[17] = 0.0;
    dPsidL2[18] = 27.0 * L0 * L3;
    dPsidL2[19] = 27.0 * L1 * L3;
    
    dPsidL3[0] = 0.0;
    dPsidL3[1] = 0.0;
    dPsidL3[2] = 0.0;
    dPsidL3[3] = 0.5 * (27.0*L3*L3 - 18.0*L3 + 2.0);
    dPsidL3[4] = 0.0;
    dPsidL3[5] = 0.0;
    dPsidL3[6] = 0.0;
    dPsidL3[7] = 0.0;
    dPsidL3[8] = 0.0;
    dPsidL3[9] = 0.0;
    dPsidL3[10] = 4.5 * L0 * (6.0*L3 - 1.0);
    dPsidL3[11] = 4.5 * L0 * (3.0*L0 - 1.0);
    dPsidL3[12] = 4.5 * L2 * (6.0*L3 - 1.0);
    dPsidL3[13] = 4.5 * L2 * (3.0*L2 - 1.0);
    dPsidL3[14] = 4.5 * L1 * (6.0*L3 - 1.0);
    dPsidL3[15] = 4.5 * L1 * (3.0*L1 - 1.0);
    dPsidL3[16] = 0.0;
    dPsidL3[17] = 27.0 * L0 * L1;
    dPsidL3[18] = 27.0 * L0 * L2;
    dPsidL3[19] = 27.0 * L1 * L2;

    for (int i = 0; i < 20; i++)
    {
        // dPsi/dzeta = dPsi/dL1 - dPsi/dL0
        dPsidzeta[i][0] = dPsidL1[i] - dPsidL0[i];
        
        // dPsi/deta = dPsi/dL2 - dPsi/dL0
        dPsidzeta[i][1] = dPsidL2[i] - dPsidL0[i];
        
        // dPsi/dxi = dPsi/dL3 - dPsi/dL0
        dPsidzeta[i][2] = dPsidL3[i] - dPsidL0[i];
    }
}
#endif
