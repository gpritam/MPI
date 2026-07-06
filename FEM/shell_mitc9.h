#ifndef SHELL_MITC9_H
#define SHELL_MITC9_H

//______________________________________________________________________________
// MITC9 degenerated shell element kernel.
//
// 9-node biquadratic quadrilateral shell, 5 DOF per node:
//     (u_x, u_y, u_z, alpha, beta)
// where alpha, beta are rotations about the nodal triad vectors V1, V2.
//
// Formulation: continuum-based (degenerated solid) shell, Bathe-style, with
// Mixed Interpolation of Tensorial Components (MITC9, Bucalem & Bathe 1993):
//   - covariant in-plane strains   e_rr, e_ss, e_rs  : assumed/tied
//   - covariant transverse shears  e_rt, e_st        : assumed/tied
// Tying points and assumed-strain interpolations follow the standard MITC9
// scheme; the tied covariant strains are then rotated into a local Cartesian
// lamina frame and combined with a plane-stress shell constitutive matrix.
//
// Depends only on: general.h, quadrature.h, shape_functions.h.
//
// Each benchmark driver includes this header (single translation unit per
// program), so file-scope state is declared static / functions inline.
//______________________________________________________________________________

#include <cmath>

#include "general.h"
#include "quadrature.h"
#include "shape_functions.h"

//==============================================================================
// 1.  Small fixed-size linear algebra helpers (plain double[3] / double[3][3]).
//==============================================================================

inline double v3dot ( const double a[3], const double b[3] )
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void v3cross ( const double a[3], const double b[3], double c[3] )
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double v3norm ( const double a[3] )
{
    return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

inline void v3normalize ( double a[3] )
{
    const double n = v3norm(a);
    if (n < 1.0e-30) print_error_message("v3normalize: zero-length vector.");
    a[0] /= n;  a[1] /= n;  a[2] /= n;
}

inline void v3copy ( const double a[3], double b[3] )
{
    b[0] = a[0];  b[1] = a[1];  b[2] = a[2];
}

// Determinant of a 3x3 matrix stored row-major as M[3][3].
inline double m3det ( const double M[3][3] )
{
    return   M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
           - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
           + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
}

// Inverse of a 3x3 matrix. Returns the determinant; aborts if singular.
inline double m3inv ( const double M[3][3], double Inv[3][3] )
{
    const double det = m3det(M);
    if (std::fabs(det) < 1.0e-30)
        print_error_message("m3inv: singular 3x3 matrix.");
    const double id = 1.0 / det;

    Inv[0][0] =  (M[1][1]*M[2][2] - M[1][2]*M[2][1]) * id;
    Inv[0][1] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) * id;
    Inv[0][2] =  (M[0][1]*M[1][2] - M[0][2]*M[1][1]) * id;
    Inv[1][0] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) * id;
    Inv[1][1] =  (M[0][0]*M[2][2] - M[0][2]*M[2][0]) * id;
    Inv[1][2] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) * id;
    Inv[2][0] =  (M[1][0]*M[2][1] - M[1][1]*M[2][0]) * id;
    Inv[2][1] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) * id;
    Inv[2][2] =  (M[0][0]*M[1][1] - M[0][1]*M[1][0]) * id;
    return det;
}

//==============================================================================
// 2.  Shell model data.
//
// A driver fills these arrays (via shell_allocate + direct writes) and then
// calls shell_build_directors() once before assembly.
//==============================================================================

namespace shell
{
    static int      NN   = 0;       // number of nodes
    static int      NE   = 0;       // number of elements
    static int      NDOF = 0;       // total DOF = 5 * NN

    static double **X     = nullptr;   // NN x 3 nodal coordinates
    static int    **conn  = nullptr;   // NE x 9 connectivity (0-based, quad9 order)
    static double  *thick = nullptr;   // NN nodal shell thickness

    static double **Vn = nullptr;      // NN x 3 nodal director (unit)
    static double **V1 = nullptr;      // NN x 3 nodal triad vector 1 (unit)
    static double **V2 = nullptr;      // NN x 3 nodal triad vector 2 (unit)

    static double  E  = 0.0;           // Young's modulus
    static double  nu = 0.0;           // Poisson ratio
    static double  shear_correction = 5.0 / 6.0;

    static const int NPE  = 9;         // nodes per element
    static const int NEDOF = 45;       // DOF per element = 5 * 9

    // Map (node, component) -> global DOF.  component: 0,1,2 = u; 3 = alpha; 4 = beta.
    inline int dof ( const int node, const int comp ) { return 5 * node + comp; }
}

// Allocate the model arrays once NN and NE are known.
inline void shell_allocate ( const int nn, const int ne )
{
    shell::NN   = nn;
    shell::NE   = ne;
    shell::NDOF = 5 * nn;

    allocate(shell::X,     nn, 3);
    allocate(shell::conn,  ne, shell::NPE);
    allocate(shell::thick, nn);
    allocate(shell::Vn,    nn, 3);
    allocate(shell::V1,    nn, 3);
    allocate(shell::V2,    nn, 3);
}

inline void shell_free ()
{
    deallocate(shell::X,     shell::NN, 3);
    deallocate(shell::conn,  shell::NE, shell::NPE);
    deallocate(shell::thick, shell::NN);
    deallocate(shell::Vn,    shell::NN, 3);
    deallocate(shell::V1,    shell::NN, 3);
    deallocate(shell::V2,    shell::NN, 3);
}

//==============================================================================
// 3.  Nodal director triads.
//
// The director Vn at a node is the area-weighted average of the unit normals
// of all elements sharing the node (each element normal from g_r x g_s at the
// node's reference location).  V1, V2 complete a right-handed orthonormal
// triad (V1 x V2 = Vn), chosen stably with respect to the global axes.
//==============================================================================

// Reference (r,s) coordinates of the 9 quad9 nodes, matching shape_quad9's
// lattice ordering: corners, edge mids, centre.
inline void shell_quad9_node_rs ( const int a, double &r, double &s )
{
    static const double R[9] = { -1,  1,  1, -1,  0,  1,  0, -1,  0 };
    static const double S[9] = { -1, -1,  1,  1, -1,  0,  1,  0,  0 };
    r = R[a];  s = S[a];
}

inline void shell_build_directors ()
{
    for (int n = 0; n < shell::NN; n++)
        for (int k = 0; k < 3; k++)
            shell::Vn[n][k] = 0.0;

    // Accumulate element normals at each node.
    for (int e = 0; e < shell::NE; e++)
    {
        for (int a = 0; a < shell::NPE; a++)
        {
            double r, s;
            shell_quad9_node_rs(a, r, s);
            const double xi[2] = { r, s };

            double N[9], dN[9][2];
            shape_quad9(xi, N, dN);

            double gr[3] = {0,0,0}, gs[3] = {0,0,0};
            for (int b = 0; b < shell::NPE; b++)
            {
                const int nb = shell::conn[e][b];
                for (int k = 0; k < 3; k++)
                {
                    gr[k] += dN[b][0] * shell::X[nb][k];
                    gs[k] += dN[b][1] * shell::X[nb][k];
                }
            }
            double nrm[3];
            v3cross(gr, gs, nrm);
            v3normalize(nrm);

            const int na = shell::conn[e][a];
            for (int k = 0; k < 3; k++)
                shell::Vn[na][k] += nrm[k];
        }
    }

    // Normalize directors and build the complementary triad.
    for (int n = 0; n < shell::NN; n++)
    {
        v3normalize(shell::Vn[n]);

        // Pick the global axis least aligned with Vn to seed V1.
        const double ax = std::fabs(shell::Vn[n][0]);
        const double ay = std::fabs(shell::Vn[n][1]);
        const double az = std::fabs(shell::Vn[n][2]);
        double seed[3] = {0,0,0};
        if      (ax <= ay && ax <= az) seed[0] = 1.0;
        else if (ay <= ax && ay <= az) seed[1] = 1.0;
        else                           seed[2] = 1.0;

        // V1 = normalize(seed x Vn) ... actually Vn x seed, then orthonormalize.
        double v1[3];
        v3cross(seed, shell::Vn[n], v1);
        v3normalize(v1);
        double v2[3];
        v3cross(shell::Vn[n], v1, v2);   // V2 = Vn x V1  => V1 x V2 = Vn
        v3normalize(v2);

        v3copy(v1, shell::V1[n]);
        v3copy(v2, shell::V2[n]);
    }
}

//==============================================================================
// 4.  Constitutive matrix and covariant strain-displacement operator.
//==============================================================================

// Plane-stress shell constitutive matrix (5x5) for engineering strains
//   [ eps11, eps22, gamma12, gamma13, gamma23 ].
inline void shell_constitutive ( double D[5][5] )
{
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            D[i][j] = 0.0;

    const double f  = shell::E / (1.0 - shell::nu * shell::nu);
    const double g  = 0.5 * (1.0 - shell::nu) * f;          // in-plane shear
    const double gs = shell::shear_correction * g;          // transverse shear

    D[0][0] = f;             D[0][1] = shell::nu * f;
    D[1][0] = shell::nu * f; D[1][1] = f;
    D[2][2] = g;
    D[3][3] = gs;
    D[4][4] = gs;
}

// Covariant base vectors g_r, g_s, g_t of element e at natural point (r,s,t).
inline void shell_base_vectors ( const int e, const double r, const double s,
                                 const double t,
                                 double gr[3], double gs[3], double gt[3] )
{
    const double xi[2] = { r, s };
    double N[9], dN[9][2];
    shape_quad9(xi, N, dN);

    for (int k = 0; k < 3; k++) { gr[k] = gs[k] = gt[k] = 0.0; }

    for (int b = 0; b < shell::NPE; b++)
    {
        const int    nb   = shell::conn[e][b];
        const double half = 0.5 * shell::thick[nb];
        for (int k = 0; k < 3; k++)
        {
            const double xfib = shell::X[nb][k] + t * half * shell::Vn[nb][k];
            gr[k] += dN[b][0] * xfib;
            gs[k] += dN[b][1] * xfib;
            gt[k] += N[b] * half * shell::Vn[nb][k];
        }
    }
}

// Covariant strain-displacement rows at natural point (r,s,t) of element e.
//   B[0..4][0..44]  : tensor covariant components  e_rr, e_ss, e_rs, e_rt, e_st
// Also returns the covariant base vectors used (needed for the lamina rotation).
inline void shell_covariant_B ( const int e, const double r, const double s,
                                const double t,
                                double B[5][45],
                                double gr[3], double gs[3], double gt[3] )
{
    const double xi[2] = { r, s };
    double N[9], dN[9][2];
    shape_quad9(xi, N, dN);

    shell_base_vectors(e, r, s, t, gr, gs, gt);

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 45; j++)
            B[i][j] = 0.0;

    for (int b = 0; b < shell::NPE; b++)
    {
        const int    nb   = shell::conn[e][b];
        const double half = 0.5 * shell::thick[nb];
        const double nr   = dN[b][0];
        const double ns   = dN[b][1];
        const double nn   = N[b];

        const double *V1 = shell::V1[nb];
        const double *V2 = shell::V2[nb];

        // Coefficient vectors of the rotational DOFs in u_,r, u_,s, u_,t.
        double cAr[3], cBr[3], cAs[3], cBs[3], cAt[3], cBt[3];
        for (int k = 0; k < 3; k++)
        {
            cAr[k] = nr * t * half * (-V2[k]);
            cBr[k] = nr * t * half * ( V1[k]);
            cAs[k] = ns * t * half * (-V2[k]);
            cBs[k] = ns * t * half * ( V1[k]);
            cAt[k] = nn * half * (-V2[k]);
            cBt[k] = nn * half * ( V1[k]);
        }

        const int d0 = 5 * b;          // local DOF block start for this node

        // e_rr = g_r . u_,r
        for (int k = 0; k < 3; k++) B[0][d0+k] += nr * gr[k];
        B[0][d0+3] += v3dot(gr, cAr);
        B[0][d0+4] += v3dot(gr, cBr);

        // e_ss = g_s . u_,s
        for (int k = 0; k < 3; k++) B[1][d0+k] += ns * gs[k];
        B[1][d0+3] += v3dot(gs, cAs);
        B[1][d0+4] += v3dot(gs, cBs);

        // e_rs = 1/2 (g_r . u_,s + g_s . u_,r)
        for (int k = 0; k < 3; k++)
            B[2][d0+k] += 0.5 * (ns * gr[k] + nr * gs[k]);
        B[2][d0+3] += 0.5 * (v3dot(gr, cAs) + v3dot(gs, cAr));
        B[2][d0+4] += 0.5 * (v3dot(gr, cBs) + v3dot(gs, cBr));

        // e_rt = 1/2 (g_r . u_,t + g_t . u_,r)
        for (int k = 0; k < 3; k++)
            B[3][d0+k] += 0.5 * nr * gt[k];
        B[3][d0+3] += 0.5 * (v3dot(gr, cAt) + v3dot(gt, cAr));
        B[3][d0+4] += 0.5 * (v3dot(gr, cBt) + v3dot(gt, cBr));

        // e_st = 1/2 (g_s . u_,t + g_t . u_,s)
        for (int k = 0; k < 3; k++)
            B[4][d0+k] += 0.5 * ns * gt[k];
        B[4][d0+3] += 0.5 * (v3dot(gs, cAt) + v3dot(gt, cAs));
        B[4][d0+4] += 0.5 * (v3dot(gs, cBt) + v3dot(gt, cBs));
    }
}

//==============================================================================
// 5.  MITC9 assumed-strain interpolation and element stiffness.
//
// All five strain components are re-interpolated from tying points to remove
// membrane locking AND transverse-shear locking.  Two ingredients make the
// element pass the distorted-mesh patch test:
//
//  (1) Fixed-frame tying.  The strains are transformed into a single, fixed
//      element-centre Cartesian frame BEFORE the tying interpolation.  A
//      constant strain tensor then has constant components there, which any
//      partition-of-unity tying reproduces exactly.
//
//  (2) Mean correction.  Pointwise constant-strain reproduction is not enough
//      for the patch test; the assumed operator must also satisfy the
//      equilibrium-consistency condition  int(B_assumed) = int(B_compatible)
//      over each element.  This is enforced by adding the constant
//          Delta = (1/V) * integral( B_compatible - B_tied ) dV
//      to the tied operator.  Delta vanishes for constant strain states (so
//      (1) is preserved) and removes only the element-average mismatch (so
//      the locking cure of the tying is preserved).
//
// Tying layout (centre-frame Cartesian strains):
//   eps11, gam13 : linear in r at +/-1/sqrt(3), quadratic in s     (2x3 points)
//   eps22, gam23 : quadratic in r, linear in s                     (3x2 points)
//   gam12        : bilinear at r,s = +/-1/sqrt(3)                  (2x2 points)
//==============================================================================

namespace shell { static bool use_tying = true; }   // false -> displacement-based

// 1D linear Lagrange weight through the two points +/- 1/sqrt(3).
inline double ty_lin ( const double x, const int idx )
{
    const double a   = 1.0 / SQRT3;
    const double sgn = (idx == 0 ? -1.0 : 1.0);
    return 0.5 * (1.0 + sgn * x / a);
}
inline double ty_lin_pt ( const int idx )
{
    const double a = 1.0 / SQRT3;
    return (idx == 0 ? -a : a);
}

// 1D quadratic Lagrange weight through the three points -sqrt(0.6), 0, +sqrt(0.6).
inline double ty_quad ( const double x, const int idx )
{
    const double b  = std::sqrt(0.6);
    const double xb = x / b;
    if (idx == 0) return 0.5 * xb * (xb - 1.0);
    if (idx == 1) return 1.0 - xb * xb;
    return 0.5 * xb * (xb + 1.0);
}
inline double ty_quad_pt ( const int idx )
{
    const double b = std::sqrt(0.6);
    return (idx == 0 ? -b : (idx == 1 ? 0.0 : b));
}

// Lamina Cartesian frame (e1,e2,e3) from the covariant base vectors:
// e3 normal to the mid-surface, e1 the in-plane projection of g_r, e2 = e3 x e1.
inline void shell_lamina_frame ( const double gr[3], const double gs[3],
                                 double e1[3], double e2[3], double e3[3] )
{
    v3cross(gr, gs, e3);  v3normalize(e3);
    const double d = v3dot(gr, e3);
    for (int k = 0; k < 3; k++) e1[k] = gr[k] - d * e3[k];
    v3normalize(e1);
    v3cross(e3, e1, e2);  v3normalize(e2);
}

// 5x5 transform T mapping tensor covariant strains
//   [ e_rr, e_ss, e_rs, e_rt, e_st ]
// to engineering Cartesian strains [ eps11, eps22, gam12, gam13, gam23 ]
// expressed in the supplied orthonormal frame (e1,e2,e3).
inline void shell_strain_transform ( const double gr[3], const double gs[3],
                                     const double gt[3],
                                     const double e1[3], const double e2[3],
                                     const double e3[3], double T[5][5] )
{
    // Contravariant base vectors = rows of inv([gr gs gt]).
    double G[3][3], Ginv[3][3];
    for (int k = 0; k < 3; k++)
    { G[k][0] = gr[k];  G[k][1] = gs[k];  G[k][2] = gt[k]; }
    m3inv(G, Ginv);

    double gR[3], gS[3], gT[3];
    for (int k = 0; k < 3; k++)
    { gR[k] = Ginv[0][k];  gS[k] = Ginv[1][k];  gT[k] = Ginv[2][k]; }

    // q[a][i] = e_a . g^i  with a in {e1,e2,e3}, i in {gR,gS,gT}.
    double q[3][3];
    const double *ea[3] = { e1, e2, e3 };
    const double *gi[3] = { gR, gS, gT };
    for (int a = 0; a < 3; a++)
        for (int i = 0; i < 3; i++)
            q[a][i] = v3dot(ea[a], gi[i]);

    // eps_ab = sum_ij q[a][i] q[b][j] e_ij , e symmetric, e_tt dropped.
    for (int row = 0; row < 5; row++)
    {
        int a, b;  double scale;
        if      (row == 0) { a = 0; b = 0; scale = 1.0; }   // eps11
        else if (row == 1) { a = 1; b = 1; scale = 1.0; }   // eps22
        else if (row == 2) { a = 0; b = 1; scale = 2.0; }   // gam12
        else if (row == 3) { a = 0; b = 2; scale = 2.0; }   // gam13
        else               { a = 1; b = 2; scale = 2.0; }   // gam23

        T[row][0] = scale * ( q[a][0]*q[b][0] );
        T[row][1] = scale * ( q[a][1]*q[b][1] );
        T[row][2] = scale * ( q[a][0]*q[b][1] + q[a][1]*q[b][0] );
        T[row][3] = scale * ( q[a][0]*q[b][2] + q[a][2]*q[b][0] );
        T[row][4] = scale * ( q[a][1]*q[b][2] + q[a][2]*q[b][1] );
    }
}

// Displacement-based Cartesian strain-displacement rows at natural point
// (r,s,t) of element e, expressed in the fixed frame (e1,e2,e3).
inline void shell_cartesian_B ( const int e, const double r, const double s,
                                const double t,
                                const double e1[3], const double e2[3],
                                const double e3[3], double B[5][45] )
{
    double cb[5][45], gr[3], gs[3], gt[3];
    shell_covariant_B(e, r, s, t, cb, gr, gs, gt);

    double T[5][5];
    shell_strain_transform(gr, gs, gt, e1, e2, e3, T);

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 45; j++)
        {
            double acc = 0.0;
            for (int k = 0; k < 5; k++) acc += T[i][k] * cb[k][j];
            B[i][j] = acc;
        }
}

// 45x45 element stiffness matrix of element e.
inline void shell_element_stiffness ( const int e, double Ke[45][45] )
{
    double D[5][5];
    shell_constitutive(D);

    for (int p = 0; p < 45; p++)
        for (int q = 0; q < 45; q++)
            Ke[p][q] = 0.0;

    // Fixed element-centre Cartesian frame used for all tying.
    double ec1[3], ec2[3], ec3[3];
    {
        double gr0[3], gs0[3], gt0[3];
        shell_base_vectors(e, 0.0, 0.0, 0.0, gr0, gs0, gt0);
        shell_lamina_frame(gr0, gs0, ec1, ec2, ec3);
    }

    int ng, nt;
    double *gx = nullptr, *gw = nullptr, *tx = nullptr, *tw = nullptr;
    get_gauss_quadrature_points(5, ng, gx, gw);    // 3-point in-plane
    get_gauss_quadrature_points(3, nt, tx, tw);    // 2-point through thickness

    const int NGP = ng * ng * nt;                  // 18 integration points
    double Btied[18][5][45];                       // tied operator per Gauss pt
    double JW[18];                                 // detJ * weight per Gauss pt

    double sumDiff[5][45];                         // integral( Bcomp - Btied )
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 45; j++)
            sumDiff[i][j] = 0.0;
    double sumJW = 0.0;

    int gp = 0;
    for (int it = 0; it < nt; it++)
    {
        const double tg = tx[it];

        // Centre-frame Cartesian strain operators at the tying points
        // (depend on tg only).  BA: eps11/gam13, BB: eps22/gam23, BC: gam12.
        double BA[6][5][45], BB[6][5][45], BC[4][5][45];
        if (shell::use_tying)
        {
            int a = 0;
            for (int i2 = 0; i2 < 2; i2++)
            for (int j2 = 0; j2 < 3; j2++, a++)
                shell_cartesian_B(e, ty_lin_pt(i2), ty_quad_pt(j2), tg,
                                  ec1, ec2, ec3, BA[a]);
            a = 0;
            for (int i2 = 0; i2 < 3; i2++)
            for (int j2 = 0; j2 < 2; j2++, a++)
                shell_cartesian_B(e, ty_quad_pt(i2), ty_lin_pt(j2), tg,
                                  ec1, ec2, ec3, BB[a]);
            a = 0;
            for (int i2 = 0; i2 < 2; i2++)
            for (int j2 = 0; j2 < 2; j2++, a++)
                shell_cartesian_B(e, ty_lin_pt(i2), ty_lin_pt(j2), tg,
                                  ec1, ec2, ec3, BC[a]);
        }

        for (int ir = 0; ir < ng; ir++)
        for (int is = 0; is < ng; is++, gp++)
        {
            const double rg = gx[ir], sg = gx[is];
            const double w  = gw[ir] * gw[is] * tw[it];

            // compatible centre-frame Cartesian strain at the Gauss point
            double Bcomp[5][45];
            shell_cartesian_B(e, rg, sg, tg, ec1, ec2, ec3, Bcomp);

            // detJ at the Gauss point
            double gr[3], gs[3], gt[3];
            shell_base_vectors(e, rg, sg, tg, gr, gs, gt);
            double G[3][3];
            for (int k = 0; k < 3; k++)
            { G[k][0] = gr[k];  G[k][1] = gs[k];  G[k][2] = gt[k]; }
            JW[gp] = m3det(G) * w;

            // assumed (tied) operator at the Gauss point
            if (shell::use_tying)
            {
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 45; j++)
                        Btied[gp][i][j] = 0.0;

                int a = 0;
                for (int i2 = 0; i2 < 2; i2++)
                for (int j2 = 0; j2 < 3; j2++, a++)
                {
                    const double wt = ty_lin(rg, i2) * ty_quad(sg, j2);
                    for (int j = 0; j < 45; j++)
                    {
                        Btied[gp][0][j] += wt * BA[a][0][j];   // eps11
                        Btied[gp][3][j] += wt * BA[a][3][j];   // gam13
                    }
                }
                a = 0;
                for (int i2 = 0; i2 < 3; i2++)
                for (int j2 = 0; j2 < 2; j2++, a++)
                {
                    const double wt = ty_quad(rg, i2) * ty_lin(sg, j2);
                    for (int j = 0; j < 45; j++)
                    {
                        Btied[gp][1][j] += wt * BB[a][1][j];   // eps22
                        Btied[gp][4][j] += wt * BB[a][4][j];   // gam23
                    }
                }
                a = 0;
                for (int i2 = 0; i2 < 2; i2++)
                for (int j2 = 0; j2 < 2; j2++, a++)
                {
                    const double wt = ty_lin(rg, i2) * ty_lin(sg, j2);
                    for (int j = 0; j < 45; j++)
                        Btied[gp][2][j] += wt * BC[a][2][j];   // gam12
                }
            }
            else
            {
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 45; j++)
                        Btied[gp][i][j] = Bcomp[i][j];
            }

            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 45; j++)
                    sumDiff[i][j] += (Bcomp[i][j] - Btied[gp][i][j]) * JW[gp];
            sumJW += JW[gp];
        }
    }

    // Mean correction enforcing  integral(B) = integral(B_compatible).
    double Delta[5][45];
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 45; j++)
            Delta[i][j] = sumDiff[i][j] / sumJW;

    // Assemble Ke with the corrected assumed-strain operator.
    for (gp = 0; gp < NGP; gp++)
    {
        double Bc[5][45];
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 45; j++)
                Bc[i][j] = Btied[gp][i][j] + Delta[i][j];

        double DB[5][45];
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 45; j++)
            {
                double acc = 0.0;
                for (int k = 0; k < 5; k++) acc += D[i][k] * Bc[k][j];
                DB[i][j] = acc;
            }

        const double f = JW[gp];
        for (int p = 0; p < 45; p++)
            for (int q = 0; q < 45; q++)
            {
                double acc = 0.0;
                for (int i = 0; i < 5; i++) acc += Bc[i][p] * DB[i][q];
                Ke[p][q] += acc * f;
            }
    }

    deallocate(gx, ng);  deallocate(gw, ng);
    deallocate(tx, nt);  deallocate(tw, nt);
}

//==============================================================================
// 6.  Global load vector, boundary conditions, HYPRE assembly and solve.
//==============================================================================

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "HYPRE_utilities.h"

namespace shell
{
    static double *F          = nullptr;   // NDOF external load vector
    static int    *fixed      = nullptr;   // NDOF : 1 if DOF is prescribed
    static double *fixed_val  = nullptr;   // NDOF : prescribed value (usually 0)
    static double *U          = nullptr;   // NDOF solution
}

// Allocate the DOF-sized arrays (call after shell_allocate has set NDOF).
inline void shell_allocate_dofs ()
{
    allocate(shell::F,         shell::NDOF);
    allocate(shell::fixed,     shell::NDOF);
    allocate(shell::fixed_val, shell::NDOF);
    allocate(shell::U,         shell::NDOF);
}

inline void shell_free_dofs ()
{
    deallocate(shell::F,         shell::NDOF);
    deallocate(shell::fixed,     shell::NDOF);
    deallocate(shell::fixed_val, shell::NDOF);
    deallocate(shell::U,         shell::NDOF);
}

// --- driver-facing boundary-condition / load helpers -------------------------

inline void shell_fix ( const int node, const int comp, const double value = 0.0 )
{
    const int d = shell::dof(node, comp);
    shell::fixed[d]     = 1;
    shell::fixed_val[d] = value;
}

inline void shell_fix_all ( const int node )            // clamp: all 5 DOFs
{
    for (int c = 0; c < 5; c++) shell_fix(node, c, 0.0);
}

inline void shell_add_point_load ( const int node, const int comp,
                                   const double value )
{
    shell::F[shell::dof(node, comp)] += value;
}

// Consistent nodal load from a uniform traction q (force per mid-surface area),
// q given in global components.  Integrated over the t = 0 mid-surface.
inline void shell_add_surface_load ( const double qx, const double qy,
                                     const double qz )
{
    const double q[3] = { qx, qy, qz };

    int ng;  double *gx = nullptr, *gw = nullptr;
    get_gauss_quadrature_points(5, ng, gx, gw);

    for (int e = 0; e < shell::NE; e++)
    {
        for (int ir = 0; ir < ng; ir++)
        for (int is = 0; is < ng; is++)
        {
            const double xi[2] = { gx[ir], gx[is] };
            double N[9], dN[9][2];
            shape_quad9(xi, N, dN);

            double gr[3] = {0,0,0}, gs[3] = {0,0,0};
            for (int b = 0; b < shell::NPE; b++)
            {
                const int nb = shell::conn[e][b];
                for (int k = 0; k < 3; k++)
                {
                    gr[k] += dN[b][0] * shell::X[nb][k];
                    gs[k] += dN[b][1] * shell::X[nb][k];
                }
            }
            double nrm[3];  v3cross(gr, gs, nrm);
            const double dA = v3norm(nrm) * gw[ir] * gw[is];

            for (int b = 0; b < shell::NPE; b++)
            {
                const int nb = shell::conn[e][b];
                for (int k = 0; k < 3; k++)
                    shell::F[shell::dof(nb, k)] += N[b] * q[k] * dA;
            }
        }
    }
    deallocate(gx, ng);  deallocate(gw, ng);
}

// --- assemble and solve ------------------------------------------------------

inline void shell_solve ( int &iters, double &res, const bool verbose = true )
{
    const int N = shell::NDOF;

    HYPRE_IJMatrix    A;   HYPRE_ParCSRMatrix parA;
    HYPRE_IJVector    b;   HYPRE_ParVector    parb;
    HYPRE_IJVector    u;   HYPRE_ParVector    paru;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N-1, 0, N-1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &u);
    HYPRE_IJVectorSetObjectType(u, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(u);

    for (int d = 0; d < N; d++)
    {
        double z = 0.0;
        HYPRE_IJVectorSetValues(u, 1, &d, &z);
        HYPRE_IJVectorSetValues(b, 1, &d, &z);
    }

    // --- element assembly with static condensation of prescribed DOFs ---
    int one = 1;
    double Ke[45][45];
    int    gdof[45];

    for (int e = 0; e < shell::NE; e++)
    {
        shell_element_stiffness(e, Ke);

        for (int a = 0; a < shell::NPE; a++)
            for (int c = 0; c < 5; c++)
                gdof[5*a + c] = shell::dof(shell::conn[e][a], c);

        for (int p = 0; p < 45; p++)
        {
            const int gp = gdof[p];
            if (shell::fixed[gp]) continue;            // prescribed row: skip

            for (int q = 0; q < 45; q++)
            {
                const int gq = gdof[q];
                const double v = Ke[p][q];
                if (v == 0.0) continue;

                if (shell::fixed[gq])
                {
                    double bv = -v * shell::fixed_val[gq];
                    HYPRE_IJVectorAddToValues(b, 1, &gp, &bv);
                }
                else
                {
                    HYPRE_IJMatrixAddToValues(A, 1, &one, &gp, &gq, &v);
                }
            }
        }
    }

    // external loads on free DOFs
    for (int d = 0; d < N; d++)
        if (!shell::fixed[d] && shell::F[d] != 0.0)
            HYPRE_IJVectorAddToValues(b, 1, &d, &shell::F[d]);

    // identity rows for prescribed DOFs
    for (int d = 0; d < N; d++)
        if (shell::fixed[d])
        {
            double uval = 1.0;
            HYPRE_IJMatrixAddToValues(A, 1, &one, &d, &d, &uval);
            HYPRE_IJVectorSetValues(b, 1, &d, &shell::fixed_val[d]);
        }

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorAssemble(u);
    HYPRE_IJMatrixGetObject(A, (void**) &parA);
    HYPRE_IJVectorGetObject(b, (void**) &parb);
    HYPRE_IJVectorGetObject(u, (void**) &paru);

    // --- PCG preconditioned by BoomerAMG (as in cantilever_*_v2.cpp) ---
    HYPRE_Solver solver, prec;
    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_PCGSetMaxIter(solver, 20000);
    HYPRE_PCGSetTol(solver, 1e-12);
    HYPRE_PCGSetTwoNorm(solver, 1);
    HYPRE_PCGSetLogging(solver, 1);
    if (verbose) HYPRE_PCGSetPrintLevel(solver, 0);

    HYPRE_BoomerAMGCreate(&prec);
    HYPRE_BoomerAMGSetOldDefault(prec);
    HYPRE_BoomerAMGSetCoarsenType(prec, 6);
    HYPRE_BoomerAMGSetRelaxType(prec, 6);
    HYPRE_BoomerAMGSetNumSweeps(prec, 1);
    HYPRE_BoomerAMGSetTol(prec, 0.0);
    HYPRE_BoomerAMGSetMaxIter(prec, 1);
    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, prec);

    HYPRE_ParCSRPCGSetup(solver, parA, parb, paru);
    HYPRE_ParCSRPCGSolve(solver, parA, parb, paru);
    HYPRE_PCGGetNumIterations(solver, &iters);
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &res);

    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_BoomerAMGDestroy(prec);

    for (int d = 0; d < N; d++)
        HYPRE_IJVectorGetValues(u, 1, &d, &shell::U[d]);

    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(u);
}

// Convenience: nodal translation component from the solution vector.
inline double shell_disp ( const int node, const int comp )
{
    return shell::U[shell::dof(node, comp)];
}

//==============================================================================
// 7.  Structured quad9 mesh generator.
//
// Builds an nx-by-ny grid of MITC9 elements over a logical unit square
// (u,v) in [0,1]^2.  The driver supplies a map (u,v) -> (x,y,z) that places
// the mid-surface in 3D.  Node ids run i + j*(2*nx+1) with the lattice index
// i in [0,2*nx], j in [0,2*ny]; node (i,j) carries (u,v) = (i/2nx, j/2ny).
//==============================================================================

namespace shell
{
    static int  gnx = 0, gny = 0;         // stored grid resolution
    static bool periodic_u = false;       // true: u-direction wraps (closed)

    // Global node id of lattice point (i,j).  When periodic_u is set, the
    // i index wraps modulo 2*gnx so the u = 0 and u = 1 edges are the same
    // nodes (a closed surface of revolution).
    inline int node_id ( const int i, const int j )
    {
        if (periodic_u)
        {
            int ii = i % (2 * gnx);
            if (ii < 0) ii += 2 * gnx;
            return ii + j * (2 * gnx);
        }
        return i + j * (2 * gnx + 1);
    }
}

// Map signature: fill xyz[3] from logical coordinates (u,v) in [0,1]^2.
typedef void (*shell_map_fn) ( double u, double v, double xyz[3] );

// Generate the structured mesh and allocate every model array.
inline void shell_build_structured ( const int nx, const int ny,
                                     shell_map_fn map )
{
    shell::gnx = nx;
    shell::gny = ny;
    shell::periodic_u = false;

    const int ni = 2 * nx + 1;
    const int nj = 2 * ny + 1;
    const int nn = ni * nj;
    const int ne = nx * ny;

    shell_allocate(nn, ne);
    shell_allocate_dofs();

    // node coordinates
    for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
        {
            const double u = double(i) / double(2 * nx);
            const double v = double(j) / double(2 * ny);
            double xyz[3];
            map(u, v, xyz);
            const int n = shell::node_id(i, j);
            shell::X[n][0] = xyz[0];
            shell::X[n][1] = xyz[1];
            shell::X[n][2] = xyz[2];
        }

    // connectivity in shape_quad9 lattice order:
    //   corners (0,0)(2,0)(2,2)(0,2), edge mids (1,0)(2,1)(1,2)(0,1), centre (1,1)
    static const int off[9][2] = {
        {0,0}, {2,0}, {2,2}, {0,2}, {1,0}, {2,1}, {1,2}, {0,1}, {1,1}
    };
    for (int ey = 0; ey < ny; ey++)
        for (int ex = 0; ex < nx; ex++)
        {
            const int e  = ex + ey * nx;
            const int i0 = 2 * ex;
            const int j0 = 2 * ey;
            for (int a = 0; a < 9; a++)
                shell::conn[e][a] =
                    shell::node_id(i0 + off[a][0], j0 + off[a][1]);
        }
}

// Closed surface of revolution: same as shell_build_structured but the
// u-direction wraps, so the u = 0 and u = 1 meridians share nodes (no seam).
// The map is sampled for u in [0,1) only.
inline void shell_build_structured_periodic ( const int nx, const int ny,
                                              shell_map_fn map )
{
    shell::gnx = nx;
    shell::gny = ny;
    shell::periodic_u = true;

    const int ni = 2 * nx;                 // unique nodes around the ring
    const int nj = 2 * ny + 1;
    const int nn = ni * nj;
    const int ne = nx * ny;

    shell_allocate(nn, ne);
    shell_allocate_dofs();

    for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
        {
            const double u = double(i) / double(2 * nx);
            const double v = double(j) / double(2 * ny);
            double xyz[3];
            map(u, v, xyz);
            const int n = shell::node_id(i, j);
            shell::X[n][0] = xyz[0];
            shell::X[n][1] = xyz[1];
            shell::X[n][2] = xyz[2];
        }

    static const int off[9][2] = {
        {0,0}, {2,0}, {2,2}, {0,2}, {1,0}, {2,1}, {1,2}, {0,1}, {1,1}
    };
    for (int ey = 0; ey < ny; ey++)
        for (int ex = 0; ex < nx; ex++)
        {
            const int e  = ex + ey * nx;
            const int i0 = 2 * ex;
            const int j0 = 2 * ey;
            for (int a = 0; a < 9; a++)
                shell::conn[e][a] =
                    shell::node_id(i0 + off[a][0], j0 + off[a][1]);
        }
}

//==============================================================================
// 8.  Visualization dump.
//
// Writes the mid-surface mesh and the nodal translation field to a plain-text
// file consumed by plot_shell_results.py (matplotlib).  Format:
//     NN <n>
//     NE <n>
//     NODES
//       x y z  ux uy uz          (one line per node)
//     ELEMS
//       n0 n1 ... n8             (one line per element, quad9 lattice order)
//==============================================================================

#include <cstdio>

inline void shell_write_viz ( const char *path )
{
    std::FILE *fp = std::fopen(path, "w");
    if (fp == nullptr)
        print_error_message("shell_write_viz: cannot open output file.");

    std::fprintf(fp, "NN %d\n", shell::NN);
    std::fprintf(fp, "NE %d\n", shell::NE);

    std::fprintf(fp, "NODES\n");
    for (int n = 0; n < shell::NN; n++)
        std::fprintf(fp, "%.10e %.10e %.10e %.10e %.10e %.10e\n",
                     shell::X[n][0], shell::X[n][1], shell::X[n][2],
                     shell::U[shell::dof(n, 0)],
                     shell::U[shell::dof(n, 1)],
                     shell::U[shell::dof(n, 2)]);

    std::fprintf(fp, "ELEMS\n");
    for (int e = 0; e < shell::NE; e++)
    {
        for (int a = 0; a < shell::NPE; a++)
            std::fprintf(fp, "%d ", shell::conn[e][a]);
        std::fprintf(fp, "\n");
    }

    std::fclose(fp);
}

#endif



