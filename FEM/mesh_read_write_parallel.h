#ifndef MESH_READ_WRITE_H
#define MESH_READ_WRITE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

#include "general.h"

//______________________________________________________________________________
// Map (d, pe) to the Gmsh MSH 2.2 element-type number for the volume element.
//
//   d  = 2 (2D) or 3 (3D)
//   pe = nodes per volume element
//
// Returns -1 if (d, pe) is not in the supported P1/P2/P3 family.
//______________________________________________________________________________
inline int gmsh_volume_type ( const int d,
                              const int pe )
{
    if (d == 2)
    {
        if (pe == 3 ) return 2;     // Tri3   (P1)
        if (pe == 6 ) return 9;     // Tri6   (P2)
        if (pe == 10) return 21;    // Tri10  (P3)
        if (pe == 4 ) return 3;     // Quad4  (Q1)
        if (pe == 9 ) return 10;    // Quad9  (Q2 full Lagrange)
        if (pe == 16) return 36;    // Quad16 (Q3 full Lagrange)
    }
    else if (d == 3)
    {
        if (pe == 4 ) return 4;     // Tet4   (P1)
        if (pe == 10) return 11;    // Tet10  (P2)
        if (pe == 20) return 29;    // Tet20  (P3)
        if (pe == 8 ) return 5;     // Hex8   (Q1)
        if (pe == 27) return 12;    // Hex27  (Q2 full Lagrange)
        if (pe == 64) return 92;    // Hex64  (Q3 full Lagrange)
    }
    return -1;
}

//______________________________________________________________________________
// Map (d, pb) to the Gmsh MSH 2.2 element-type number for the (d-1)-dimensional
// boundary element. For 2D the boundary is a Line; for 3D it's a Tri or Quad,
// uniquely determined by pb because no two boundary families share a pb in the
// P1/P2/P3 range we support.
//
// Returns -1 if (d, pb) is not supported.
//______________________________________________________________________________
inline int gmsh_boundary_type ( const int d,
                                const int pb )
{
    if (d == 2)
    {
        if (pb == 2) return 1;      // Line2 (P1 edge)
        if (pb == 3) return 8;      // Line3 (P2 edge)
        if (pb == 4) return 26;     // Line4 (P3 edge)
    }
    else if (d == 3)
    {
        if (pb == 3 ) return 2;     // Tri3   face
        if (pb == 6 ) return 9;     // Tri6   face
        if (pb == 10) return 21;    // Tri10  face
        if (pb == 4 ) return 3;     // Quad4  face
        if (pb == 9 ) return 10;    // Quad9  face
        if (pb == 16) return 36;    // Quad16 face
    }
    return -1;
}

//______________________________________________________________________________
// Number of nodes per Gmsh element type (only the types we consume; everything
// else returns 0 and the caller treats the element as "skip but still advance
// the token stream by its node list"). The caller uses this to know how many
// node ids to read off an $Elements line before deciding to keep or skip it.
//______________________________________________________________________________
inline int gmsh_nodes_per_type ( const int type )
{
    switch (type)
    {
        case 1 : return 2;     // Line2
        case 2 : return 3;     // Tri3
        case 3 : return 4;     // Quad4
        case 4 : return 4;     // Tet4
        case 5 : return 8;     // Hex8
        case 8 : return 3;     // Line3
        case 9 : return 6;     // Tri6
        case 10: return 9;     // Quad9
        case 11: return 10;    // Tet10
        case 12: return 27;    // Hex27
        case 15: return 1;     // Point
        case 21: return 10;    // Tri10
        case 26: return 4;     // Line4
        case 29: return 20;    // Tet20
        case 36: return 16;    // Quad16
        case 92: return 64;    // Hex64
        default: return 0;
    }
}

//______________________________________________________________________________
// Priority of a Dirichlet tag: lower index in `dirichlet_tags` -> higher
// priority. Returns -1 if `tag` is not in the list.
//
// Convention used by read_gmsh_mesh: when a node is hit by multiple Dirichlet
// boundary elements, the tag with the highest priority wins. The caller orders
// `dirichlet_tags` so that the most restrictive tag comes first (typically
// the "all components fixed" tag before any partial-component tag).
//______________________________________________________________________________
inline int dirichlet_priority ( const int  tag,
                                const int *dirichlet_tags,
                                const int  n_dirichlet_tags )
{
    for (int i = 0; i < n_dirichlet_tags; i++)
        if (dirichlet_tags[i] == tag)
            return n_dirichlet_tags - i;
    return -1;
}

//______________________________________________________________________________
// Read a Gmsh MSH 2.2 ASCII mesh and split its boundary entities into a
// Dirichlet node tag array and a Neumann boundary element array.
//
// PARAMETERS
//   path                [in ]  Path to the input .msh file.
//   d                   [in ]  Spatial dimension (2 or 3).
//   pe                  [in ]  Number of nodes per volume element.
//   pb                  [in ]  Number of nodes per boundary element.
//   dirichlet_tags      [in ]  Array of physical tag IDs to treat as Dirichlet.
//                              Earlier entries take priority on a shared node.
//   n_dirichlet_tags    [in ]  Length of `dirichlet_tags`.
//   neumann_tags        [in ]  Array of physical tag IDs to treat as Neumann.
//   n_neumann_tags      [in ]  Length of `neumann_tags`.
//   Nv                  [out]  Number of vertices.
//   Ne                  [out]  Number of volume elements of the (d, pe) type.
//   Nb                  [out]  Number of Neumann boundary elements found.
//   Vertices            [out]  Allocated Nv x d coordinate array.
//   Elements            [out]  Allocated Ne x pe node-id array (0-based,
//                              Gmsh native ordering).
//   BoundaryElements    [out]  Allocated Nb x pb node-id array (0-based,
//                              Gmsh native ordering) for Neumann elements only.
//   VertexTags          [out]  Allocated array of length Nv. Each entry is the
//                              winning Dirichlet tag for that node, or -1 if
//                              the node is not on any Dirichlet boundary.
//   BoundaryTags        [out]  Allocated array of length Nb. Neumann tag of
//                              each entry of `BoundaryElements`.
//
// All five output arrays are allocated by this function with general.h's
// allocate(); the caller releases them via deallocate(). Allocations with
// zero size are skipped (the corresponding pointer is left null).
//______________________________________________________________________________
inline void read_gmsh_mesh ( const char  *path,
                             const int    d,
                             const int    pe,
                             const int    pb,
                             const int   *dirichlet_tags,
                             const int    n_dirichlet_tags,
                             const int   *neumann_tags,
                             const int    n_neumann_tags,
                             int         &Nv,
                             int         &Ne,
                             int         &Nb,
                             double    **&Vertices,
                             int       **&Elements,
                             int       **&BoundaryElements,
                             int        *&VertexTags,
                             int        *&BoundaryTags )
{
    const int vol_type = gmsh_volume_type(d, pe);
    const int bdy_type = gmsh_boundary_type(d, pb);

    if (vol_type < 0)
        print_error_message("read_gmsh_mesh: unsupported (d, pe).");

    if (bdy_type < 0)
        print_error_message("read_gmsh_mesh: unsupported (d, pb).");

    //----------------------------------------------------------------
    // Pass 1: count Nv, Ne, Nb (so we can size the output arrays).
    //----------------------------------------------------------------
    std::ifstream fin1(path);

    if (!fin1)
        print_error_message("read_gmsh_mesh: cannot open input file.");

    Nv = 0;
    Ne = 0;
    Nb = 0;

    char tok[64];

    while (fin1 >> tok)
    {
        if (std::strcmp(tok, "$Nodes") == 0)
        {
            int n;
            fin1 >> n;
            Nv = n;

            for (int k = 0; k < n; k++)
            {
                int id;
                double x, y, z;
                fin1 >> id >> x >> y >> z;
            }
        }
        else if (std::strcmp(tok, "$Elements") == 0)
        {
            int total;
            fin1 >> total;

            for (int k = 0; k < total; k++)
            {
                int id, type, ntags;
                fin1 >> id >> type >> ntags;

                int phys = 0;

                for (int t = 0; t < ntags; t++)
                {
                    int v;
                    fin1 >> v;

                    if (t == 0)
                        phys = v;
                }

                const int npe = gmsh_nodes_per_type(type);

                for (int t = 0; t < npe; t++)
                {
                    int v;
                    fin1 >> v;
                }

                if (type == vol_type)
                {
                    Ne++;
                }
                else if (type == bdy_type)
                {
                    for (int i = 0; i < n_neumann_tags; i++)
                    {
                        if (neumann_tags[i] == phys)
                        {
                            Nb++;
                            break;
                        }
                    }
                }
            }
        }
    }

    fin1.close();

    //----------------------------------------------------------------
    // Allocate output arrays.
    //----------------------------------------------------------------
    Vertices         = nullptr;
    Elements         = nullptr;
    BoundaryElements = nullptr;
    VertexTags       = nullptr;
    BoundaryTags     = nullptr;

    if (Nv > 0)
    {
        allocate(Vertices,   Nv, d);
        allocate(VertexTags, Nv);

        for (int i = 0; i < Nv; i++)
            VertexTags[i] = -1;
    }

    if (Ne > 0)
        allocate(Elements, Ne, pe);

    if (Nb > 0)
    {
        allocate(BoundaryElements, Nb, pb);
        allocate(BoundaryTags,     Nb);
    }

    //----------------------------------------------------------------
    // Pass 2: fill Vertices, Elements, BoundaryElements, VertexTags,
    // BoundaryTags. Walk Dirichlet-tagged boundary elements to mark
    // VertexTags with priority resolution.
    //----------------------------------------------------------------
    std::ifstream fin2(path);

    if (!fin2)
        print_error_message("read_gmsh_mesh: cannot reopen input file.");

    int ie = 0;
    int ib = 0;

    int nodes_buf[256];   // any Gmsh element we read has <= 64 nodes

    while (fin2 >> tok)
    {
        if (std::strcmp(tok, "$Nodes") == 0)
        {
            int n;
            fin2 >> n;

            for (int k = 0; k < n; k++)
            {
                int id;
                double x, y, z;
                fin2 >> id >> x >> y >> z;

                if (d >= 1) Vertices[id-1][0] = x;
                if (d >= 2) Vertices[id-1][1] = y;
                if (d >= 3) Vertices[id-1][2] = z;
            }
        }
        else if (std::strcmp(tok, "$Elements") == 0)
        {
            int total;
            fin2 >> total;

            for (int k = 0; k < total; k++)
            {
                int id, type, ntags;
                fin2 >> id >> type >> ntags;

                int phys = 0;

                for (int t = 0; t < ntags; t++)
                {
                    int v;
                    fin2 >> v;

                    if (t == 0)
                        phys = v;
                }

                const int npe = gmsh_nodes_per_type(type);

                for (int t = 0; t < npe; t++)
                {
                    int v;
                    fin2 >> v;
                    nodes_buf[t] = v - 1;
                }

                if (type == vol_type)
                {
                    for (int t = 0; t < pe; t++)
                        Elements[ie][t] = nodes_buf[t];

                    ie++;
                }
                else if (type == bdy_type)
                {
                    int is_neumann = 0;

                    for (int i = 0; i < n_neumann_tags; i++)
                    {
                        if (neumann_tags[i] == phys)
                        {
                            is_neumann = 1;
                            break;
                        }
                    }

                    if (is_neumann)
                    {
                        for (int t = 0; t < pb; t++)
                            BoundaryElements[ib][t] = nodes_buf[t];

                        BoundaryTags[ib] = phys;

                        ib++;
                    }
                    else
                    {
                        const int new_prio = dirichlet_priority(phys,
                                                                dirichlet_tags,
                                                                n_dirichlet_tags);

                        if (new_prio >= 0)
                        {
                            for (int t = 0; t < pb; t++)
                            {
                                const int n       = nodes_buf[t];
                                const int cur_tag = VertexTags[n];
                                const int cur_pr  = (cur_tag == -1)
                                                        ? -1
                                                        : dirichlet_priority(cur_tag,
                                                                             dirichlet_tags,
                                                                             n_dirichlet_tags);

                                if (new_prio > cur_pr)
                                    VertexTags[n] = phys;
                            }
                        }
                    }
                }
            }
        }
    }

    fin2.close();
}

//______________________________________________________________________________
// Write a Gmsh MSH 2.2 file that copies the input mesh verbatim through
// "$EndElements" and appends one $NodeData block holding a 3-component vector
// field per node. In 2D the third component is padded with zero.
//
// PARAMETERS
//   in_path     [in]  Source .msh file (must be MSH 2.2 ASCII).
//   out_path    [in]  Destination .msh file.
//   view_name   [in]  Name shown by Gmsh for the data view.
//   d           [in]  Spatial dimension of the displacement field (1, 2, or 3).
//   Nv          [in]  Number of vertices in the source mesh.
//   u           [in]  Nv x d displacement array.
//______________________________________________________________________________
inline void write_gmsh_solution ( const char *in_path,
                                  const char *out_path,
                                  const char *view_name,
                                  const int   d,
                                  const int   Nv,
                                  double    **u )
{
    std::ifstream fin(in_path);

    if (!fin)
        print_error_message("write_gmsh_solution: cannot reopen input mesh.");

    std::ofstream fout(out_path);

    if (!fout)
        print_error_message("write_gmsh_solution: cannot open output file.");

    char line[4096];

    while (fin.getline(line, sizeof(line)))
    {
        fout << line << '\n';

        if (std::strcmp(line, "$EndElements") == 0)
            break;
    }

    fin.close();

    fout << "$NodeData\n";
    fout << "1\n\"" << view_name << "\"\n";
    fout << "1\n0.0\n";
    fout << "3\n0\n3\n" << Nv << '\n';
    fout << std::scientific << std::setprecision(15);

    for (int n = 0; n < Nv; n++)
    {
        fout << (n + 1);

        for (int c = 0; c < 3; c++)
        {
            const double v = (c < d) ? u[n][c] : 0.0;
            fout << ' ' << v;
        }

        fout << '\n';
    }

    fout << "$EndNodeData\n";
    fout.close();
}
#endif
