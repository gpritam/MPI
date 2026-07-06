//==============================================================================
// To run:
//    mpic++ -O3 -I/usr/include/hdf5/serial read_mesh_p3_tet.cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lmetis -o run
//    mpirun -np 4 ./run
//==============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <string>
#include <vector>

#include "hdf5.h"

#include "general.h"
#include "metis.h"

const char *mesh_path = "cantilever_p3_tetrahedral.msh";

// Tag 12: Neumann
// Tag 13: Dirichlet for all components
// Tag 14: Dirichlet for w component only

const int dirichlet_tags[] = { 13, 14 };
const int neumann_tags[]   = { 12 };

constexpr int n_dirichlet_tags = sizeof(dirichlet_tags) / sizeof(dirichlet_tags[0]);
constexpr int n_neumann_tags   = sizeof(neumann_tags)   / sizeof(neumann_tags[0]);

static const int d  = 3;     // spatial dimension
static const int pe = 20;    // nodes per volume element (P3 tet)
static const int pb = 10;    // nodes per boundary face   (P3 tri)

const int vol_type = 29;     // Tet20
const int bdy_type = 21;     // Tri10

// A node that lies on one or more tagged boundaries owns a SingleVertexTags
// record. The two std::vector<int> members store every physical tag that the
// node carries, split by boundary condition kind:
//
//   dirichlet_tags : list of Dirichlet physical tags on this node
//   neumann_tags   : list of Neumann   physical tags on this node
//
// Both counts are given by the corresponding vector's size(). At most one of
// the two lists is populated on any given vertex: a vertex on a Dirichlet
// boundary cannot simultaneously carry a Neumann tag, and vice versa. The
// two-list layout is used simply so the same record type serves both cases.

struct SingleVertexTags
{
    std::vector<int> dirichlet_tags;
    std::vector<int> neumann_tags;
};

// Returns true if 'v' appears in arr[0..n-1].
inline bool is_tag_in_list ( const int *arr, const int n, const int v )
{
    for (int i = 0; i < n; i++)
        if (arr[i] == v)
            return true;
    
    return false;
}

//______________________________________________________________________________
// HDF5 file operation mode used by write_vertices_global() and
// read_vertices_global(). One value is passed on every call:
//
//   HDF5_OPEN  : open the file (truncate if it already exists on write,
//                open read-only on read), then perform the operation.
//   HDF5_CONT  : the file is already open — just perform the operation.
//   HDF5_CLOSE : perform the operation, then close the dataset and file.
//
// The two functions share the file/dataset handles below across successive
// HDF5_CONT calls, so a whole stream of writes (or reads) reduces to a
// single open + single close pair.
//______________________________________________________________________________
enum HDF5FileOp
{
    HDF5_OPEN  = 1,
    HDF5_CONT  = 2,
    HDF5_CLOSE = 3
};

// Small file-name-keyed cache so the four helpers can juggle several HDF5
// files simultaneously across HDF5_CONT calls. Linear scan is plenty at
// this scale (typical use holds at most 4 open files at once).
struct HDF5CacheEntry
{
    char    name[256];
    hid_t   file_id;
    hid_t   dset_id;
    hsize_t extent;
};

#define HDF5_CACHE_MAX 8
static HDF5CacheEntry hdf5_cache[HDF5_CACHE_MAX];
static int            hdf5_cache_count = 0;

static HDF5CacheEntry* hdf5_cache_lookup ( const char *name )
{
    for (int i = 0; i < hdf5_cache_count; i++)
        if (std::strcmp(hdf5_cache[i].name, name) == 0)
            return &hdf5_cache[i];
    return nullptr;
}

static HDF5CacheEntry* hdf5_cache_add ( const char *name )
{
    if (hdf5_cache_count >= HDF5_CACHE_MAX)
        print_error_message("HDF5 helper cache full — increase HDF5_CACHE_MAX.");

    HDF5CacheEntry *e = &hdf5_cache[hdf5_cache_count++];
    std::strncpy(e->name, name, sizeof(e->name) - 1);
    e->name[sizeof(e->name) - 1] = '\0';
    e->file_id = -1;
    e->dset_id = -1;
    e->extent  = 0;
    return e;
}

static void hdf5_cache_evict ( HDF5CacheEntry *e )
{
    const int idx = static_cast<int>(e - hdf5_cache);
    hdf5_cache[idx] = hdf5_cache[--hdf5_cache_count];
}

//______________________________________________________________________________
// Write the coordinate of a single node into an HDF5 file at row 'node_number'.
// 'op' selects one of the three modes documented above.
//______________________________________________________________________________
inline void write_vertices_global ( const char *file_name,
                                    const int    op,
                                    const int    node_number,
                                    const double x,
                                    const double y,
                                    const double z )
{
    HDF5CacheEntry *c = hdf5_cache_lookup(file_name);
    if (c == nullptr) c = hdf5_cache_add(file_name);

    if (op == HDF5_OPEN)
    {
        c->file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t init_dims[2]  = { 0, 3 };
        hsize_t max_dims[2]   = { H5S_UNLIMITED, 3 };
        hsize_t chunk_dims[2] = { 1024, 3 };

        hid_t dspace   = H5Screate_simple(2, init_dims, max_dims);
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk_dims);

        c->dset_id = H5Dcreate2(c->file_id, "/vertices", H5T_NATIVE_DOUBLE,
                                dspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);

        H5Pclose(plist_id);
        H5Sclose(dspace);
        c->extent = 0;
    }

    if (static_cast<hsize_t>(node_number + 1) > c->extent)
    {
        hsize_t new_dims[2] = { static_cast<hsize_t>(node_number + 1), 3 };
        H5Dset_extent(c->dset_id, new_dims);
        c->extent = node_number + 1;
    }

    hid_t   file_space = H5Dget_space(c->dset_id);
    hsize_t offset[2]  = { static_cast<hsize_t>(node_number), 0 };
    hsize_t count[2]   = { 1, 3 };
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    hsize_t mem_dims[2] = { 1, 3 };
    hid_t   mem_space   = H5Screate_simple(2, mem_dims, nullptr);

    const double buf[3] = { x, y, z };
    H5Dwrite(c->dset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, buf);

    H5Sclose(mem_space);
    H5Sclose(file_space);

    if (op == HDF5_CLOSE)
    {
        H5Dclose(c->dset_id);
        H5Fclose(c->file_id);
        hdf5_cache_evict(c);
    }
}

//______________________________________________________________________________
// Read one node coordinate from row 'node_number' of the HDF5 file created by
// write_vertices_global(). 'op' selects one of the three modes documented at
// the top of this section. x, y, z are returned by reference.
//______________________________________________________________________________
inline void read_vertices_global ( const char *file_name,
                                   const int    op,
                                   const int    node_number,
                                   double      &x,
                                   double      &y,
                                   double      &z )
{
    HDF5CacheEntry *c = hdf5_cache_lookup(file_name);
    if (c == nullptr) c = hdf5_cache_add(file_name);

    if (op == HDF5_OPEN)
    {
        c->file_id = H5Fopen (file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
        c->dset_id = H5Dopen2(c->file_id, "/vertices", H5P_DEFAULT);
    }

    hid_t   file_space = H5Dget_space(c->dset_id);
    hsize_t offset[2]  = { static_cast<hsize_t>(node_number), 0 };
    hsize_t count[2]   = { 1, 3 };
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    hsize_t mem_dims[2] = { 1, 3 };
    hid_t   mem_space   = H5Screate_simple(2, mem_dims, nullptr);

    double buf[3];
    H5Dread(c->dset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, buf);

    x = buf[0];
    y = buf[1];
    z = buf[2];

    H5Sclose(mem_space);
    H5Sclose(file_space);

    if (op == HDF5_CLOSE)
    {
        H5Dclose(c->dset_id);
        H5Fclose(c->file_id);
        hdf5_cache_evict(c);
    }
}

//______________________________________________________________________________
// Write a single integer 'i' at row 'n' of an HDF5 file. 'op' selects one of
// the three modes documented alongside HDF5FileOp. Multiple files may be open
// simultaneously — each file is tracked by its own entry in the cache.
//______________________________________________________________________________
inline void write_vertices_position_global ( const char *file_name,
                                             const int    op,
                                             const int    n,
                                             const int    i )
{
    HDF5CacheEntry *c = hdf5_cache_lookup(file_name);
    if (c == nullptr) c = hdf5_cache_add(file_name);

    if (op == HDF5_OPEN)
    {
        c->file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t init_dims[1]  = { 0 };
        hsize_t max_dims[1]   = { H5S_UNLIMITED };
        hsize_t chunk_dims[1] = { 1024 };

        hid_t dspace   = H5Screate_simple(1, init_dims, max_dims);
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_dims);

        c->dset_id = H5Dcreate2(c->file_id, "/positions", H5T_NATIVE_INT,
                                dspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);

        H5Pclose(plist_id);
        H5Sclose(dspace);
        c->extent = 0;
    }

    if (static_cast<hsize_t>(n + 1) > c->extent)
    {
        hsize_t new_dims[1] = { static_cast<hsize_t>(n + 1) };
        H5Dset_extent(c->dset_id, new_dims);
        c->extent = n + 1;
    }

    hid_t   file_space = H5Dget_space(c->dset_id);
    hsize_t offset[1]  = { static_cast<hsize_t>(n) };
    hsize_t count[1]   = { 1 };
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    hsize_t mem_dims[1] = { 1 };
    hid_t   mem_space   = H5Screate_simple(1, mem_dims, nullptr);

    const int buf = i;
    H5Dwrite(c->dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &buf);

    H5Sclose(mem_space);
    H5Sclose(file_space);

    if (op == HDF5_CLOSE)
    {
        H5Dclose(c->dset_id);
        H5Fclose(c->file_id);
        hdf5_cache_evict(c);
    }
}

//______________________________________________________________________________
// Read the integer stored at row 'i' of an HDF5 file written by
// write_vertices_position_global(). 'op' selects one of the three modes.
// The result is returned in 'vertex_position' (by reference).
//______________________________________________________________________________
inline void read_vertices_position_global ( const char *file_name,
                                            const int    op,
                                            const int    i,
                                            int         &vertex_position )
{
    HDF5CacheEntry *c = hdf5_cache_lookup(file_name);
    if (c == nullptr) c = hdf5_cache_add(file_name);

    if (op == HDF5_OPEN)
    {
        c->file_id = H5Fopen (file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
        c->dset_id = H5Dopen2(c->file_id, "/positions", H5P_DEFAULT);
    }

    hid_t   file_space = H5Dget_space(c->dset_id);
    hsize_t offset[1]  = { static_cast<hsize_t>(i) };
    hsize_t count[1]   = { 1 };
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    hsize_t mem_dims[1] = { 1 };
    hid_t   mem_space   = H5Screate_simple(1, mem_dims, nullptr);

    int buf;
    H5Dread(c->dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &buf);
    vertex_position = buf;

    H5Sclose(mem_space);
    H5Sclose(file_space);

    if (op == HDF5_CLOSE)
    {
        H5Dclose(c->dset_id);
        H5Fclose(c->file_id);
        hdf5_cache_evict(c);
    }
}

//______________________________________________________________________________
// 
//______________________________________________________________________________
int main ( int argc, char *argv[] )
{
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int Nv_global = 0;
    int Ne_global = 0;
    int Nb_global = 0;

    int **Elements_global         = nullptr;
    int **BoundaryElements_global = nullptr;
    int  *BoundaryTags_global     = nullptr;
    int  *rank_of_vertices        = nullptr;
    int  *rank_of_elements        = nullptr;

    SingleVertexTags **VertexTags_global = nullptr;

    int Nv = 0;
    int Nv_private = 0;
    int Ne = 0;
    int Ne_private = 0;
    int Nb = 0;

    double **Vertices = nullptr;

    int **Elements         = nullptr;
    int **BoundaryElements = nullptr;
    int  *BoundaryTags     = nullptr;

    SingleVertexTags **VertexTags = nullptr;

    int *GlobalVertexNumbers = nullptr;

    int Nsender_ranks;
    int Nreceive_buffer;
    int *Sender_ranks;
    int *Receive_counts;
    int *Cumulative_receive_counts;
    int *Receive_buffer_local_node_numbers;

    int Nreceiver_ranks;
    int Nsend_buffer;
    int *Receiver_ranks;
    int *Send_counts;
    int *Cumulative_send_counts;
    int *Send_buffer_local_node_numbers;

    if (rank == 0)
    {
        char gmsh_string[64];
        int  gmsh_int[256];

        // Pass 1: count Nv_global, Ne_global, Nb_global
        std::ifstream gmsh_read_first_pass(mesh_path);

        if (!gmsh_read_first_pass)
            print_error_message("read_gmsh_mesh: cannot open input file.");

        while (gmsh_read_first_pass >> gmsh_string)
        {
            if (std::strcmp(gmsh_string, "$Nodes") == 0)
            {
                gmsh_read_first_pass >> Nv_global;

                std::string node_line;
                std::getline(gmsh_read_first_pass, node_line);   // discard rest of the count line

                for (int k = 0; k < Nv_global; ++k)
                    std::getline(gmsh_read_first_pass, node_line);
            }
            else if (std::strcmp(gmsh_string, "$Elements") == 0)
            {
                int Ne_total;

                gmsh_read_first_pass >> Ne_total;

                for (int k = 0; k < Ne_total; k++)
                {
                    int id, type, ntags, physical, v;

                    gmsh_read_first_pass >> id >> type >> ntags;

                    for (int t = 0; t < ntags; t++)
                    {
                        if (t)
                            gmsh_read_first_pass >> v;
                        else
                            gmsh_read_first_pass >> physical;
                    }

                    // Skip the rest of the line: we don't need the node IDs in pass 1.
                    std::string rest_of_line;
                    std::getline(gmsh_read_first_pass, rest_of_line);

                    if (type == vol_type)
                        Ne_global++;
                    else if (type == bdy_type)
                        if (is_tag_in_list(neumann_tags, n_neumann_tags, physical))
                            Nb_global++;
                }
            }
        }
        
        gmsh_read_first_pass.close();

        // Allocate output arrays
        if (Nv_global > 0)
        {
            allocate(VertexTags_global, Nv_global);

            for (int i = 0; i < Nv_global; ++i)
                VertexTags_global[i] = nullptr;
        }

        if (Ne_global > 0)
            allocate(Elements_global, Ne_global, pe);
        
        if (Nb_global > 0)
        {
            allocate(BoundaryElements_global, Nb_global, pb);
            allocate(BoundaryTags_global,     Nb_global);
        }

        // Pass 2: fill arrays
        std::ifstream gmsh_read_second_pass(mesh_path);

        if (!gmsh_read_second_pass)
            print_error_message("read_gmsh_mesh: cannot reopen input file.");

        int  ie = 0;
        int  ib = 0;

        while (gmsh_read_second_pass >> gmsh_string)
        {
            if (std::strcmp(gmsh_string, "$Nodes") == 0)
            {
                int n;
                gmsh_read_second_pass >> n;

                // HDF5_OPEN on the first call truncates any pre-existing file
                // for us, so no explicit std::remove() is needed.
                for (int k = 0; k < n; k++)
                {
                    int    id;
                    double x, y, z;

                    gmsh_read_second_pass >> id >> x >> y >> z;

                    const int op = (k == 0)     ? HDF5_OPEN
                                 : (k == n - 1) ? HDF5_CLOSE
                                                : HDF5_CONT;

                    write_vertices_global("vertices_global.h5", op, id - 1, x, y, z);
                }
            }
            else if (std::strcmp(gmsh_string, "$Elements") == 0)
            {
                int Ne_total;
                
                gmsh_read_second_pass >> Ne_total;

                for (int k = 0; k < Ne_total; k++)
                {
                    int id, type, ntags, physical, v;

                    gmsh_read_second_pass >> id >> type >> ntags;

                    for (int t = 0; t < ntags; t++)
                    {
                        if (t)
                            gmsh_read_second_pass >> v;
                        else
                            gmsh_read_second_pass >> physical;
                    }

                    const int npe = (type == vol_type) ? pe : (type == bdy_type) ? pb : 0;

                    if (npe)
                    {
                        for (int t = 0; t < npe; t++)
                        {
                            gmsh_read_second_pass >> v;
                            
                            gmsh_int[t] = v - 1;
                        }
                    }
                    else{
                        // Skip the rest of the line: we don't need the node IDs.
                        std::string rest_of_line;
                        std::getline(gmsh_read_second_pass, rest_of_line);
                    }

                    if (type == vol_type)
                    {
                        for (int t = 0; t < pe; t++)
                            Elements_global[ie][t] = gmsh_int[t];

                        ie++;
                    }
                    else if (type == bdy_type)
                    {
                        if (is_tag_in_list(neumann_tags, n_neumann_tags, physical))
                        {
                            for (int t = 0; t < pb; t++)
                                BoundaryElements_global[ib][t] = gmsh_int[t];

                            BoundaryTags_global[ib] = physical;
                            ib++;
                        }
                        else
                        {
                            if (is_tag_in_list(dirichlet_tags, n_dirichlet_tags, physical))
                            {
                                for (int t = 0; t < pb; t++)
                                {
                                    const int n = gmsh_int[t];

                                    // Allocate the per-node record on first hit.
                                    // The empty vector inside is default-constructed.
                                    if (VertexTags_global[n] == nullptr)
                                        VertexTags_global[n] = new SingleVertexTags;

                                    // Append 'physical' if not already in the list.
                                    bool already = false;

                                    for (int j = 0; j < static_cast<int>(VertexTags_global[n]->dirichlet_tags.size()); j++)
                                    {
                                        if (VertexTags_global[n]->dirichlet_tags[j] == physical)
                                        {
                                            already = true;
                                            break;
                                        }
                                    }

                                    if (!already)
                                        VertexTags_global[n]->dirichlet_tags.push_back(physical);
                                }
                            }
                        }
                    }
                }
            }
        }

        gmsh_read_second_pass.close();

        //--------------------------------------------------------------------------
        // Report
        //--------------------------------------------------------------------------
        int Ndn_global = 0;

        for (int n = 0; n < Nv_global; n++)
        {
            if (VertexTags_global[n] == nullptr)
                continue;
            
            Ndn_global++;
        }

        std::cout << "Read mesh: '" << mesh_path << "'" << std::endl;
        std::cout << "  Nv_global             = " << Nv_global   << std::endl;
        std::cout << "  Ne_global             = " << Ne_global   << std::endl;
        std::cout << "  Nb_global             = " << Nb_global   << std::endl;
        std::cout << "  Ndn_global            = " << Ndn_global   << std::endl;

        int Nv_local[size];
        int Ne_local[size];

        allocate(rank_of_vertices, Nv_global);
        allocate(rank_of_elements, Ne_global);

        // Partition the volume mesh into 'size' parts with METIS_PartMeshNodal.
        // Build CSR-style (eptr, eind) connectivity from Elements_global; METIS
        // returns both an element partition (epart) and a node partition (npart).
        {
            idx_t ne_metis = Ne_global;
            idx_t nn_metis = Nv_global;
            idx_t nparts   = size;

            idx_t *eptr = nullptr;
            idx_t *eind = nullptr;
            allocate(eptr, Ne_global + 1);
            allocate(eind, Ne_global * pe);

            for (int e = 0; e < Ne_global; e++)
            {
                eptr[e] = e * pe;
                for (int j = 0; j < pe; j++)
                    eind[e * pe + j] = Elements_global[e][j];
            }
            eptr[Ne_global] = Ne_global * pe;

            idx_t *epart_metis = nullptr;
            idx_t *npart_metis = nullptr;
            allocate(epart_metis, Ne_global);
            allocate(npart_metis, Nv_global);

            idx_t objval = 0;
            METIS_PartMeshNodal(&ne_metis, &nn_metis,
                                eptr, eind,
                                nullptr, nullptr,   // vwgt, vsize
                                &nparts,
                                nullptr,            // tpwgts
                                nullptr,            // options
                                &objval,
                                epart_metis, npart_metis);

            for (int i = 0; i < Nv_global; i++)
                rank_of_vertices[i] = static_cast<int>(npart_metis[i]);

            // Take METIS's jointly-optimised element partition as-is.
            for (int i = 0; i < Ne_global; i++)
                rank_of_elements[i] = static_cast<int>(epart_metis[i]);

            // Get Ne_local[size]
            for (int r = 0; r < size; r++)
            {
                int cnt_elements = 0;
                for (int i = 0; i < Ne_global; i++)
                    if (rank_of_elements[i] == r)
                        cnt_elements++;
                
                Ne_local[r] = cnt_elements;
            }

            deallocate(eptr,        Ne_global + 1);
            deallocate(eind,        Ne_global * pe);
            deallocate(epart_metis, Ne_global);
            deallocate(npart_metis, Nv_global);
        }

        // Renumber global nodes such that nodes are contiguous on ranks.
        // As nodes are changing position in the 'Vertices' array, accordinly modify the following:
        // 1. VertexTags_global
        // 2. rank_of_vertices
        // 3. Elements_global
        // 4. BoundaryElements_global
        {
            // New vertex positions
            int n = 0;
            double x, y, z;

            for (int r = 0; r < size; r++)
            {
                int cnt_nodes = 0;
                for (int i = 0; i < Nv_global; i++)
                {
                    if (rank_of_vertices[i] == r)
                    {
                        const int op = (n == 0) ? HDF5_OPEN
                         : (n == Nv_global - 1) ? HDF5_CLOSE
                                                : HDF5_CONT;

                        read_vertices_global("vertices_global.h5", op, i, x, y, z);
                        write_vertices_global("vertices_global_new.h5", op, n, x, y, z);

                        write_vertices_position_global("vertex_position_global.h5", op, n, i);
                        write_vertices_position_global("inverse_vertex_position_global.h5", op, i, n);

                        cnt_nodes++;
                        n++;
                    }
                }
                Nv_local[r] = cnt_nodes;
            }

            // vertices_global_new.h5 now holds the reordered coordinates.
            // Replace the original coordinate file with it so downstream
            // code can keep using the "vertices_global.h5" name.
            if (std::remove("vertices_global.h5") != 0)
                print_error_message("failed to remove vertices_global.h5");

            if (std::rename("vertices_global_new.h5", "vertices_global.h5") != 0)
                print_error_message("failed to rename vertices_global_new.h5 -> vertices_global.h5");

            // Modify rank_of_vertices
            n = 0;

            for (int r = 0; r < size; r++)
            {
                for (int i = 0; i < Nv_local[r]; ++i)
                {
                    rank_of_vertices[n] = r;
                    n++;
                }
            }

            // Modify VertexTags_global
            // Rearrange VertexTags_global so that entry i corresponds to the vertex whose
            // original global index is vertex_position_global[i]. The transformation is
            // done through a small scratch file that is deleted at the end.
            //
            // Format per line (tab-separated):
            //     nullptr  -> "-1"
            //     non-null -> "<kind>  <count>  <tag_1>  <tag_2>  ...  <tag_count>"
            //
            // where <kind> = 0 marks a Dirichlet-tagged vertex and <kind> = 1 marks a
            // Neumann-tagged vertex. Because a vertex holds either Dirichlet tags or
            // Neumann tags but never both, a single kind marker is enough for the reader 
            // to decide which vector inside SingleVertexTags the following <count> tag 
            // values should be pushed into.
            {
                const char *scratch_path = "vertex_tags_reorder.txt";

                std::ofstream fout(scratch_path);

                if (!fout)
                    print_error_message("VertexTags rearrange: cannot open scratch file for write!");

                for (int i = 0; i < Nv_global; i++)
                {
                    const int op = (i == 0) ? HDF5_OPEN
                         : (i == Nv_global - 1) ? HDF5_CLOSE
                                                : HDF5_CONT;
                    int vertex_position;
                    read_vertices_position_global("vertex_position_global.h5", op, i, vertex_position);

                    SingleVertexTags *v = VertexTags_global[vertex_position];

                    if (v == nullptr)
                    {
                        fout << -1 << std::endl;
                        continue;
                    }

                    const int nd = static_cast<int>(v->dirichlet_tags.size());
                    const int nn = static_cast<int>(v->neumann_tags.size());

                    // Exactly one of nd, nn is non-zero for a non-null vertex.
                    fout << (nd ? 0 : 1);

                    if (nd)
                    {
                        fout << "\t" << nd;

                        for (int j = 0; j < nd; j++)
                            fout << "\t" << v->dirichlet_tags[j];
                    }
                    else
                    {
                        fout << "\t" << nn;

                        for (int j = 0; j < nn; j++)
                            fout << "\t" << v->neumann_tags[j];
                    }

                    fout << std::endl;
                }

                fout.close();

                // vertex_position_global.h5 is no longer needed downstream.
                std::remove("vertex_position_global.h5");

                // free the existing SingleVertexTags records
                for (int i = 0; i < Nv_global; i++)
                {
                    if (VertexTags_global[i] != nullptr)
                    {
                        delete VertexTags_global[i];
                        VertexTags_global[i] = nullptr;
                    }
                }

                // reload in the new order
                std::ifstream fin(scratch_path);

                if (!fin)
                    print_error_message("VertexTags rearrange: cannot reopen scratch file for read!");

                for (int i = 0; i < Nv_global; i++)
                {
                    int tag_type, tag_number, tag;

                    fin >> tag_type;

                    if (tag_type == -1)
                        continue;

                    VertexTags_global[i] = new SingleVertexTags;

                    // read <tag_number> then <tag_number> tag values.
                    fin >> tag_number >> tag;

                    std::vector<int> &vec = ( tag_type == 0 ? VertexTags_global[i]->dirichlet_tags
                                                            : VertexTags_global[i]->neumann_tags);
                    vec.push_back(tag);

                    for (int j = 1; j < tag_number; j++)
                    {
                        fin >> tag;
                        vec.push_back(tag);
                    }
                }

                fin.close();

                // delete the scratch file
                std::remove(scratch_path);
            }

            int vertex_position;
            read_vertices_position_global("inverse_vertex_position_global.h5", HDF5_OPEN, 0, vertex_position);

            // Modify Elements_global
            for (int i = 0; i < Ne_global; ++i)
            {
                for (int j = 0; j < pe; ++j)
                {
                    read_vertices_position_global("inverse_vertex_position_global.h5", HDF5_CONT, Elements_global[i][j], vertex_position);
                    Elements_global[i][j] = vertex_position;
                }
            }

            // Modify BoundaryElements_global
            for (int i = 0; i < Nb_global; ++i)
            {
                for (int j = 0; j < pb; ++j)
                {
                    read_vertices_position_global("inverse_vertex_position_global.h5", HDF5_CONT, BoundaryElements_global[i][j], vertex_position);
                    BoundaryElements_global[i][j] = vertex_position;
                }
            }

            read_vertices_position_global("inverse_vertex_position_global.h5", HDF5_CLOSE, 0, vertex_position);

            // inverse_vertex_position_global.h5 is no longer needed downstream.
            std::remove("inverse_vertex_position_global.h5");
        }

        // Get rank of elements

        // rank_of_elements is now populated directly from METIS's epart output
        // right after the METIS_PartMeshNodal call above, so the "min-rank-of-
        // vertex" derivation is no longer needed.

        /*// First strategy
        
        allocate(rank_of_elements, Ne_global);
        
        for (int i = 0; i < Ne_global; ++i)
            rank_of_elements[i] = -1;

        for (int r = 0; r < size; ++r)
        {
            for (int i = 0; i < Ne_global; ++i)
            {
                if (rank_of_elements[i] == -1)
                {
                    for (int j = 0; j < pe; ++j)
                    {
                        if (rank_of_vertices[Elements_global[i][j]] == r)
                        {
                            rank_of_elements[i] = r;
                            break;
                        }
                    }
                }
            }
        }*/
        
        /*// Alternative strategy
        
        allocate(rank_of_elements, Ne_global);
        
        for (int i = 0; i < Ne_global; ++i)
        {
            int best = size;
        
            for (int j = 0; j < pe; ++j)
            {
                const int r = rank_of_vertices[Elements_global[i][j]];
        
                if (r < best)
                    best = r;
            }
        
            rank_of_elements[i] = best;
        }*/

        // Report nodes and elements per rank.
        std::cout << "  METIS partition into " << size << " parts:" << std::endl;

        for (int r = 0; r < size; r++)
        {
            int cnt_nodes = 0;
            int cnt_elements = 0;

            for (int i = 0; i < Nv_global; i++)
                if (rank_of_vertices[i] == r)
                    cnt_nodes++;
            
            for (int i = 0; i < Ne_global; i++)
                if (rank_of_elements[i] == r)
                    cnt_elements++;

            std::cout << "    rank " << r << ": " << cnt_nodes << " nodes, " << cnt_elements << " elements" << std::endl;
        }

        int Ne_extended[size];
        int Nv_extended[size];
        int Nb_extended[size];

        int *Selected_vertices          = nullptr;
        int *Selected_elements          = nullptr;
        int *Selected_boundary_elements = nullptr;

        allocate(Selected_vertices,          Nv_global);
        allocate(Selected_elements,          Ne_global);
        allocate(Selected_boundary_elements, Nb_global);

        for (int r = (size-1); r >= 0; r--)
        {
            for (int i = 0; i < Nv_global; ++i)
                Selected_vertices[i] = -1;

            for (int i = 0; i < Ne_global; ++i)
                Selected_elements[i] = -1;

            for (int i = 0; i < Nb_global; ++i)
                Selected__boundary_elements[i] = -1;

            // Get extended element list for each rank
            // Mark the extra elements with 0 that does not belong to rank r
            int cnt_elements = 0;

            for (int i = 0; i < Ne_global; ++i)
            {
                if (rank_of_elements[i] != r)
                {
                    for (int j = 0; j < pe; ++j)
                    {
                        if (rank_of_vertices[Elements_global[i][j]] == r)
                        {
                            Selected_elements[i] = 0;
                            cnt_elements++;
                            break;
                        }
                    }
                }
            }

            Ne_extended[r] = Ne_local[r] + cnt_elements;

            allocate(Elements,Ne_extended[r],pe);

            int cnt = 0;
            cnt_elements = Ne_local[r];

            for (int i = 0; i < Ne_global; ++i)
            {
                if (rank_of_elements[i] == r)
                {
                    for (int j = 0; j < pe; ++j)
                        Elements[cnt][j] = Elements_global[i][j];

                    cnt++;
                }

                if (Selected_elements[i] == 0)
                {
                    for (int j = 0; j < pe; ++j)
                        Elements[cnt_elements][j] = Elements_global[i][j];

                    cnt_elements++;
                }
            }

            // Get extended node list for each rank
            // Get global node numbers in the extended node list
            for (int i = 0; i < Ne_global; ++i)
            {
                if ( (Selected_elements[i] == 0) || (rank_of_elements[i] == r) )
                {
                    for (int j = 0; j < pe; ++j)
                    {
                        int vrtx = Elements_global[i][j];

                        if (rank_of_vertices[vrtx] != r)
                            Selected_vertices[vrtx] = 0;
                    }
                }
            }

            cnt = 0;
            int cnt_vertices = Nv_local[r];

            for (int i = 0; i < Nv_global; ++i)
            {
                if (rank_of_vertices[i] == r)
                {
                    Selected_vertices[i] = cnt;
                    cnt++;
                }

                if (Selected_vertices[i] == 0)
                {
	                Selected_vertices[i] = cnt_vertices;
                    cnt_vertices++;
                }
            }

            Nv_extended[r] = cnt_vertices;

            allocate(Vertices,Nv_extended[r],3);

            cnt = 0;
            cnt_vertices = 0;

            double x, y, z;
            read_vertices_global("vertices_global.h5", HDF5_OPEN, 0, x, y, z);

            for (int i = 0; i < Nv_global; ++i)
            {
                if (rank_of_vertices[i] == r)
                {
                    read_vertices_global("vertices_global.h5", HDF5_CONT, i, x, y, z);

                    Vertices[cnt][0] = x;
                    Vertices[cnt][1] = y;
                    Vertices[cnt][2] = z;

                    cnt++;
                }
                else
                {
                    if (Selected_vertices[i] != -1)
                    {
                        read_vertices_global("vertices_global.h5", HDF5_CONT, i, x, y, z);

                        Vertices[cnt_vertices][0] = x;
                        Vertices[cnt_vertices][1] = y;
                        Vertices[cnt_vertices][2] = z;

                        cnt_vertices++;
                    }
                }
            }

            read_vertices_global("vertices_global.h5", HDF5_CLOSE, 0, x, y, z);

            allocate(GlobalVertexNumbers,Nv_extended[r]);

            for (int i = 0; i < Nv_global; ++i)
                if (Selected_vertices[i] != -1)
                    GlobalVertexNumbers[Selected_vertices[i]] = i;

            // Get local 'VertexTags' for each rank
            if (r == 0)
            {
                if (Nv > 0)
                {
                    allocate(VertexTags, Nv);

                    for (int i = 0; i < Nv; ++i)
                        VertexTags[i] = nullptr;

                    int cnt_vertex_tags = 0;

                    for (int i = 0; i < Nv_global; ++i)
                    {
                        if (Selected_vertices[i] != -1)
                        {
                            if (VertexTags_global[i] != nullptr)
                            {
                                VertexTags[cnt_vertex_tags] = new SingleVertexTags;

                                for (int j = 0; j < static_cast<int>(VertexTags_global[i]->dirichlet_tags.size()); ++j)
                                    VertexTags[cnt_vertex_tags]->dirichlet_tags.push_back(VertexTags_global[i]->dirichlet_tags[j]);

                                for (int j = 0; j < static_cast<int>(VertexTags_global[i]->neumann_tags.size()); ++j)
                                    VertexTags[cnt_vertex_tags]->neumann_tags.push_back(VertexTags_global[i]->neumann_tags[j]);
                            }

                            cnt_vertex_tags++;
                        }
                    }
                }
            }
            else
            {
                int N_node_tags = 0;
                int *node_types;
                int *number_of_tags;
                int *tag_list;

                for (int i = 0; i < Nv_global; ++i)
                    if (Selected_vertices[i] != -1)
                        if (VertexTags_global[i] != nullptr)
                            N_node_tags += ( (static_cast<int>(VertexTags_global[i]->dirichlet_tags.size()) + (static_cast<int>(VertexTags_global[i]->neumann_tags.size()) );

                // Send 'N_node_tags' to rank r


                if (N_node_tags != 0)
                {
                    allocate(node_types,Nv);
                    allocate(number_of_tags,Nv);
                    allocate(tag_list,N_node_tags);

                    int cnt_vertex_tags = 0;
                    int cnt_tags = 0;

                    for (int i = 0; i < Nv_global; ++i)
                    {
                        if (Selected_vertices[i] != -1)
                        {
                            if (VertexTags_global[i] != nullptr)
                            {
                                bool dirichlet_node;
                                int ntags;

                                dirichlet_node = (static_cast<int>(VertexTags_global[i]->dirichlet_tags.size()) > 0 ? true : false);

                                if (dirichlet_node)
                                {
                                    ntags = (static_cast<int>(VertexTags_global[i]->dirichlet_tags.size());

                                    node_types[cnt_vertex_tags] = 1;
                                    number_of_tags[cnt_vertex_tags] = ntags;

                                    for (int j = 0; j < ntags; ++j)
                                    {
                                        tag_list[cnt_tags] = VertexTags_global[i]->dirichlet_tags[j];
                                        cnt_tags++;
                                    }
                                }
                                else
                                {
                                    ntags = (static_cast<int>(VertexTags_global[i]->neumann_tags.size());
                                    
                                    node_types[cnt_vertex_tags] = 2;
                                    number_of_tags[cnt_vertex_tags] = ntags;

                                    for (int j = 0; j < ntags; ++j)
                                    {
                                        tag_list[cnt_tags] = VertexTags_global[i]->neumann_tags[j];
                                        cnt_tags++;
                                    }
                                }
                            }
                            else
                            {
                                node_types[cnt_vertex_tags] = 0;
                                number_of_tags[cnt_vertex_tags] = 0;
                            }

                            cnt_vertex_tags++;
                        }
                    }

                    // Send 'N_node_tags' to rank r
                    // Send 'node_types' array to rank r
                    // Send 'number_of_tags' array to rank r
                    // Send 'tag_list' array to rank r





                    deallocate(node_types,Nv);
                    deallocate(number_of_tags,Nv);
                    deallocate(tag_list,N_node_tags);
                }
            }

            // Get local 'BoundaryTags' for each rank
            int cnt_boundary_elements = 0;

            for (int i = 0; i < Nb_global; ++i)
            {
                bool selected = true;

                for (int j = 0; j < pb; ++j)
                {
                    int vrtx = BoundaryElements_global[i][j];

                    if (Selected_vertices[vrtx] == -1)
                    {
                        selected = false;
                        break;
                    }
                }

                if (selected)
                {
                    Selected_boundary_elements[i] = 0;
                    cnt_boundary_elements++;
                }
            }

            Nb_extended[r] = cnt_boundary_elements;

            allocate(BoundaryElements,Nb_extended[r],pb);
            allocate(BoundaryTags,Nb_extended[r]);

            cnt_boundary_elements = 0;

            for (int i = 0; i < Nb_global; ++i)
            {
                if (Selected_boundary_elements[i] == 0)
                {
                    for (int j = 0; j < pb; ++j)
                    {
                        int vrtx = BoundaryElements_global[i][j];
                        BoundaryElements[cnt_boundary_elements][j] = Selected_vertices[vrtx];
                    }

                    BoundaryTags[cnt_boundary_elements] = BoundaryTags_global[i];
                    cnt_boundary_elements++;
                }
            }

            // Modify 'Elements' to introduce local node numbers
            for (int i = 0; i < Ne_extended[r]; ++i)
            {
                for (int j = 0; j < pe; ++j)
                {
                    int vrtx = Elements[i][j];

                    Elements[i][j] = Selected_vertices[vrtx];
                }
            }

            // Get 'Nsender_ranks' (how many processors will send data to rank r)
            // Get 'Sender_ranks[Nsender_ranks]'
            // Get 'Receive_counts[Nsender_ranks]'
            // Get 'Cumulative_receive_counts[Nsender_ranks]'
            // Get 'Nreceive_buffer' = Cumulative_receive_counts[Nsender_ranks-1]
            // Get 'Receive_buffer_local_node_numbers[Nreceive_buffer]'

            int *rank_wise_counter;

            allocate(rank_wise_counter,size);

            Nreceive_buffer = 0;

            for (int i = 0; i < size; ++i)
                rank_wise_counter[i] = 0;

            for (int i = 0; i < Nv_global; ++i)
            {
                if (Selected_vertices[i] != -1)
                {
                    int node_belongs_to = rank_of_vertices[i];

                    if (node_belongs_to != r)
                    {
                        rank_wise_counter[node_belongs_to]++;
                        Nreceive_buffer++;
                    }
                }
            }

            Nsender_ranks = 0;

            for (int i = 0; i < size; ++i)
                if (rank_wise_counter[i])
                    Nsender_ranks++;
            
            allocate(Sender_ranks,Nsender_ranks);
            allocate(Receive_counts,Nsender_ranks);
            allocate(Cumulative_receive_counts,Nsender_ranks);
            allocate(Receive_buffer_local_node_numbers,Nreceive_buffer);

            Nsender_ranks = 0;

            for (int i = 0; i < size; ++i)
            {
                if (rank_wise_counter[i])
                {
                    Sender_ranks[Nsender_ranks] = i;
                    Receive_counts[Nsender_ranks] = rank_wise_counter[i];
                    Cumulative_receive_counts[Nsender_ranks] = Receive_counts[Nsender_ranks] + (Nsender_ranks == 0 ? 0 : Cumulative_receive_counts[Nsender_ranks-1]);
                    Nsender_ranks++;
                }
            }

            for (int i = 0; i < Nsender_ranks; ++i)
                Send_counts[i] = 0;

            for (int i = 0; i < Nv_global; ++i)
            {
                if (Selected_vertices[i] != -1)
                {
                    int node_belongs_to = rank_of_vertices[i];

                    if (node_belongs_to != r)
                    {
                        for (int j = 0; j < Nsender_ranks; ++j)
                        {
                            if (node_belongs_to == Sender_ranks[j])
                            {
                                int start_index = (j == 0 ? 0 : Cumulative_receive_counts[j-1]);
                                Send_buffer_local_node_numbers[start_index+Send_counts[j]] = Selected_vertices[i]; 
                                Send_counts[j]++;
                                break;
                            }
                        }
                    }

                    if (node_belongs_to != r)
                    {
                        rank_wise_counter[node_belongs_to]++;
                        Nreceive_buffer++;
                    }
                }
            }

            deallocate(rank_wise_counter,size);

            if (r != 0)
            {
                deallocate(Elements,Ne_extended[r],pe);
                deallocate(Vertices,Nv_extended[r],3);
                deallocate(BoundaryElements,Nb_extended[r],pb);
                deallocate(BoundaryTags,Nb_extended[r]);

                deallocate(GlobalVertexNumbers,Nv_extended[r]);

                deallocate(Sender_ranks,Nsender_ranks);
                deallocate(Receive_counts,Nsender_ranks);
                deallocate(Cumulative_receive_counts,Nsender_ranks);
                deallocate(Receive_buffer_local_node_numbers,Nreceive_buffer);
            }
        }

        if (Selected_vertices          != nullptr) deallocate(Selected_vertices, Nv_global);
        if (Selected_elements          != nullptr) deallocate(Selected_elements, Ne_global);
        if (Selected_boundary_elements != nullptr) deallocate(Selected_boundary_elements, Nb_global);
    }
    else
    {
        // Receive data from rank 0



    }

    // Prepare 'Send_buffer'





    // Initialize 'Receive_buffer_i' for receiving integers
    // Initialize 'Receive_buffer_d' for receiving doubles



    // Initialize 'Send_buffer_i' for sending integers
    // Initialize 'Send_buffer_d' for sending doubles




    // Deallocation on every rank
    if (rank_of_vertices        != nullptr) deallocate(rank_of_vertices,        Nv_global);
    if (rank_of_elements        != nullptr) deallocate(rank_of_elements,        Ne_global);

    if (Elements_global         != nullptr) deallocate(Elements_global,         Ne_global, pe);
    if (BoundaryElements_global != nullptr) deallocate(BoundaryElements_global, Nb_global, pb);
    if (BoundaryTags_global     != nullptr) deallocate(BoundaryTags_global,     Nb_global);

    if (VertexTags_global != nullptr)
    {
        for (int n = 0; n < Nv_global; n++)
        {
            if (VertexTags_global[n] == nullptr)
                continue;
            
            // The std::vector inside is destroyed by 'delete', releasing its storage.
            delete VertexTags_global[n];
        }

        deallocate(VertexTags_global, Nv_global);
    }

    if (Vertices            != nullptr) deallocate(Vertices,            Nv, 3);
    if (Elements            != nullptr) deallocate(Elements,            Ne, pe);
    if (BoundaryElements    != nullptr) deallocate(BoundaryElements,    Nb, pb);
    if (BoundaryTags        != nullptr) deallocate(BoundaryTags,        Nb);
    if (GlobalVertexNumbers != nullptr) deallocate(GlobalVertexNumbers, Nv);

    if (Sender_ranks                      != nullptr) deallocate(Sender_ranks,                      Nsender_ranks);
    if (Receive_counts                    != nullptr) deallocate(Receive_counts,                    Nsender_ranks);
    if (Cumulative_receive_counts         != nullptr) deallocate(Cumulative_receive_counts,         Nsender_ranks);
    if (Receive_buffer_local_node_numbers != nullptr) deallocate(Receive_buffer_local_node_numbers, Nreceive_buffer);

    if (VertexTags != nullptr)
    {
        for (int n = 0; n < Nv; n++)
        {
            if (VertexTags[n] == nullptr)
                continue;
            
            // The std::vector inside is destroyed by 'delete', releasing its storage.
            delete VertexTags[n];
        }

        deallocate(VertexTags, Nv);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}
