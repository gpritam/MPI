#include "Helmholtz2D_MPI.h"

int Rank, size;

MPI_Request *request;
MPI_Status *status;

MPI_Datatype BottomTop, LeftRight, BottomTopInt, LeftRightInt;

int Nx, Ny, Npx, Npy, offset = 1;

double xl, xr, yl, yr, Lx, Ly, *x, *y, dx, dy;

double alphal, betal, alphar, betar, Lambda;
double alphab, betab, alphat, betat;

double *leftBC, *rightBC, *bottomBC, *topBC, **rhsPDE, residual;

bool PeriodicBottomTop, PeriodicLeftRight;

int I, J, BottomID, TopID, LeftID, RightID;

//_______________________________________________________________________________
// Hypre variables
// Index : Variable for global numbering of unknowns
//_______________________________________________________________________________
int **Index, iterations, solverID = 0;

HYPRE_IJMatrix Ahypre;
HYPRE_ParCSRMatrix Acsr;

HYPRE_IJVector xhypre, bhypre;
HYPRE_ParVector xcsr, bcsr;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SetSubDomains ( int &Npx, 
                     int &Npy, 
                     int &Nx, 
                     int &Ny )
{
    bool reverse = false;
    double ratio = Lx/Ly;
    
    Npx = size;
    Npy = 1;
    
    if (Lx < Ly)
    {
        reverse = true;
        ratio = Ly/Lx;
    }
    
    int imax = floor(sqrt(double(size)));
    double distance = Absolute(ratio-size);
    
    for (int i{2}; i <= imax; ++i)
    {
        if (size % i == 0)
        {
            if (Absolute(ratio-size/i) < distance)
            {
                distance = Absolute(ratio-size/i);
                Npx = size/i;
                Npy = i;
            }
        }
    }
    
    if (reverse)
        Swap(Npx,Npy);
    
    Nx = ceil(double(Nx)/Npx);
    Ny = ceil(double(Ny)/Npy);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( int Nx, 
                     int Ny, 
                     double xl, 
                     double xr, 
                     double yl, 
                     double yr, 
                     bool PeriodicBottomTop, 
                     bool PeriodicLeftRight )
{
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    ::PeriodicBottomTop = PeriodicBottomTop;
    ::PeriodicLeftRight = PeriodicLeftRight;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    
    SetSubDomains(Npx,Npy,Nx,Ny);
    
    ::Nx = Nx;
    ::Ny = Ny;
    
    if ( (Nx < 3) || (Ny < 3) )
    {
        std::cout << "Number of points in each direction has to be at least 3" << std::endl;
        
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    // Find neighboring ranks
    I = Rank % Npx;
    J = Rank / Npx;
    
    LeftID   = (I == 0       ? (PeriodicLeftRight ? (J*Npx + Npx - 1) : -1) : (J*Npx     + I - 1));
    RightID  = (I == (Npx-1) ? (PeriodicLeftRight ? (J*Npx)           : -1) : (J*Npx     + I + 1));
    BottomID = (J == 0       ? (PeriodicBottomTop ? ((Npy-1)*Npx + I) : -1) : ((J-1)*Npx + I));
    TopID    = (J == (Npy-1) ? (PeriodicBottomTop ? I                 : -1) : ((J+1)*Npx + I));
    
    // Variables for MPI communication
    MPI_Type_vector(offset, Ny, (Ny+2*offset), MPI_DOUBLE, &LeftRight);
    MPI_Type_commit(&LeftRight);
    
    MPI_Type_vector((Nx+2*offset), offset, (Ny+2*offset), MPI_DOUBLE, &BottomTop);
    MPI_Type_commit(&BottomTop);
    
    MPI_Type_vector(offset, Ny, (Ny+2*offset), MPI_INT, &LeftRightInt);
    MPI_Type_commit(&LeftRightInt);
    
    MPI_Type_vector((Nx+2*offset), offset, (Ny+2*offset), MPI_INT, &BottomTopInt);
    MPI_Type_commit(&BottomTopInt);
    
    request = new (std::nothrow) MPI_Request [8];
    status  = new (std::nothrow) MPI_Status  [8];
    
    // Find 'Index' matrix
    Allocate(Index,Nx+2*offset,Ny+2*offset,offset,offset);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            Index[i][j] = Rank*Nx*Ny + j*Nx + i;
    
    ExchangeData(Index);
    
    // Grid points and grid metric
    dx = Lx/(PeriodicLeftRight ? Nx*Npx : (Nx*Npx-1.0));
    dy = Ly/(PeriodicBottomTop ? Ny*Npy : (Ny*Npy-1.0));
    
    Allocate(x,Nx+2*offset,offset);
    Allocate(y,Ny+2*offset,offset);
    
    for (int i{-offset}; i < (Nx + offset); ++i)
        x[i] = xl + (Nx*I + i)*dx;
    
    for (int j{-offset}; j < (Ny + offset); ++j)
        y[j] = yl + (Ny*J + j)*dy;
    
    Allocate(bottomBC,Nx);
    Allocate(topBC,Nx);
    Allocate(leftBC,Ny);
    Allocate(rightBC,Ny);
    Allocate(rhsPDE,Nx,Ny);
    
    int start = Index[0][0], end = Index[Nx-1][Ny-1];
    
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, start, end, start, end, &Ahypre);
    HYPRE_IJMatrixSetObjectType(Ahypre, HYPRE_PARCSR);
    
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, start, end, &xhypre);
    HYPRE_IJVectorSetObjectType(xhypre, HYPRE_PARCSR);
    
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, start, end, &bhypre);
    HYPRE_IJVectorSetObjectType(bhypre, HYPRE_PARCSR);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void DestroyPDEPlan ()
{
    HYPRE_IJMatrixDestroy(Ahypre);
    HYPRE_IJVectorDestroy(xhypre);
    HYPRE_IJVectorDestroy(bhypre);
    
    Deallocate(bottomBC,Nx);
    Deallocate(topBC,Nx);
    Deallocate(leftBC,Ny);
    Deallocate(rightBC,Ny);
    Deallocate(rhsPDE,Nx,Nx);
    
    Deallocate(Index,Nx+2*offset,Ny+2*offset,offset,offset);
    
    Deallocate(x,Nx+2*offset,offset);
    Deallocate(y,Ny+2*offset,offset);
    
    delete [] request;
    delete [] status;
}

//_______________________________________________________________________________
// Set A matrix
//_______________________________________________________________________________
void SetA ()
{
    HYPRE_IJMatrixInitialize(Ahypre);
    
    int n = 5, row, columns[n];
    
    double values[n], cx0 = 1.0/dx, cx1 = 1.0/(dx*dx), cy0 = 1.0/dy, cy1 = 1.0/(dy*dy);
    
    for (int j{}; j < Ny; j++)
    {
        for (int i{}; i < Nx; i++)
        {
            if ((j == 0) && (BottomID == -1))
            {
                columns[0] = Index[i][j];
                columns[1] = Index[i][j+1];
                columns[2] = Index[i][j+2];
                
                values[0] = 1.5*alphab*cy0 + betab;
                values[1] = -2.0*alphab*cy0;
                values[2] = 0.5*alphab*cy0;
                
                n = 3;
            }
            else if ((j == (Ny-1)) && (TopID == -1))
            {
                columns[0] = Index[i][j];
                columns[1] = Index[i][j-1];
                columns[2] = Index[i][j-2];
                
                values[0] = 1.5*alphat*cy0 + betat;
                values[1] = -2.0*alphat*cy0;
                values[2] = 0.5*alphat*cy0;
                
                n = 3;
            }
            else
            {
                if ((i == 0) && (LeftID == -1))
                {
                    columns[0] = Index[i][j];
                    columns[1] = Index[i+1][j];
                    columns[2] = Index[i+2][j];
                    
                    values[0] = 1.5*alphal*cx0 + betal;
                    values[1] = -2.0*alphal*cx0;
                    values[2] = 0.5*alphal*cx0;
                    
                    n = 3;
                }
                else if ((i == (Nx-1)) && (RightID == -1))
                {
                    columns[0] = Index[i][j];
                    columns[1] = Index[i-1][j];
                    columns[2] = Index[i-2][j];
                    
                    values[0] = 1.5*alphar*cx0 + betar;
                    values[1] = -2.0*alphar*cx0;
                    values[2] = 0.5*alphar*cx0;
                    
                    n = 3;
                }
                else
                {
                    columns[0] = Index[i-1][j];
                    columns[1] = Index[i+1][j];
                    columns[2] = Index[i][j-1];
                    columns[3] = Index[i][j+1];
                    columns[4] = Index[i][j];
                    
                    values[0] = cx1;
                    values[1] = cx1;
                    values[2] = cy1;
                    values[3] = cy1;
                    values[4] = -2.0*(cx1+cy1) - Lambda;
                    
                    n = 5;
                }
            }
            
            row = Index[i][j];
            
            HYPRE_IJMatrixSetValues(Ahypre, 1, &n, &row, columns, values);
        }
    }
    
    HYPRE_IJMatrixAssemble(Ahypre);
    HYPRE_IJMatrixGetObject(Ahypre, (void**)&Acsr);
}

//_______________________________________________________________________________
// Construct RHS
//_______________________________________________________________________________
void Setb ()
{
    double value;
    
    HYPRE_IJVectorInitialize(bhypre);
    
    for (int j{}; j < Ny; ++j)
    {
        for (int i{}; i < Nx; ++i)
        {
            if ((j == 0) && (BottomID == -1))
                value = bottomBC[i];
            else if ((j == (Ny-1)) && (TopID == -1))
                value = topBC[i];
            else
            {
                if ((i == 0) && (LeftID == -1))
                    value = leftBC[j];
                else if ((i == (Nx-1)) && (RightID == -1))
                    value = rightBC[j];
                else
                    value = rhsPDE[i][j];
            }
            
            HYPRE_IJVectorSetValues(bhypre, 1, &Index[i][j], &value);
        }
    }
    
    HYPRE_IJVectorAssemble(bhypre);
    HYPRE_IJVectorGetObject(bhypre,(void **)&bcsr);
}

//_______________________________________________________________________________
// Solve the PDE
//_______________________________________________________________________________
void SolvePDE ( double **&U )
{
    SetA();
    Setb();
    
    HYPRE_IJVectorInitialize(xhypre);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            HYPRE_IJVectorSetValues(xhypre, 1, &Index[i][j], &U[i][j]);
    
    HYPRE_IJVectorAssemble(xhypre);
    HYPRE_IJVectorGetObject(xhypre, (void **)&xcsr);
    
    // Solve Ax = b
    HypreSolve(Acsr,bcsr,xcsr,iterations,residual,solverID,1E-10,100,false);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            HYPRE_IJVectorGetValues(xhypre, 1, &Index[i][j], &U[i][j]);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    ExchangeData(U);
}

//_______________________________________________________________________________
// First derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset)
//_______________________________________________________________________________
double FirstDerivativeX ( double **U, 
                          const int i,
                          const int j )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-0.5*U[i+2][j]+2.0*U[i+1][j]-1.5*U[i][j]) : (rightBoundary ? (1.5*U[i][j]-2.0*U[i-1][j]+0.5*U[i-2][j]) : 0.5*(U[i+1][j]-U[i-1][j])))/dx;
}

//_______________________________________________________________________________
// First derivative with respect to y (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset)
//_______________________________________________________________________________
double FirstDerivativeY ( double **U, 
                          const int i, 
                          const int j )
{
    bool bottomBoundary = ((BottomID  == -1) && (j == 0)  ? true : false);
    bool topBoundary    = ((TopID == -1) && (j == (Ny-1)) ? true : false);
    
    return (bottomBoundary ? (-0.5*U[i][j+2]+2.0*U[i][j+1]-1.5*U[i][j]) : (topBoundary ? (1.5*U[i][j]-2.0*U[i][j-1]+0.5*U[i][j-2]) : 0.5*(U[i][j+1]-U[i][j-1])))/dy;
}

//_______________________________________________________________________________
// Second derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset)
//_______________________________________________________________________________
double SecondDerivativeX ( double **U, 
                           const int i, 
                           const int j )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-U[i+3][j]+4.0*U[i+2][j]-5.0*U[i+1][j]+2.0*U[i][j])/dx : (rightBoundary ? (2.0*U[i][j]-5.0*U[i-1][j]+4.0*U[i-2][j]-U[i-3][j])/dx : (U[i+1][j]-2.0*U[i][j]+U[i-1][j])))/(dx*dx);
}

//_______________________________________________________________________________
// Second derivative with respect to y (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset)
//_______________________________________________________________________________
double SecondDerivativeY ( double **U, 
                           const int i, 
                           const int j )
{
    bool bottomBoundary = ((BottomID  == -1) && (j == 0)  ? true : false);
    bool topBoundary    = ((TopID == -1) && (j == (Ny-1)) ? true : false);
    
    return (bottomBoundary ? (-U[i][j+3]+4.0*U[i][j+2]-5.0*U[i][j+1]+2.0*U[i][j])/dy : (topBoundary ? (2.0*U[i][j]-5.0*U[i][j-1]+4.0*U[i][j-2]-U[i][j-3])/dy : (U[i][j+1]-2.0*U[i][j]+U[i][j-1])))/(dy*dy);
}

//_______________________________________________________________________________
// MPI communication for double data
//_______________________________________________________________________________
void ExchangeData ( double **&U )
{
    int comm = 0;
    
    // Left face
    if (LeftID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                U[-i0-1][i1] = U[Nx-i0-1][i1];
    }
    else
    {
        if (LeftID != -1)
        {
            MPI_Isend(&U[0][0],       1, LeftRight, LeftID, 3, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][0], 1, LeftRight, LeftID, 1, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Right face
    if (RightID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                U[Nx+i0][i1] = U[i0][i1];
    }
    else
    {
        if (RightID != -1)
        {
            MPI_Isend(&U[Nx-offset][0], 1, LeftRight, RightID, 1, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[Nx][0],        1, LeftRight, RightID, 3, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
    
    comm = 0;
    
    // Bottom face
    if (BottomID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                U[i0][-i1-1] = U[i0][Ny-i1-1];
    }
    else
    {
        if (BottomID != -1)
        {
            MPI_Isend(&U[-offset][0],       1, BottomTop, BottomID, 0, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset], 1, BottomTop, BottomID, 2, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Top face
    if (TopID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                U[i0][Ny+i1] = U[i0][i1];
    }
    else
    {
        if (TopID != -1)
        {
            MPI_Isend(&U[-offset][Ny-offset], 1, BottomTop, TopID, 2, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][Ny],        1, BottomTop, TopID, 0, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
}

//_______________________________________________________________________________
// MPI communication for int data
//_______________________________________________________________________________
void ExchangeData ( int **&U )
{
    int comm = 0;
    
    // Left face
    if (LeftID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                U[-i0-1][i1] = U[Nx-i0-1][i1];
    }
    else
    {
        if (LeftID != -1)
        {
            MPI_Isend(&U[0][0],       1, LeftRightInt, LeftID, 3, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][0], 1, LeftRightInt, LeftID, 1, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Right face
    if (RightID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                U[Nx+i0][i1] = U[i0][i1];
    }
    else
    {
        if (RightID != -1)
        {
            MPI_Isend(&U[Nx-offset][0], 1, LeftRightInt, RightID, 1, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[Nx][0],        1, LeftRightInt, RightID, 3, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
    
    comm = 0;
    
    // Bottom face
    if (BottomID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                U[i0][-i1-1] = U[i0][Ny-i1-1];
    }
    else
    {
        if (BottomID != -1)
        {
            MPI_Isend(&U[-offset][0],       1, BottomTopInt, BottomID, 0, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset], 1, BottomTopInt, BottomID, 2, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Top face
    if (TopID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                U[i0][Ny+i1] = U[i0][i1];
    }
    else
    {
        if (TopID != -1)
        {
            MPI_Isend(&U[-offset][Ny-offset], 1, BottomTopInt, TopID, 2, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][Ny],        1, BottomTopInt, TopID, 0, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double **U, 
                    double (*ExactFunction)(const double, const double) )
{
    double Error, LocalError = 0.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            LocalError = Maximum(LocalError,Absolute(U[i][j]-ExactFunction(x[i],y[j])));
    
    if (size > 1)
        MPI_Allreduce(&LocalError, &Error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    else
        Error = LocalError;
    
    if (Rank == 0)
        std::cout << "Error = " << Error << std::endl << std::endl;
}

//_______________________________________________________________________________
// Write file in Tecplot format
//_______________________________________________________________________________
void WriteFile ( double **U, 
                 double (*ExactFunction)(const double, const double) )
{
    char *s;
    
    Allocate(s,200);
    
    int Nwx = (RightID == -1 ? Nx : (Nx+1));
    int Nwy = (TopID   == -1 ? Ny : (Ny+1));
    
    #ifdef TECPLOT
    sprintf(s,"Output/Solution-%04d.tec",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
        
    FileWrite << "TITLE = \"Helmholtz solver\"" << std::endl;
    FileWrite << "Variables = \"X\",\"Y\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nwx << ", J = " << Nwy << ", DATAPACKING=POINT" << std::endl;
    
    for (int j{}; j < Nwy; ++j)
        for (int i{}; i < Nwx; ++i)
            FileWrite << x[i] << "\t" << y[j] << "\t" << U[i][j] << "\t" << ExactFunction(x[i],y[j]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution-%04d.vtk",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "# vtk DataFile Version 3.1" << std::endl;
    FileWrite << "Helmholtz solver" << std::endl;
    FileWrite << "ASCII" << std::endl;
    FileWrite << "DATASET STRUCTURED_GRID" << std::endl;
    FileWrite << "DIMENSIONS " << Nwx << " " << Nwy << " " << 1 << std::endl; 
    FileWrite << "POINTS " << Nwx*Nwy << " FLOAT" << std::endl;
    
    for (int j{}; j < Nwy; ++j)
        for (int i{}; i < Nwx; ++i)
            FileWrite << x[i] << "\t" << y[j] << "\t" << 0.0 << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nwx*Nwy << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int j{}; j < Nwy; ++j)
        for (int i{}; i < Nwx; ++i)
            FileWrite << U[i][j] << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::out);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        Output << "!NBLOCKS " << size << std::endl;
        
        Output.close();
    }
    
    if (Rank == 0)
    {
        std::ofstream Output("Output/Field.visit", std::ios::app);
        
        if ( !Output )
            ErrorMessage("Output file couldnot be opened!");
        
        for (int i{}; i < size; ++i)
        {
            sprintf(s,"Solution-%04d.vtk",i);
            
            Output << s << std::endl;
        }
        
        Output.close();
    }
    #endif
    
    Deallocate(s,200);
}
