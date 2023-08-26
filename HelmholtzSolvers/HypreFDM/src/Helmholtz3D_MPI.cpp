#include "Helmholtz3D_MPI.h"

int Rank, size;

MPI_Datatype BackFront, BottomTop, LeftRight, BackFrontInt, BottomTopInt, LeftRightInt;

MPI_Request *request;
MPI_Status *status;

int Nx, Ny, Nz, Npx, Npy, Npz, offset = 1;

double xl, xr, yl, yr, zl, zr, Lx, Ly, Lz, *x, *y, *z, dx, dy, dz;

double alphal, betal, alphar, betar, Lambda;
double alphab, betab, alphat, betat;
double alphaback, betaback, alphaf, betaf;

double **leftBC, **rightBC, **bottomBC, **topBC, **backBC, **frontBC, ***rhsPDE, residual;

bool PeriodicBackFront, PeriodicBottomTop, PeriodicLeftRight;

int I, J, K, BackID, FrontID, BottomID, TopID, LeftID, RightID;

//_______________________________________________________________________________
// Hypre variables
// Index : Variable for global numbering of unknowns
//_______________________________________________________________________________
int ***Index, iterations, solverID = 0;

HYPRE_IJMatrix Ahypre;
HYPRE_ParCSRMatrix Acsr;

HYPRE_IJVector xhypre, bhypre;
HYPRE_ParVector xcsr, bcsr;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void SetSubDomains ( int &Npx, 
                     int &Npy, 
                     int &Npz, 
                     int &Nx, 
                     int &Ny, 
                     int &Nz )
{
    bool reverse = false;
    double ratio = Lx/(Ly+Lz);
    
    Npx = size;
    
    int Npyz = 1;
    
    if (Lx < (Ly+Lz))
    {
        reverse = true;
        ratio = (Ly+Lz)/Lx;
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
                Npyz = i;
            }
        }
    }
    
    if (reverse)
        Swap(Npx,Npyz);
    
    reverse = false;
    ratio = Ly/Lz;
    
    Npy = Npyz;
    Npz = 1;
    
    if (Ly < Lz)
    {
        reverse = true;
        ratio = Lz/Ly;
    }
    
    imax = floor(sqrt(double(Npyz)));
    distance = Absolute(ratio-Npyz);
    
    for (int i{2}; i <= imax; ++i)
    {
        if (Npyz % i == 0)
        {
            if (Absolute(ratio-Npyz/i) < distance)
            {
                distance = Absolute(ratio-Npyz/i);
                Npy = Npyz/i;
                Npz = i;
            }
        }
    }
    
    if (reverse)
        Swap(Npy,Npz);
    
    Nx = ceil(double(Nx)/Npx);
    Ny = ceil(double(Ny)/Npy);
    Nz = ceil(double(Nz)/Npz);
}

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( int Nx, 
                     int Ny, 
                     int Nz, 
                     double xl, 
                     double xr, 
                     double yl, 
                     double yr, 
                     double zl, 
                     double zr, 
                     bool PeriodicBackFront,
                     bool PeriodicBottomTop,
                     bool PeriodicLeftRight )
{
    ::xl = xl;
    ::xr = xr;
    ::yl = yl;
    ::yr = yr;
    ::zl = zl;
    ::zr = zr;
    ::PeriodicBackFront = PeriodicBackFront;
    ::PeriodicBottomTop = PeriodicBottomTop;
    ::PeriodicLeftRight = PeriodicLeftRight;
    
    Lx = (xr-xl);
    Ly = (yr-yl);
    Lz = (zr-zl);
    
    SetSubDomains(Npx,Npy,Npz,Nx,Ny,Nz);
    
    ::Nx = Nx;
    ::Ny = Ny;
    ::Nz = Nz;
    
    if ( (Nx < 3) || (Ny < 3) || (Nz < 3) )
    {
        std::cout << "Number of points in each direction has to be at least 3" << std::endl;
        
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    I = (Rank % (Npx*Npy)) % Npx;
    J = (Rank % (Npx*Npy)) / Npx;
    K =  Rank / (Npx*Npy);
    
    LeftID =   (I == 0       ? (PeriodicLeftRight ? (K*Npx*Npy + J*Npx + (Npx-1))   : -1) : (K*Npx*Npy        + J*Npx     + I - 1));
    RightID =  (I == (Npx-1) ? (PeriodicLeftRight ? (K*Npx*Npy + J*Npx)             : -1) : (K*Npx*Npy        + J*Npx     + I + 1));
    BottomID = (J == 0       ? (PeriodicBottomTop ? (K*Npx*Npy + (Npy-1)*Npx + I)   : -1) : (K*Npx*Npy        + (J-1)*Npx + I));
    TopID =    (J == (Npy-1) ? (PeriodicBottomTop ? (K*Npx*Npy + I)                 : -1) : (K*Npx*Npy        + (J+1)*Npx + I));
    BackID =   (K == 0       ? (PeriodicBackFront ? ((Npz-1)*Npx*Npy + J*Npx + I)   : -1) : ((K-1)*Npx*Npy + J*Npx        + I));
    FrontID =  (K == (Npz-1) ? (PeriodicBackFront ? (J*Npx + I)                     : -1) : ((K+1)*Npx*Npy + J*Npx        + I));
    
    // Variables for MPI communication
    int *BlockLength, *Displacement;
    
    Allocate(BlockLength,(Nx+2*offset)*(Ny+2*offset));
    Allocate(Displacement,(Nx+2*offset)*(Ny+2*offset));
    
    // LeftRight
    for (int i{}; i < Ny*offset; ++i)
        BlockLength[i] = Nz;
    
    for (int i{}; i < offset; ++i)
        for (int j{}; j < Ny; ++j)
            Displacement[i*Ny+j] = i*(Nz+2*offset)*(Ny+2*offset) + j*(Nz+2*offset);
    
    MPI_Type_indexed(Ny*offset, BlockLength, Displacement, MPI_DOUBLE, &LeftRight);
    MPI_Type_commit(&LeftRight);
    
    MPI_Type_indexed(Ny*offset, BlockLength, Displacement, MPI_INT, &LeftRightInt);
    MPI_Type_commit(&LeftRightInt);
    
    // BottomTop
    for (int i{}; i < (Nx+2*offset)*offset; ++i)
        BlockLength[i] = Nz;
    
    for (int i{}; i < (Nx+2*offset); ++i)
        for (int j{}; j < offset; ++j)
            Displacement[i*offset+j] = i*(Nz+2*offset)*(Ny+2*offset) + j*(Nz+2*offset);
    
    MPI_Type_indexed((Nx+2*offset)*offset, BlockLength, Displacement, MPI_DOUBLE, &BottomTop);
    MPI_Type_commit(&BottomTop);
    
    MPI_Type_indexed((Nx+2*offset)*offset, BlockLength, Displacement, MPI_INT, &BottomTopInt);
    MPI_Type_commit(&BottomTopInt);
    
    // BackFront
    for (int i{}; i < (Nx+2*offset)*(Ny+2*offset); ++i)
        BlockLength[i] = offset;
    
    for (int i{}; i < (Nx+2*offset); ++i)
        for (int j{}; j < (Ny+2*offset); ++j)
            Displacement[i*Ny+j] = i*(Nz+2*offset)*(Ny+2*offset) + j*(Nz+2*offset);
    
    MPI_Type_indexed((Nx+2*offset)*(Ny+2*offset), BlockLength, Displacement, MPI_DOUBLE, &BackFront);
    MPI_Type_commit(&BackFront);
    
    MPI_Type_indexed((Nx+2*offset)*(Ny+2*offset), BlockLength, Displacement, MPI_INT, &BackFrontInt);
    MPI_Type_commit(&BackFrontInt);
    
    Deallocate(BlockLength,(Nx+2*offset)*(Ny+2*offset));
    Deallocate(Displacement,(Nx+2*offset)*(Ny+2*offset));
    
    request = new (std::nothrow) MPI_Request [12];
    status  = new (std::nothrow) MPI_Status  [12];
    
    Allocate(Index,Nx+2*offset,Ny+2*offset,Nz+2*offset,offset,offset,offset);
    Allocate(x,Nx+2*offset,offset);
    Allocate(y,Ny+2*offset,offset);
    Allocate(z,Nz+2*offset,offset);
    
    for (int k{}; k < Nz; ++k)
        for (int j{}; j < Ny; ++j)
            for (int i{}; i < Nx; ++i)
                Index[i][j][k] = Rank*Nx*Ny*Nz + k*Nx*Ny + j*Nx + i;
    
    ExchangeData(Index);
    
    // Grid points and grid metric
    dx = Lx/(PeriodicLeftRight ? Nx*Npx : (Nx*Npx-1.0));
    dy = Ly/(PeriodicBottomTop ? Ny*Npy : (Ny*Npy-1.0));
    dz = Lz/(PeriodicBackFront ? Nz*Npz : (Nz*Npz-1.0));
    
    for (int i{-offset}; i < (Nx + offset); ++i)
        x[i] = xl + (Nx*I + i)*dx;
    
    for (int j{-offset}; j < (Ny + offset); ++j)
        y[j] = yl + (Ny*J + j)*dy;
    
    for (int k{-offset}; k < (Nz + offset); ++k)
        z[k] = zl + (Nz*K + k)*dz;
    
    Allocate(rhsPDE,Nx,Ny,Nz);
    Allocate(backBC,Nx,Ny);
    Allocate(frontBC,Nx,Ny);
    Allocate(bottomBC,Nx,Nz);
    Allocate(topBC,Nx,Nz);
    Allocate(leftBC,Ny,Nz);
    Allocate(rightBC,Ny,Nz);
    
    int start = Index[0][0][0], end = Index[Nx-1][Ny-1][Nz-1];
    
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
    
    Deallocate(backBC,Nx,Ny);
    Deallocate(frontBC,Nx,Ny);
    Deallocate(bottomBC,Nx,Nz);
    Deallocate(topBC,Nx,Nz);
    Deallocate(leftBC,Ny,Nz);
    Deallocate(rightBC,Ny,Nz);
    Deallocate(rhsPDE,Nx,Ny,Nz);
    
    Deallocate(Index,Nx+2*offset,Ny+2*offset,Nz+2*offset,offset,offset,offset);
    Deallocate(x,Nx+2*offset,offset);
    Deallocate(y,Ny+2*offset,offset);
    Deallocate(z,Nz+2*offset,offset);
    
    delete [] request;
    delete [] status;
}

//_______________________________________________________________________________
// Set A matrix
//_______________________________________________________________________________
void SetA ()
{
    HYPRE_IJMatrixInitialize(Ahypre);
    
    int n = 7, row, columns[n];
    
    double values[n], cx0 = 1.0/dx, cx1 = 1.0/(dx*dx), cy0 = 1.0/dy, cy1 = 1.0/(dy*dy), cz0 = 1.0/dz, cz1 = 1.0/(dz*dz);
    
    for (int k{}; k < Nz; ++k)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                if ((k == 0) && (BackID == -1))
                {
                    columns[0] = Index[i][j][k];
                    columns[1] = Index[i][j][k+1];
                    columns[2] = Index[i][j][k+2];
                    
                    values[0] = 1.5*alphaback*cz0 + betaback;
                    values[1] = -2.0*alphaback*cz0;
                    values[2] = 0.5*alphaback*cz0;
                    
                    n = 3;
                }
                else if ((k == (Nz-1)) && (FrontID == -1))
                {
                    columns[0] = Index[i][j][k];
                    columns[1] = Index[i][j][k-1];
                    columns[2] = Index[i][j][k-2];
                    
                    values[0] = 1.5*alphaf*cz0 + betaf;
                    values[1] = -2.0*alphaf*cz0;
                    values[2] = 0.5*alphaf*cz0;
                    
                    n = 3;
                }
                else
                {
                    if ((j == 0) && (BottomID == -1))
                    {
                        columns[0] = Index[i][j][k];
                        columns[1] = Index[i][j+1][k];
                        columns[2] = Index[i][j+2][k];
                        
                        values[0] = 1.5*alphab*cy0 + betab;
                        values[1] = -2.0*alphab*cy0;
                        values[2] = 0.5*alphat*cy0;
                        
                        n = 3;
                    }
                    else if ((j == (Ny-1)) && (TopID == -1))
                    {
                        columns[0] = Index[i][j][k];
                        columns[1] = Index[i][j-1][k];
                        columns[2] = Index[i][j-2][k];
                        
                        values[0] = 1.5*alphat*cy0 + betat;
                        values[1] = -2.0*alphat*cy0;
                        values[2] = 0.5*alphat*cy0;
                        
                        n = 3;
                    }
                    else
                    {
                        if ((i == 0) && (LeftID == -1))
                        {
                            columns[0] = Index[i][j][k];
                            columns[1] = Index[i+1][j][k];
                            columns[2] = Index[i+2][j][k];
                            
                            values[0] = 1.5*alphal*cx0 + betal;
                            values[1] = -2.0*alphal*cx0;
                            values[2] = 0.5*alphal*cx0;
                            
                            n = 3;
                        }
                        else if ((i == (Nx-1)) && (RightID == -1))
                        {
                            columns[0] = Index[i][j][k];
                            columns[1] = Index[i-1][j][k];
                            columns[2] = Index[i-2][j][k];
                            
                            values[0] = 1.5*alphar*cx0 + betar;
                            values[1] = -2.0*alphar*cx0;
                            values[2] = 0.5*alphar*cx0;
                            
                            n = 3;
                        }
                        else
                        {
                            columns[0] = Index[i-1][j][k];
                            columns[1] = Index[i+1][j][k];
                            columns[2] = Index[i][j-1][k];
                            columns[3] = Index[i][j+1][k];
                            columns[4] = Index[i][j][k-1];
                            columns[5] = Index[i][j][k+1];
                            columns[6] = Index[i][j][k];
                            
                            values[0] = cx1;
                            values[1] = cx1;
                            values[2] = cy1;
                            values[3] = cy1;
                            values[4] = cz1;
                            values[5] = cz1;
                            values[6] = -2.0*(cx1 + cy1 + cz1) - Lambda;
                            
                            n = 7;
                        }
                    }
                }
                
                row = Index[i][j][k];
                
                HYPRE_IJMatrixSetValues(Ahypre, 1, &n, &row, columns, values);
            }
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
    
    for (int k{}; k < Nz; ++k)
    {
        for (int j{}; j < Ny; ++j)
        {
            for (int i{}; i < Nx; ++i)
            {
                if ((k == 0) && (BackID == -1))
                    value = backBC[i][j];
                else if ((k == (Nz-1)) && (FrontID == -1))
                    value = frontBC[i][j];
                else
                {
                    if ((j == 0) && (BottomID == -1))
                        value = bottomBC[i][k];
                    else if ((j == (Ny-1)) && (TopID == -1))
                        value = topBC[i][k];
                    else
                    {
                        if ((i == 0) && (LeftID == -1))
                            value = leftBC[j][k];
                        else if ((i == (Nx-1)) && (RightID == -1))
                            value = rightBC[j][k];
                        else
                            value = rhsPDE[i][j][k];
                    }
                }
                
                HYPRE_IJVectorSetValues(bhypre, 1, &Index[i][j][k], &value);
            }
        }
    }
    
    HYPRE_IJVectorAssemble(bhypre);
    HYPRE_IJVectorGetObject(bhypre,(void **)&bcsr);
}

//_______________________________________________________________________________
// Solve the PDE
//_______________________________________________________________________________
void SolvePDE ( double ***&U )
{
    SetA();
    Setb();
    
    HYPRE_IJVectorInitialize(xhypre);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                HYPRE_IJVectorSetValues(xhypre, 1, &Index[i][j][k], &U[i][j][k]);
    
    HYPRE_IJVectorAssemble(xhypre);
    HYPRE_IJVectorGetObject(xhypre, (void **)&xcsr);
    
    // Solve Ax = b
    HypreSolve(Acsr,bcsr,xcsr,iterations,residual,solverID,1E-10,100,false);
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                HYPRE_IJVectorGetValues(xhypre, 1, &Index[i][j][k], &U[i][j][k]);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    ExchangeData(U);
}

//_______________________________________________________________________________
// First derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double FirstDerivativeX ( double ***U, 
                          const int i, 
                          const int j, 
                          const int k )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-0.5*U[i+2][j][k]+2.0*U[i+1][j][k]-1.5*U[i][j][k]) : (rightBoundary ? (1.5*U[i][j][k]-2.0*U[i-1][j][k]+0.5*U[i-2][j][k]) : 0.5*(U[i+1][j][k]-U[i-1][j][k])))/dx;
}

//_______________________________________________________________________________
// First derivative with respect to y (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double FirstDerivativeY ( double ***U, 
                          const int i, 
                          const int j, 
                          const int k )
{
    bool bottomBoundary  = ((BottomID  == -1) && (j == 0)  ? true : false);
    bool topBoundary     = ((TopID == -1) && (j == (Ny-1)) ? true : false);
    
    return (bottomBoundary ? (-0.5*U[i][j+2][k]+2.0*U[i][j+1][k]-1.5*U[i][j][k]) : (topBoundary ? (1.5*U[i][j][k]-2.0*U[i][j-1][k]+0.5*U[i][j-2][k]) : 0.5*(U[i][j+1][k]-U[i][j-1][k])))/dy;
}

//_______________________________________________________________________________
// First derivative with respect to z (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double FirstDerivativeZ ( double ***U, 
                          const int i, 
                          const int j, 
                          const int k )
{
    bool backBoundary  = ((BackID  == -1) && (k == 0)      ? true : false);
    bool frontBoundary = ((FrontID == -1) && (k == (Nz-1)) ? true : false);
    
    return (backBoundary ? (-0.5*U[i][j][k+2]+2.0*U[i][j][k+1]-1.5*U[i][j][k]) : (frontBoundary ? (1.5*U[i][j][k]-2.0*U[i][j][k-1]+0.5*U[i][j][k-2]) : 0.5*(U[i][j][k+1]-U[i][j][k-1])))/dz;
}

//_______________________________________________________________________________
// Second derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double SecondDerivativeX ( double ***U, 
                           const int i, 
                           const int j, 
                           const int k )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-U[i+3][j][k]+4.0*U[i+2][j][k]-5.0*U[i+1][j][k]+2.0*U[i][j][k])/dx : (rightBoundary ? (2.0*U[i][j][k]-5.0*U[i-1][j][k]+4.0*U[i-2][j][k]-U[i-3][j][k])/dx : (U[i+1][j][k]-2.0*U[i][j][k]+U[i-1][j][k])))/(dx*dx);
}

//_______________________________________________________________________________
// Second derivative with respect to y (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double SecondDerivativeY ( double ***U, 
                           const int i, 
                           const int j, 
                           const int k )
{
    bool bottomBoundary  = ((BottomID  == -1) && (j == 0)  ? true : false);
    bool topBoundary     = ((TopID == -1) && (j == (Ny-1)) ? true : false);
    
    return (bottomBoundary ? (-U[i][j+3][k]+4.0*U[i][j+2][k]-5.0*U[i][j+1][k]+2.0*U[i][j][k])/dy : (topBoundary ? (2.0*U[i][j][k]-5.0*U[i][j-1][k]+4.0*U[i][j-2][k]-U[i][j-3][k])/dy : (U[i][j+1][k]-2.0*U[i][j][k]+U[i][j-1][k])))/(dy*dy);
}

//_______________________________________________________________________________
// Second derivative with respect to z (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : (Nx+2*offset) x (Ny+2*offset) x (Nz+2*offset)
//_______________________________________________________________________________
double SecondDerivativeZ ( double ***U, 
                           const int i, 
                           const int j, 
                           const int k )
{
    bool backBoundary  = ((BackID  == -1) && (k == 0)      ? true : false);
    bool frontBoundary = ((FrontID == -1) && (k == (Nz-1)) ? true : false);
    
    return (backBoundary ? (-U[i][j][k+3]+4.0*U[i][j][k+2]-5.0*U[i][j][k+1]+2.0*U[i][j][k])/dz : (frontBoundary ? (2.0*U[i][j][k]-5.0*U[i][j][k-1]+4.0*U[i][j][k-2]-U[i][j][k-3])/dz : (U[i][j][k+1]-2.0*U[i][j][k]+U[i][j][k-1])))/(dz*dz);
}

//_______________________________________________________________________________
// MPI communication for double data
//_______________________________________________________________________________
void ExchangeData ( double ***&U )
{
    int comm = 0;
    
    // Left face
    if (LeftID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[-i0-1][i1][i2] = U[Nx-i0-1][i1][i2];
    }
    else
    {
        if (LeftID != -1)
        {
            MPI_Isend(&U[0][0][0],       1, LeftRight, LeftID, 5, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][0][0], 1, LeftRight, LeftID, 4, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Right face
    if (RightID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[Nx+i0][i1][i2] = U[i0][i1][i2];
    }
    else
    {
        if (RightID != -1)
        {
            MPI_Isend(&U[Nx-offset][0][0], 1, LeftRight, RightID, 4, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[Nx][0][0],        1, LeftRight, RightID, 5, MPI_COMM_WORLD, &request[comm+1]);
            
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
                for (int i2{}; i2 < Nz; ++i2)
                    U[i0][-i1-1][i2] = U[i0][Ny-i1-1][i2];
    }
    else
    {
        if (BottomID != -1)
        {
            MPI_Isend(&U[-offset][0][0],       1, BottomTop, BottomID, 3, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][0], 1, BottomTop, BottomID, 1, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Top face
    if (TopID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[i0][Ny+i1][i2] = U[i0][i1][i2];
    }
    else
    {
        if (TopID != -1)
        {
            MPI_Isend(&U[-offset][Ny-offset][0], 1, BottomTop, TopID, 1, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][Ny][0],        1, BottomTop, TopID, 3, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
    
    comm = 0;
    
    // Back face
    if (BackID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{-offset}; i1 < (Ny+offset); ++i1)
                for (int i2{}; i2 < offset; ++i2)
                    U[i0][i1][-i2-1] = U[i0][i1][Nz-i2-1];
    }
    else
    {
        if (BackID != -1)
        {
            MPI_Isend(&U[-offset][-offset][0],       1, BackFront, BackID, 0, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][-offset], 1, BackFront, BackID, 2, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Front face
    if (FrontID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{-offset}; i1 < (Ny+offset); ++i1)
                for (int i2{}; i2 < offset; ++i2)
                    U[i0][i1][Nz+i2] = U[i0][i1][i2];
    }
    else
    {
        if (FrontID != -1)
        {
            MPI_Isend(&U[-offset][-offset][Nz-offset], 1, BackFront, FrontID, 2, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][Nz],        1, BackFront, FrontID, 0, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
}

//_______________________________________________________________________________
// MPI communication for int data
//_______________________________________________________________________________
void ExchangeData ( int ***&U )
{
    int comm = 0;
    
    // Left face
    if (LeftID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[-i0-1][i1][i2] = U[Nx-i0-1][i1][i2];
    }
    else
    {
        if (LeftID != -1)
        {
            MPI_Isend(&U[0][0][0],       1, LeftRightInt, LeftID, 5, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][0][0], 1, LeftRightInt, LeftID, 4, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Right face
    if (RightID == Rank)
    {
        for (int i0{}; i0 < offset; ++i0)
            for (int i1{}; i1 < Ny; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[Nx+i0][i1][i2] = U[i0][i1][i2];
    }
    else
    {
        if (RightID != -1)
        {
            MPI_Isend(&U[Nx-offset][0][0], 1, LeftRightInt, RightID, 4, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[Nx][0][0],        1, LeftRightInt, RightID, 5, MPI_COMM_WORLD, &request[comm+1]);
            
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
                for (int i2{}; i2 < Nz; ++i2)
                    U[i0][-i1-1][i2] = U[i0][Ny-i1-1][i2];
    }
    else
    {
        if (BottomID != -1)
        {
            MPI_Isend(&U[-offset][0][0],       1, BottomTopInt, BottomID, 3, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][0], 1, BottomTopInt, BottomID, 1, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Top face
    if (TopID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{}; i1 < offset; ++i1)
                for (int i2{}; i2 < Nz; ++i2)
                    U[i0][Ny+i1][i2] = U[i0][i1][i2];
    }
    else
    {
        if (TopID != -1)
        {
            MPI_Isend(&U[-offset][Ny-offset][0], 1, BottomTopInt, TopID, 1, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][Ny][0],        1, BottomTopInt, TopID, 3, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
    
    comm = 0;
    
    // Back face
    if (BackID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{-offset}; i1 < (Ny+offset); ++i1)
                for (int i2{}; i2 < offset; ++i2)
                    U[i0][i1][-i2-1] = U[i0][i1][Nz-i2-1];
    }
    else
    {
        if (BackID != -1)
        {
            MPI_Isend(&U[-offset][-offset][0],       1, BackFrontInt, BackID, 0, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][-offset], 1, BackFrontInt, BackID, 2, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Front face
    if (FrontID == Rank)
    {
        for (int i0{-offset}; i0 < (Nx+offset); ++i0)
            for (int i1{-offset}; i1 < (Ny+offset); ++i1)
                for (int i2{}; i2 < offset; ++i2)
                    U[i0][i1][Nz+i2] = U[i0][i1][i2];
    }
    else
    {
        if (FrontID != -1)
        {
            MPI_Isend(&U[-offset][-offset][Nz-offset], 1, BackFrontInt, FrontID, 2, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset][-offset][Nz],        1, BackFrontInt, FrontID, 0, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double ***U, 
                    double (*ExactFunction)(const double, const double, const double) )
{
    double Error, LocalError = 0.0;
    
    for (int i{}; i < Nx; ++i)
        for (int j{}; j < Ny; ++j)
            for (int k{}; k < Nz; ++k)
                LocalError = Maximum(LocalError,Absolute(U[i][j][k]-ExactFunction(x[i],y[j],z[k])));
    
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
void WriteFile ( double ***U, 
                 double (*ExactFunction)(const double, const double, const double) )
{
    char *s;
    
    Allocate(s,200);
    
    int Nwx = (RightID == -1 ? Nx : (Nx+1));
    int Nwy = (TopID   == -1 ? Ny : (Ny+1));
    int Nwz = (FrontID == -1 ? Nz : (Nz+1));
    
    #ifdef TECPLOT
    sprintf(s,"Output/Solution-%04d.tec",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    FileWrite << "TITLE = \"Helmholtz solver\"" << std::endl;
    FileWrite << "Variables = \"X\",\"Y\",\"Z\",\"U\",\"Exact\"" << std::endl;
    FileWrite << "Zone I = " << Nwx << ", J = " << Nwy << ", K = " << Nwz << ", DATAPACKING=POINT" << std::endl;
    
    for (int k{}; k < Nwz; ++k)
        for (int j{}; j < Nwy; ++j)
            for (int i{}; i < Nwx; ++i)
                FileWrite << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << U[i][j][k] << "\t" << ExactFunction(x[i],y[j],z[k]) << std::endl;
    
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
    FileWrite << "DIMENSIONS " << Nwx << " " << Nwy << " " << Nwz << std::endl; 
    FileWrite << "POINTS " << Nwx*Nwy*Nwz << " FLOAT" << std::endl;
    
    for (int k{}; k < Nwz; ++k)
        for (int j{}; j < Nwy; ++j)
            for (int i{}; i < Nwx; ++i)
                FileWrite << x[i] << "\t" << y[j] << "\t" << z[k] << std::endl;
    
    FileWrite << std::endl << "POINT_DATA " << Nwx*Nwy*Nwz << std::endl;
    FileWrite << "SCALARS Phi float" << std::endl << "LOOKUP_TABLE default" << std::endl;
    
    for (int k{}; k < Nwz; ++k)
        for (int j{}; j < Nwy; ++j)
            for (int i{}; i < Nwx; ++i)
                FileWrite << U[i][j][k] << std::endl;
    
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
