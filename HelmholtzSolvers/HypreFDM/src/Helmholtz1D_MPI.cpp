#include "Helmholtz1D_MPI.h"

int Rank, size;

MPI_Request *request;
MPI_Status *status;

int Nx, offset = 1;

double xl, xr, Lx, *x, dx;

double alphal, betal, alphar, betar, Lambda;

double leftBC, rightBC, *rhsPDE, residual;

bool PeriodicLeftRight;

int LeftID, RightID;

//_______________________________________________________________________________
// Hypre variables
// Index : Variable for global numbering of unknowns
//_______________________________________________________________________________
int *Index, iterations, solverID = 0;

HYPRE_IJMatrix Ahypre;
HYPRE_ParCSRMatrix Acsr;

HYPRE_IJVector xhypre, bhypre;
HYPRE_ParVector xcsr, bcsr;

//________________________________________________________________________________________________
// 
//________________________________________________________________________________________________
void CreatePDEPlan ( int Nx, 
                     double xl, 
                     double xr, 
                     bool PeriodicLeftRight )
{
    Nx = ceil(double(Nx)/size);
    
    ::Nx = Nx;
    ::xl = xl;
    ::xr = xr;
    ::PeriodicLeftRight = PeriodicLeftRight;
    
    Lx = (xr-xl);
    
    if (Nx < 3)
    {
        std::cout << "Minimum number of points should be 3" << std::endl;
        
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    LeftID =   (Rank == 0        ? (PeriodicLeftRight ? (size-1) : -1) : (Rank-1));
    RightID =  (Rank == (size-1) ? (PeriodicLeftRight ? 0         : -1) : (Rank+1));
    
    request = new (std::nothrow) MPI_Request [4];
    status  = new (std::nothrow) MPI_Status  [4];
    
    Allocate(Index,Nx+2*offset,offset);
    Allocate(x,Nx+2*offset,offset);
    
    for (int i{}; i < Nx; ++i)
        Index[i] = Nx*Rank + i;
    
    if (LeftID != -1)
        for (int i{1}; i <= offset; ++i)
            Index[-i] = Nx*(LeftID + 1) - i;
    
    if (RightID != -1)
        for (int i{0}; i < offset; ++i)
            Index[Nx+i] = Nx*RightID + i;
    
    // Grid points and grid metric
    dx = Lx/(PeriodicLeftRight ? Nx*size : (Nx*size-1.0));
    
    for (int i{-offset}; i < (Nx + offset); ++i)
        x[i] = xl + (Nx*Rank + i)*dx;
    
    Allocate(rhsPDE,Nx);
    
    int start = Index[0], end = Index[Nx-1];
    
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
    
    Deallocate(rhsPDE,Nx);
    
    Deallocate(Index,Nx+2*offset,offset);
    Deallocate(x,Nx+2*offset,offset);
    
    delete [] request;
    delete [] status;
}

//_______________________________________________________________________________
// Set A matrix
//_______________________________________________________________________________
void SetA ()
{
    HYPRE_IJMatrixInitialize(Ahypre);
    
    int n = 3, row, columns[n];
    
    double values[n], c0 = 1.0/dx, c1 = 1.0/(dx*dx);
    
    for (int i{}; i < Nx; ++i)
    {
        if ((i == 0) && (LeftID == -1))
        {
            columns[0] = Index[i];
            columns[1] = Index[i+1];
            columns[2] = Index[i+2];
            
            values[0] = 1.5*alphal*c0 + betal;
            values[1] = -2.0*alphal*c0;
            values[2] = 0.5*alphal*c0;
            
            n = 3;
        }
        else if ((i == (Nx-1)) && (RightID == -1))
        {
            columns[0] = Index[i];
            columns[1] = Index[i-1];
            columns[2] = Index[i-2];
            
            values[0] = 1.5*alphar*c0 + betar;
            values[1] = -2.0*alphar*c0;
            values[2] = 0.5*alphar*c0;
            
            n = 3;
        }
        else
        {
            columns[0] = Index[i-1];
            columns[1] = Index[i];
            columns[2] = Index[i+1];
            
            values[0] = c1;
            values[1] = -2.0*c1 - Lambda;
            values[2] = c1;
            
            n = 3;
        }
        
        row = Index[i];
        
        HYPRE_IJMatrixSetValues(Ahypre, 1, &n, &row, columns, values);
    }
    
    HYPRE_IJMatrixAssemble(Ahypre);
    HYPRE_IJMatrixGetObject(Ahypre, (void**)&Acsr);
}

//_______________________________________________________________________________
// Construct RHS
//_______________________________________________________________________________
void Setb ()
{
    bool leftBoundary, rightBoundary;
    
    double value;
    
    HYPRE_IJVectorInitialize(bhypre);
    
    for (int i{}; i < Nx; ++i)
    {
        leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
        rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
        
        value = (leftBoundary ? leftBC : (rightBoundary ? rightBC : rhsPDE[i]));
        
        HYPRE_IJVectorSetValues(bhypre, 1, &Index[i], &value);
    }
    
    HYPRE_IJVectorAssemble(bhypre);
    HYPRE_IJVectorGetObject(bhypre,(void **)&bcsr);
}

//_______________________________________________________________________________
// Solve the PDE
//_______________________________________________________________________________
void SolvePDE ( double *&U )
{
    SetA();
    Setb();
    
    HYPRE_IJVectorInitialize(xhypre);
    
    for (int i{}; i < Nx; ++i)
        HYPRE_IJVectorSetValues(xhypre, 1, &Index[i], &U[i]);
    
    HYPRE_IJVectorAssemble(xhypre);
    HYPRE_IJVectorGetObject(xhypre, (void **)&xcsr);
    
    // Solve Ax = b
    HypreSolve(Acsr,bcsr,xcsr,iterations,residual,solverID,1E-10,100,false);
    
    for (int i{}; i < Nx; ++i)
        HYPRE_IJVectorGetValues(xhypre, 1, &Index[i], &U[i]);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    ExchangeData(U);
}

//_______________________________________________________________________________
// First derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : Nx+2*offset
//_______________________________________________________________________________
double FirstDerivativeX ( double *U, 
                          const int i )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-0.5*U[i+2]+2.0*U[i+1]-1.5*U[i]) : (rightBoundary ? (1.5*U[i]-2.0*U[i-1]+0.5*U[i-2]) : 0.5*(U[i+1]-U[i-1])))/dx;
}

//_______________________________________________________________________________
// Second derivative with respect to x (WARNING! You can use this function only if 
// the data is padded from the neighboring ranks.)
// 
// U : Nx+2*offset
//_______________________________________________________________________________
double SecondDerivativeX ( double *U, 
                           const int i )
{
    bool leftBoundary  = ((LeftID  == -1) && (i == 0)      ? true : false);
    bool rightBoundary = ((RightID == -1) && (i == (Nx-1)) ? true : false);
    
    return (leftBoundary ? (-U[i+3]+4.0*U[i+2]-5.0*U[i+1]+2.0*U[i])/dx : (rightBoundary ? (2.0*U[i]-5.0*U[i-1]+4.0*U[i-2]-U[i-3])/dx : (U[i+1]-2.0*U[i]+U[i-1])))/(dx*dx);
}

//_______________________________________________________________________________
// MPI communication
//_______________________________________________________________________________
void ExchangeData ( double *&U )
{
    int comm = 0;
    
    // Left face
    if (LeftID == Rank)
        for (int i{1}; i <= offset; ++i)
            U[-i] = U[Nx-i];
    else
    {
        if (LeftID != -1)
        {
            MPI_Isend(&U[0],       offset, MPI_DOUBLE, LeftID, 3, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[-offset], offset, MPI_DOUBLE, LeftID, 1, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    // Right face
    if (RightID == Rank)
        for (int i{0}; i < offset; ++i)
            U[Nx+i] = U[i];
    else
    {
        if (RightID != -1)
        {
            MPI_Isend(&U[Nx-offset], offset, MPI_DOUBLE, RightID, 1, MPI_COMM_WORLD, &request[comm]);
            MPI_Irecv(&U[Nx],        offset, MPI_DOUBLE, RightID, 3, MPI_COMM_WORLD, &request[comm+1]);
            
            comm += 2;
        }
    }
    
    for (int i{}; i < comm; ++i)
        MPI_Wait(&request[i],&status[i]);
}

//_______________________________________________________________________________
// Compute error
//_______________________________________________________________________________
void ComputeError ( double *U, 
                    double (*ExactFunction)(const double) )
{
    double Error, LocalError = 0.0;
    
    for (int i{}; i < Nx; ++i)
        LocalError = Maximum(LocalError,Absolute(U[i]-ExactFunction(x[i])));
    
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
void WriteFile ( double *U, 
                 double (*ExactFunction)(const double) )
{
    char *s;
    
    Allocate(s,200);
    
    int Nwx = (RightID == -1 ? Nx : (Nx+1));
    
    #ifdef TECPLOT
    sprintf(s,"Output/Solution-%04d.tec",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    for (int i = 0; i < Nwx; i++)
        FileWrite << x[i] << "\t" << U[i] << "\t" << ExactFunction(x[i]) << std::endl;
    
    FileWrite.close();
    #endif
    
    #ifdef VISIT
    sprintf(s,"Output/Solution-%04d.curve",Rank);
    
    std::ofstream FileWrite(s, std::ios::out);
    FileWrite.flags( std::ios::dec | std::ios::fixed );
    FileWrite.precision(8);
    
    if ( !FileWrite )
        ErrorMessage("Output file couldnot be opened!");
    
    if (Rank == 0)
        FileWrite << "#Computed" << std::endl;
    
    for (int i = 0; i < Nwx; i++)
        FileWrite << x[i] << "\t" << U[i] << std::endl;
    
    FileWrite.close();
    #endif
    
    Deallocate(s,200);
}
