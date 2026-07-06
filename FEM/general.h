#ifndef GENERAL_H
#define GENERAL_H

#include <iostream>
#include <new>
#include <cstdlib>
#include <utility>
#include <mpi.h>

#define SQRT2 1.4142135623730950
#define SQRT3 1.7320508075688772
#define PI 3.1415926535897932

//________________________________________________________________________________________________
// This function returns absolute value
//________________________________________________________________________________________________
template <class Type> Type absolute ( Type m )
{
	return (m >= Type(0) ? m : -m);
}

//________________________________________________________________________________________________
// This function swaps two values via a temporary. Works for any copy/move-
// assignable Type, unlike the arithmetic trick which fails for non-numeric
// types and overflows for large or floating-point inputs.
//________________________________________________________________________________________________
template <class Type> void swap ( Type &a,
								  Type &b )
{
	Type t = a;
	a = b;
	b = t;
}

//_______________________________________________________________________________
// Print an error message and exit the program. Marked inline so the header can
// be included from multiple translation units without a multi-definition link
// error.
//_______________________________________________________________________________
inline void print_error_message ( const char* message )
{
	std::cout << message << std::endl;

	std::exit(1);
}

//_______________________________________________________________________________
// Maximum
//_______________________________________________________________________________
template <class Type> inline Type maximum ( Type a0,
											Type a1 )
{
        return (a0 > a1 ? a0 : a1);
}

template <class Type> inline Type maximum ( Type a0,
											Type a1,
											Type a2 )
{
	return (a0 > a1 ? maximum(a0,a2) 
					: maximum(a1,a2));
}

//_______________________________________________________________________________
// Minimum
//_______________________________________________________________________________
template <class Type> inline Type minimum ( Type a0,
											Type a1 )
{
	return (a0 < a1 ? a0 : a1);
}

template <class Type> inline Type minimum ( Type a0,
											Type a1,
											Type a2 )
{
	return (a0 < a1 ? minimum(a0,a2) 
					: minimum(a1,a2));
}

//_______________________________________________________________________________
// Allocate an 1D array in contiguous memory locations. 
// Indexing starts with m[-offset1].
//_______________________________________________________________________________
template <class Type> void allocate ( Type *&m,
									  const int d1,
									  const int offset1 = 0 )
{
	m = new (std::nothrow) Type [d1];
	
	if (m == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
		m[i] = Type();
	
	m += offset1;
}

//_______________________________________________________________________________
// Deallocate an 1D array which has been allocated at contiguous memory locations. 
// Indexing starts with m[-offset1].
//_______________________________________________________________________________
template <class Type> void deallocate ( Type *&m,
										const int d1,
										const int offset1 = 0 )
{
	m -= offset1;
	
	delete [] m;
	
	m = nullptr;
}

//____________________________________________________________________________________
// Allocate a 2D matrix in contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2].
//____________________________________________________________________________________
template <class Type> void allocate ( Type **&m,
						              const int d1,
						              const int d2,
						              const int offset1 = 0,
						              const int offset2 = 0 )
{
	m = new (std::nothrow) Type* [d1];
	
	if (m == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	m[0] = new (std::nothrow) Type [d1*d2];
	
	if (m[0] == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
		m[i] = &(m[0][d2*i]);
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			m[i][j] = Type();
	
	m += offset1;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] += offset2;
}

//____________________________________________________________________________________
// Deallocate a 2D matrix which has been allocated at contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2].
//____________________________________________________________________________________
template <class Type> void deallocate ( Type **&m,
										 const int d1,
										 const int d2,
										 const int offset1 = 0,
										 const int offset2 = 0 )
{
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] -= offset2;
	
	m -= offset1;
	
	delete [] m[0];
	
	delete [] m;
	
	m = nullptr;
}

//____________________________________________________________________________________
// Allocate a 3D array in contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3].
//____________________________________________________________________________________
template <class Type> void allocate ( Type ***&m,
									  const int d1,
									  const int d2,
									  const int d3,
									  const int offset1 = 0,
									  const int offset2 = 0,
									  const int offset3 = 0 )
{
	m = new (std::nothrow) Type** [d1];
	
	if (m == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
	{
      m[i] = new (std::nothrow) Type* [d2];
		
		if (m[i] == 0)
			print_error_message("Error: Memory can not be allocated!");
	}
	
    m[0][0] = new (std::nothrow) Type [d1*d2*d3];
	
	if (m[0][0] == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			m[i][j] = &(m[0][0][i*d2*d3+j*d3]);
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			for (int k = 0; k < d3; k++)
				m[i][j][k] = Type();
	
	m += offset1;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] += offset2;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] += offset3;
}

//____________________________________________________________________________________
// Deallocate a 3D array which has been allocated at contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3].
//____________________________________________________________________________________
template <class Type> void deallocate ( Type ***&m,
										 const int d1,
										 const int d2,
										 const int d3,
										 const int offset1 = 0,
										 const int offset2 = 0,
										 const int offset3 = 0 )
{
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] -= offset3;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] -= offset2;
	
	m -= offset1;
	
	delete [] m[0][0];
	
	for (int i = 0; i < d1; i++)
		delete [] m[i];
	
	delete [] m;
	
	m = nullptr;
}

//____________________________________________________________________________________
// Allocate a 4D array in contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3][-offset4].
//____________________________________________________________________________________
template <class Type> void allocate ( Type ****&m,
									  const int d1,
									  const int d2,
									  const int d3,
									  const int d4,
									  const int offset1 = 0,
									  const int offset2 = 0,
									  const int offset3 = 0,
									  const int offset4 = 0 )
{
	m = new (std::nothrow) Type*** [d1];
	
	if (m == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
	{
		m[i] = new (std::nothrow) Type** [d2];
		
		if (m[i] == 0)
			print_error_message("Error: Memory can not be allocated!");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = new (std::nothrow) Type* [d3];
			
			if (m[i][j] == 0)
				print_error_message("Error: Memory can not be allocated!");
		}
	}
	
	m[0][0][0] = new (std::nothrow) Type [d1*d2*d3*d4];
	
	if (m[0][0][0] == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			for (int k = 0; k < d3; k++)
				m[i][j][k] = &(m[0][0][0][i*d2*d3*d4 + j*d3*d4 + k*d4]);
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			for (int k = 0; k < d3; k++)
				for (int l = 0; l < d4; l++)
					m[i][j][k][l] = Type();
	
	m += offset1;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] += offset2;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] += offset3;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				m[i][j][k] += offset4;
}

//____________________________________________________________________________________
// Deallocate a 4D array which has been allocated at contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3][-offset4].
//____________________________________________________________________________________
template <class Type> void deallocate ( Type ****&m,
										 const int d1,
										 const int d2,
										 const int d3,
										 const int d4,
										 const int offset1 = 0,
										 const int offset2 = 0,
										 const int offset3 = 0,
										 const int offset4 = 0 )
{
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				m[i][j][k] -= offset4;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] -= offset3;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] -= offset2;
	
	m -= offset1;
	
	delete [] m[0][0][0];
	
	for (int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
			delete [] m[i][j];
		
		delete [] m[i];
	}
	
	delete [] m;
	
	m = nullptr;
}

//____________________________________________________________________________________
// Allocate a 5D array in contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3][-offset4][-offset5].
//____________________________________________________________________________________
template <class Type> void allocate ( Type *****&m,
									  const int d1,
									  const int d2,
									  const int d3,
									  const int d4,
									  const int d5,
									  const int offset1 = 0,
									  const int offset2 = 0,
									  const int offset3 = 0,
									  const int offset4 = 0,
									  const int offset5 = 0 )
{
	m = new (std::nothrow) Type**** [d1];
	
	if (m == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
	{
		m[i] = new (std::nothrow) Type*** [d2];
		
		if (m[i]==0)
			print_error_message("Error: Memory can not be allocated!");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = new (std::nothrow) Type** [d3];
			
			if (m[i][j] == 0)
				print_error_message("Error: Memory can not be allocated!");

			for (int k = 0; k < d3; k++)
			{
				m[i][j][k] = new (std::nothrow) Type* [d4];
				
				if (m[i][j][k] == 0)
					print_error_message("Error: Memory can not be allocated!");
			}
		}
	}
	
	m[0][0][0][0] = new (std::nothrow) Type [d1*d2*d3*d4*d5];
	
	if (m[0][0][0][0] == 0)
		print_error_message("Error: Memory can not be allocated!");
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			for (int k = 0; k < d3; k++)
				for (int l = 0; l < d4; l++)
					m[i][j][k][l] = &(m[0][0][0][0][i*d2*d3*d4*d5 + j*d3*d4*d5 + k*d4*d5 + l*d5]);
	
	for (int i = 0; i < d1; i++)
		for (int j = 0; j < d2; j++)
			for (int k = 0; k < d3; k++)
				for (int l = 0; l < d4; l++)
					for (int p = 0; p < d5; p++)
						m[i][j][k][l][p] = Type();
	
	m += offset1;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] += offset2;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] += offset3;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				m[i][j][k] += offset4;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				for (int l = -offset4; l < (d4-offset4); l++)
					m[i][j][k][l] += offset5;
}

//____________________________________________________________________________________
// Deallocate a 5D array which has been allocated at contiguous memory locations. 
// Indexing starts with m[-offset1][-offset2][-offset3][-offset4][-offset5].
//____________________________________________________________________________________
template <class Type> void deallocate ( Type *****&m,
										 const int d1,
										 const int d2,
										 const int d3,
										 const int d4,
										 const int d5,
										 const int offset1 = 0,
										 const int offset2 = 0,
										 const int offset3 = 0,
										 const int offset4 = 0,
										 const int offset5 = 0 )
{
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				for (int l = -offset4; l < (d4-offset4); l++)
					m[i][j][k][l] -= offset5;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			for (int k = -offset3; k < (d3-offset3); k++)
				m[i][j][k] -= offset4;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		for (int j = -offset2; j < (d2-offset2); j++)
			m[i][j] -= offset3;
	
	for (int i = -offset1; i < (d1-offset1); i++)
		m[i] -= offset2;
	
	m -= offset1;
	
	delete [] m[0][0][0][0];
	
	for (int i = 0; i < d1; i++)
	{
		for (int j = 0; j < d2; j++)
		{
			for (int k = 0; k < d3; k++)
				delete [] m[i][j][k];
			
			delete [] m[i][j];
		}
		  
		delete [] m[i];
	}
	
	delete [] m;
	
	m = nullptr;
}
#endif