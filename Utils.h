/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#pragma once
#include <stdarg.h>
#include <stdio.h>
#include <string>

using std::string;

/**************************
 * Utility functions
 **************************
*/
//--------------------------------------------------
//! MPI Rank 0 will print to stdout
void print(const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to stdout in red color
void print_error(const char format[],...);
//--------------------------------------------------
//! MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...);
//--------------------------------------------------
//! Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime();
//--------------------------------------------------
//! Call MPI_Finalize and exit (with error)
void exit_mpi();
//--------------------------------------------------
//! Check for NAN
template <class T>
inline int m2c_isnan(const T& t) {return (t != t);}
//--------------------------------------------------
