
#ifndef MATH_TOOL_h
#define MATH_TOOL_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>

#include "ALLOCATION.h"
#include "MPI_CONTROL.h"

/*--------------------------------------------------------------------------------*/

extern int PSRF(STRUCT_MPI *Mpi, double ***Chains, int Gener, int Num_Par, \
    double *Statis_R);

extern int Chains_STD(STRUCT_MPI *Mpi, int Begin_Gener, int End_Gener, \
    int Num_Par, double ***Chains, double *STD, double *Mean);

extern int Chains_STD_Single(int Num_Chain, int Begin_Gener, int End_Gener, \
    int Num_Par, double ***Chains, double *STD, double *Mean);

/*--------------------------------------------------------------------------------*/

#endif /* MATH_TOOL_h */
