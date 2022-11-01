
#ifndef DREAM_h
#define DREAM_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>

#include "ALLOCATION.h"
#include "RANDOM_NUMBER.h"
#include "GEMC.h"
#include "LIKELIHOOD.h"
#include "MPI_CONTROL.h"
#include "SORT.h"
#include "MATH_TOOL.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Dream{
  double ***Chains;
  double **Likelihood;
}STRUCT_DREAM;

typedef struct Struct_Cr{
  int Num_Cr;
  double *Cr, *Prob, *Delta, *Delta_tot, *Delta_sum;
  int *counts, *counts_tot, *counts_sum;
}STRUCT_CR;

/*--------------------------------------------------------------------------------*/

extern int DREAM(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, STRUCT_DREAM *Dream, \
    STRUCT_ATOM *Atom, STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par);

extern int Chain_Init_Dream(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_DREAM *Dream, STRUCT_ATOM *Atom, STRUCT_OBSERVATION *Observation, \
    STRUCT_PAR *Str_Par);

extern int GEMC2DREAM(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, \
    STRUCT_DREAM *Dream, STRUCT_PAR *Str_Par, STRUCT_OBSERVATION *Observation, \
    STRUCT_ATOM *Atom);

extern int Dream_Sample(STRUCT_MPI *Mpi, STRUCT_CR *Cr, STRUCT_PAR *Str_Par, \
    int Num_Pair, double **Chains, int Indx_Chain, int indx_Gener, \
    double *Sample);

extern int Sample_PairNum(int MaxNum_Pair, STRUCT_MPI *Mpi);

extern int Init_Cr(STRUCT_CR *Cr);

extern int Cr_Prob(STRUCT_CR *Cr, STRUCT_MPI *Mpi);

extern int Cr_distance(double ***Chains, int Num_Chain, int Num_Par, \
    int Indx_Chain, int Indx_Gener, int Indx_Cr, STRUCT_CR *Cr);

extern int Sample_Cr(STRUCT_CR *Cr, STRUCT_MPI *Mpi);

extern int Dream_Dim(STRUCT_MPI *Mpi, STRUCT_CR *Cr, STRUCT_PAR *Str_Par, \
    int Indx_Cr, int *Jump_dim);

extern int Dream_Diff(STRUCT_MPI *Mpi, double **Chains, int Indx_Chain, \
    int Num_Pair, int Num_Jump, int *Jump_dim, double *Diff);

extern int Rm_Outlierchain(STRUCT_MPI *Mpi, STRUCT_DREAM *Dream, int Num_Par, \
    int Gener);

/*--------------------------------------------------------------------------------*/

#endif /* DREAM_h */
