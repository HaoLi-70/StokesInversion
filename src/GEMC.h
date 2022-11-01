
#ifndef GEMC_h
#define GEMC_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>

#include "ALLOCATION.h"
#include "RANDOM_NUMBER.h"
#include "LIKELIHOOD.h"
#include "MPI_CONTROL.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Gemc{
  double **Chains;
  double *Likelihood;
  double *Best_samples;
  double Best_likelihood;
}STRUCT_GEMC;

/*--------------------------------------------------------------------------------*/

extern int Chain_Init_Gemc(STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, STRUCT_ATOM *Atom, \
    STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par);

extern int BESTFIT_CHAIN(STRUCT_MPI *Mpi, double *Likelihood, int Num_Chains, \
    double **Chains, int Num_Par, double *Bestliklihood, double *Bestsample);

extern int GEMC(STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, STRUCT_ATOM *Atom, \
    STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par);

extern int Bounds_Enforce(double *Sample, STRUCT_PAR *Str_Par);

extern double Chain_Varance(STRUCT_MPI *Mpi, double *Likelihood, int Num_Chains);

/*--------------------------------------------------------------------------------*/

#endif /* GEMC_h */
