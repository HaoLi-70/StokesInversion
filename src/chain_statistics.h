
#pragma once

/*--------------------------------------------------------------------------------*/

#include "parallel_runtime.h"

/*--------------------------------------------------------------------------------*/

extern int PSRF(STRUCT_MPI *mpi, double ***chains, int generation, 
    const int *param_indices, int nparams, int chain_nparams,
    int history_size, double *rhat);

extern int Chains_STD(STRUCT_MPI *mpi, int begin_generation, 
    int end_generation, int nparams, double ***chains, int history_size, 
    double *stddev, double *mean);

extern int Chains_STD_Single(int nchains, int begin_generation,
    int end_generation, int nparams, double ***chains, int history_size,
    double *stddev, double *mean);

/*--------------------------------------------------------------------------------*/
