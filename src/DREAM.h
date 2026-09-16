
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>

#include "allocation.h"
#include "parallel_runtime.h"
#include "me_solver.h"

/*--------------------------------------------------------------------------------*/

#define DREAM_HISTORY_MAX 1501

/*--------------------------------------------------------------------------------*/

typedef enum Sample_Output_Mode{
    SAMPLE_OUTPUT_NONE,
    SAMPLE_OUTPUT_MAGNETIC,
    SAMPLE_OUTPUT_ALL
}SAMPLE_OUTPUT_MODE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Dream{

    STRUCT_TENSOR chains;
    STRUCT_MATRIX likelihood;
    
    int nchains, nparams, burnin_generations, sampling_generations;
    int max_pairs, ncr, history_size;
    int retained_generations;
    double max_memory_gb, max_sample_file_gb;
    bool update_crossover, update_proposal_noise;
    SAMPLE_OUTPUT_MODE sample_output;

    int *jump_dims;
    double *diff, *scale_noise, *additive_noise;

    //CR
    double *crossover, *probabilities, *delta, *delta_total, *delta_sum;
    int *counts, *counts_total, *counts_sum;

}STRUCT_DREAM;

/*--------------------------------------------------------------------------------*/

struct Struct_Profile_IO;
struct Struct_Subset;

/*--------------------------------------------------------------------------------*/

extern int DREAM(STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params, 
    STRUCT_STK *stokes, struct Struct_Profile_IO *input,
    const struct Struct_Subset *subset);

extern int INIT_DREAM(STRUCT_MPI *mpi, STRUCT_DREAM *dream, 
    STRUCT_PARA *params, STRUCT_STK *stokes);

extern int Free_Dream(STRUCT_DREAM *dream);

/*--------------------------------------------------------------------------------*/
