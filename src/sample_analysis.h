
#pragma once

/*--------------------------------------------------------------------------------*/

#include "dream.h"

/*--------------------------------------------------------------------------------*/

extern int SAMPLE_ANALY(STRUCT_MPI *mpi, STRUCT_PARA *params, 
    STRUCT_DREAM *dream, STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset);

extern int SAMPLE_FILE_INIT(STRUCT_PROFILE_IO *input, STRUCT_PARA *params, 
    STRUCT_DREAM *dream, STRUCT_MPI *mpi);

extern int SAMPLE_FILE_CLOSE(STRUCT_PROFILE_IO *input);

extern int SAMPLE_RECORD_RESET(STRUCT_PROFILE_IO *input, STRUCT_DREAM *dream, 
    STRUCT_MPI *mpi, const STRUCT_SUBSET *subset);

extern int SAMPLE_WRITE_BLOCK(STRUCT_PROFILE_IO *input, STRUCT_DREAM *dream,
    STRUCT_PARA *params, STRUCT_MPI *mpi, const STRUCT_SUBSET *subset,
    int first_generation, int generation_count);

/*--------------------------------------------------------------------------------*/
