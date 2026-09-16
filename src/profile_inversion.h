#pragma once

/*--------------------------------------------------------------------------------*/

#include "dream.h"

/*--------------------------------------------------------------------------------*/

extern bool INVERSION_RESULT_FINITE(const double *values, int count);

extern int INVERT_PROFILE(STRUCT_MPI *mpi, STRUCT_PROFILE_IO *input,
    STRUCT_STK *stokes, STRUCT_PARA *params, STRUCT_DREAM *dream,
    STRUCT_SUBSET *subset);

/*--------------------------------------------------------------------------------*/
