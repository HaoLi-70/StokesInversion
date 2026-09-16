
#pragma once

/*--------------------------------------------------------------------------------*/

#include "me_solver.h"
#include "input_reader.h"
#include "dream.h"

/*--------------------------------------------------------------------------------*/

extern int INIT_INV(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_PARA *params, STRUCT_DREAM *dream, STRUCT_MPI *mpi, 
    STRUCT_SUBSET *subset);

/*--------------------------------------------------------------------------------*/
