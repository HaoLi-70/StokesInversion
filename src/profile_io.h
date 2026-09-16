
#pragma once

/*--------------------------------------------------------------------------------*/

#include "me_solver.h"
#include "parallel_runtime.h"
#include "input_reader.h"
#include "dream.h"

/*--------------------------------------------------------------------------------*/

extern int READ_WAVELENGTH(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_MPI *mpi);

extern int READ_PROFILE(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_SUBSET *subset);

extern int PIXEL_ADVANCE(STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset);

extern int CACHE_INIT(STRUCT_PROFILE_IO *input, STRUCT_MPI *mpi, 
    const STRUCT_PARA *params, const STRUCT_DREAM *dream, 
    const STRUCT_STK *stokes);

extern int WRITE_RESULT(STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset,
    bool mark_complete);

extern void MODEL_PARAMETER_NAME(const STRUCT_PARA *params, int index, 
    char *name, size_t size);

extern int CLOSE_FILES(void);

/*--------------------------------------------------------------------------------*/
