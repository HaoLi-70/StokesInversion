#include "profile_inversion.h"
#include "logger.h"
#include "sample_analysis.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Synchronized DREAM inversion and result preparation for one observed
        Stokes profile.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

bool INVERSION_RESULT_FINITE(const double *values, int count){

    /*######################################################################
      Purpose:
        Check that every value in a result buffer is finite.
      Input parameters:
        values, values to inspect.
        count, number of values.
      Return:
        true when the buffer is valid and finite; otherwise false.
    ######################################################################*/

    if(!values || count < 0) return false;
    for(int i=0; i<count; i++){
      if(!isfinite(values[i])) return false;
    }

    return true;
}

/*--------------------------------------------------------------------------------*/

int INVERT_PROFILE(STRUCT_MPI *mpi, STRUCT_PROFILE_IO *input,
    STRUCT_STK *stokes, STRUCT_PARA *params, STRUCT_DREAM *dream,
    STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Run the synchronized DREAM inversion and sample analysis for one
        profile while preserving zero output buffers on failure.
      Input parameters:
        mpi, MPI island and random-number state.
        input, profile I/O and output buffers.
        stokes, observed Stokes profile and synthesis workspace.
        params, model parameter definitions and bounds.
        dream, DREAM controls and chain state.
        subset, coordinates of the current profile.
      Output parameters:
        input, populated result and asymmetric-error buffers on success.
      Return:
        0 on success; 1 when a profile is intentionally skipped; a negative
        value for initialization, sampling, analysis, or MPI failure.
    ######################################################################*/

    if(!mpi || !input || !stokes || !params || !dream || !subset
        || !input->result_buffer || !input->error_buffer) return -1;

    const int result_count = params->nmodel+1;
    const int error_count = params->nmodel*2;
    memset(input->result_buffer, 0,
        (size_t)result_count*sizeof(*input->result_buffer));
    memset(input->error_buffer, 0,
        (size_t)error_count*sizeof(*input->error_buffer));

    int status = SAMPLE_RECORD_RESET(input, dream, mpi, subset);
    if(status == 0) status = INIT_DREAM(mpi, dream, params, stokes);
    if(MPI_SYNC_STATUS(mpi->island_comm, &status) != 0) return -1;

    if(status == 0){
      int dream_status = DREAM(mpi, dream, params, stokes, input, subset);
      status = dream_status < 0 ? dream_status : 0;
      if(MPI_SYNC_STATUS(mpi->island_comm, &status) != 0) return -1;
    }

    if(status == 0){
      status = SAMPLE_ANALY(mpi, params, dream, input, subset);
      if(status != 0) status = -1;
      if(MPI_SYNC_STATUS(mpi->island_comm, &status) != 0) return -1;
    }

    if(status == 0 && (!INVERSION_RESULT_FINITE(input->result_buffer, result_count)
        || !INVERSION_RESULT_FINITE(input->error_buffer, error_count))) status = -1;
    if(MPI_SYNC_STATUS(mpi->island_comm, &status) != 0) return -1;

    if(status < 0){
      memset(input->result_buffer, 0,
          (size_t)result_count*sizeof(*input->result_buffer));
      memset(input->error_buffer, 0,
          (size_t)error_count*sizeof(*input->error_buffer));
      if(mpi->island_rank == 0){
        LOG_ERROR(ERR_LVL_WARNING, "INVERT_PROFILE", status == -2
            ? "non-finite value detected during DREAM inversion.\n"
            : "profile inversion or sample analysis failed.\n");
      }
    }

    return status;
}

/*--------------------------------------------------------------------------------*/
