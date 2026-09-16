
#include "inversion_init.h"
#include "logger.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Validated initialization of inversion geometry, dynamic multi-line
        model layout, DREAM memory, numerical tables, and profile buffers.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/


static bool SIZE_PRODUCT(size_t left, size_t right, size_t *product){

    /*######################################################################
      Purpose:
        Multiply two allocation dimensions without overflowing size_t.
      Output parameters:
        product, multiplication result when representable.
      Return:
        true on success; false for a null output or overflow.
    ######################################################################*/

    if(!product || (right != 0 && left > SIZE_MAX/right)) return false;
    *product = left*right;
    return true;
}

/*--------------------------------------------------------------------------------*/

int INIT_INV(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_PARA *params, STRUCT_DREAM *dream, STRUCT_MPI *mpi, 
    STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Initialize model, profile, cache, and work buffers for inversion.
      Input parameters:
        input, runtime configuration and I/O buffers.
        stokes, wavelength grid and Stokes-profile state.
        params, model parameters, line data, and inversion bounds.
        mpi, MPI island and rank information.
        subset, profile subset assigned to one inversion task.
      Output parameters:
        input, with cache geometry and I/O buffers initialized.
        stokes, with synthesis ranges and Faddeeva data initialized.
        params, with line data and parameter metadata initialized.
        subset, with the initial profile coordinates and size.
      Return:
        0 on success; -1 for invalid configuration, memory-limit rejection,
          allocation failure, or spectral-line initialization failure.
    ######################################################################*/

    if(!input || !stokes || !params || !dream || !mpi || !subset) return -1;
    if(mpi->world_size < 1 || mpi->world_rank < 0
        || mpi->world_rank >= mpi->world_size || mpi->nislands < 1
        || stokes->nw < 1 || !stokes->wavelength || params->nlines < 1
        || !params->lines || !params->regions
        || dream->burnin_generations < 2 || dream->sampling_generations < 2
        || !isfinite(dream->max_memory_gb) || dream->max_memory_gb <= 0.0
        || input->nx < 1 || input->ny < 1) return -1;
    if(input->sol_box[0][0] < 0
        || input->sol_box[0][1] < input->sol_box[0][0]
        || input->sol_box[0][1] >= input->nx || input->sol_box[1][0] < 0
        || input->sol_box[1][1] < input->sol_box[1][0]
        || input->sol_box[1][1] >= input->ny) return -1;
    if((stokes->noise && stokes->noise_mode != NOISE_GLOBAL)
        || (!stokes->noise && stokes->noise_mode == NOISE_GLOBAL)
        || stokes->synthetic || input->result_buffer
        || input->error_buffer || stokes->faddeeva.exp_a2n2) return -1;

    const int report_rank = mpi->nislands > 1 ? 1 : 0;
    const bool report_init = mpi->world_rank == report_rank;

    #define INIT_VERBOSE(lvl, fmt, ...)                                  \
      do{                                                                \
        if(report_init && mpi->verbose_level >= (lvl)){                \
          snprintf(message_buffer, sizeof(message_buffer), fmt,          \
              ##__VA_ARGS__);                                            \
          LOG_WRITE(message_buffer, true, true);                         \
        }                                                                \
      }while(0)

    if(Model_Layout_Init(params, stokes, mpi) != 0) return -1;
    if(params->nmodel < 1 || params->nmodel > MAX_MODEL_PARAMS
        || params->npar < 1) return -1;
    for(int iregion=0; iregion<params->nregions; iregion++){
      const STRUCT_REGION *region = &params->regions[iregion];
      INIT_VERBOSE(3, " region %d spans wavelength indices %d:%d " 
          "(%.8g to %.8g).", iregion, region->iw_begin, region->iw_end, 
          region->wavelength_min, region->wavelength_max);
    }
    long double estimated_gib = 0.0L;
    if(!mpi->is_master){
      if(mpi->total_chains < 2) return -1;
      int max_generations = dream->burnin_generations
          > dream->sampling_generations ? dream->burnin_generations
          : dream->sampling_generations;
      int history_size = max_generations < DREAM_HISTORY_MAX
          ? max_generations : DREAM_HISTORY_MAX;
      long double estimated_bytes = (long double)history_size 
          *(long double)mpi->total_chains*(long double)(params->npar+1)
          *sizeof(double);
      estimated_gib = estimated_bytes/(1024.0L*1024.0L*1024.0L);
      long double memory_limit = (long double)dream->max_memory_gb 
          *1024.0L*1024.0L*1024.0L;
      if(estimated_bytes > memory_limit){
        if(report_init){
          snprintf(message_buffer, sizeof(message_buffer), 
              "Estimated DREAM memory per rank is %.6Lf GiB, exceeding " 
              "max_dream_memory_gb=%.3g. Reduce chains/free parameters or " 
              "increase the configured limit.\n", 
              estimated_gib, dream->max_memory_gb);
          LOG_ERROR(ERR_LVL_WARNING, "INIT_INV", message_buffer);
        }
        return -1;
      }
    }
    int line_status = Spectral_Lines_Init(stokes, params);
    if(line_status != 0){
      if(mpi->world_rank == 0){
        const char *error = line_status == -4 
            ? "Every configured spectral line must intersect at least one " 
              "observed wavelength region.\n" 
            : "Failed to initialize spectral-line wavelength ranges.\n";
        LOG_ERROR(ERR_LVL_ERROR, "INIT_INV", error);
      }
      return -1;
    }
    for(int iline=0; iline<params->nlines; iline++){
      const STRUCT_MELINE *line = &params->lines[iline];
      INIT_VERBOSE(3, " line %.8g active on wavelength indices %d:%d.",
          line->wavelength0, line->iw_begin, line->iw_end);
    }
    INIT_VERBOSE(2, " Inversion setup: %d wavelength points, %d region%s, "
        "%d line%s, %d free parameter%s, DREAM history %.6Lf GiB/rank.",
        stokes->nw, params->nregions, params->nregions == 1 ? "" : "s",
        params->nlines, params->nlines == 1 ? "" : "s", params->npar,
        params->npar == 1 ? "" : "s", estimated_gib);
    input->cache_header.nx = input->sol_box[0][1]
        -input->sol_box[0][0]+1;
    input->cache_header.ny = input->sol_box[1][1]
        -input->sol_box[1][0]+1;
    long long pixel_count = (long long)input->cache_header.nx
        *input->cache_header.ny;
    if(pixel_count < 1 || pixel_count > INT_MAX) return -1;
    input->counts = (int)pixel_count;

    size_t nbuf = 1;


    input->cache_header.ncache = input->counts;

    INIT_VERBOSE(4, " rank %d, total pixel number: %d, "
        "profiles per inversion: 1, buffer size: %zu.", mpi->world_rank,
        input->counts, nbuf);

    size_t profile_count, profile_bytes, result_count, error_count;
    if(!SIZE_PRODUCT((size_t)stokes->nw, 4U, &profile_count)
        || !SIZE_PRODUCT(profile_count, sizeof(double), &profile_bytes)
        || !SIZE_PRODUCT((size_t)params->nmodel+1U, nbuf, &result_count)
        || !SIZE_PRODUCT((size_t)params->nmodel, 2U, &error_count)
        || !SIZE_PRODUCT(error_count, nbuf, &error_count)) return -1;

    double *new_profile = NULL;
    double *new_noise = NULL;
    if(!stokes->profile){
      new_profile = malloc(profile_bytes);
    }
    if(stokes->noise_mode == NOISE_PER_PIXEL && !stokes->noise){
      new_noise = malloc(profile_bytes);
    }
    double *new_result = calloc(result_count, sizeof(*new_result));
    double *new_error = calloc(error_count, sizeof(*new_error));
    if((!stokes->profile && !new_profile)
        || (stokes->noise_mode == NOISE_PER_PIXEL
          && !stokes->noise && !new_noise)
        || !new_result || !new_error){
      free(new_profile);
      free(new_noise);
      free(new_result);
      free(new_error);
      return -1;
    }
    if(new_profile) stokes->profile = new_profile;
    if(new_noise) stokes->noise = new_noise;
    input->result_nparams = params->nmodel;
    input->result_buffer = new_result;
    input->error_buffer = new_error;

    subset->processed = 0;
    subset->coord[0] = input->sol_box[0][0];
    subset->coord[1] = input->sol_box[1][0];
    if(mpi->is_master) return 0;

    if(Faddeeva_init(&stokes->faddeeva) != 0) return -1;
    stokes->synthetic = malloc(profile_bytes);
    if(!stokes->synthetic) return -1;

    return 0;
}

#undef INIT_VERBOSE

/*--------------------------------------------------------------------------------*/
