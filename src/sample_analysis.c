
#include "sample_analysis.h"
#include "chain_statistics.h"
#include "logger.h"
#include "profile_io.h"
#include "sorting.h"
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

#define SAMPLE_HEADER_INTS (MAX_MODEL_PARAMS+23)
#define SAMPLE_MAGIC "SMPL"
#define SAMPLE_MAGIC_SIZE 4
#define SAMPLE_META_OFFSET (7+MAX_MODEL_PARAMS)

/*--------------------------------------------------------------------------------*/

static bool SAMPLE_PARAM_SELECTED(const STRUCT_DREAM *dream, 
    const STRUCT_PARA *params, int model_idx){

    /*######################################################################
      Purpose:
        Determine whether one model parameter is selected for sample output.
      Input parameters:
        dream, configured sample-output mode.
        params, model-parameter metadata.
        model_idx, model parameter index to test.
      Return:
        true when the parameter should be written; false otherwise.
    ######################################################################*/

    return dream->sample_output == SAMPLE_OUTPUT_ALL 
        || (dream->sample_output == SAMPLE_OUTPUT_MAGNETIC 
          && params->kind[model_idx] <= PARAM_B2);
}

/*--------------------------------------------------------------------------------*/

static long long SAMPLE_RECORD_OFFSET(const STRUCT_PROFILE_IO *input, 
    const STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Calculate the byte offset of one pixel record in the sample file.
      Input parameters:
        input, sample-file layout and inversion-box metadata.
        subset, coordinates of the current pixel.
      Return:
        Byte offset from the beginning of the sample file.
    ######################################################################*/

    long long pixel_index = (long long)(subset->coord[1] 
        -input->sol_box[1][0])*input->cache_header.nx 
        +subset->coord[0]-input->sol_box[0][0];

    return input->sample_header_size+pixel_index*input->sample_record_size;

}

/*--------------------------------------------------------------------------------*/

int SAMPLE_FILE_INIT(STRUCT_PROFILE_IO *input, STRUCT_PARA *params, 
    STRUCT_DREAM *dream, STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
        Create or validate the shared fixed-record sample file.
      Input parameters:
        input, paths, image dimensions, and cache-reuse state.
        params, model layout and parameter-selection metadata.
        dream, sample-output and history controls.
        mpi, MPI communicators and chain dimensions.
      Output parameters:
        input, populated sample-file layout and open file handle.
      Return:
        0 on success; -1 for an invalid layout or I/O failure.
    ######################################################################*/

    if(!input || !params || !dream || !mpi) return -1;
    if(dream->sample_output == SAMPLE_OUTPUT_NONE) return 0;
    if(sizeof(int) != 4 || sizeof(float) != 4 || sizeof(uint64_t) != 8){
      return -1;
    }

    int noutput = 0;
    int indices[MAX_MODEL_PARAMS];
    for(int imodel=0; imodel<params->nmodel; imodel++){
      if(params->inv[imodel] 
          && SAMPLE_PARAM_SELECTED(dream, params, imodel)){
        indices[noutput++] = imodel;
      }
    }

    int max_chains = mpi->total_chains;
    MPI_Allreduce(MPI_IN_PLACE, &max_chains, 1, MPI_INT, MPI_MAX, 
        mpi->world_comm);
    int sample_generations = dream->sampling_generations;
    if(max_chains < 1 || sample_generations < 1
        || max_chains > INT_MAX/sample_generations){
      return -1;
    }
    input->sample_max_count = sample_generations*max_chains;
    input->sample_header_size = (long long)SAMPLE_HEADER_INTS*sizeof(int);
    if(noutput > 0 && input->sample_max_count 
        > (LLONG_MAX-2LL*(long long)sizeof(int)) 
          /((long long)noutput*(long long)sizeof(float))){
      return -1;
    }
    input->sample_record_size = 2*(long long)sizeof(int) 
        +(long long)noutput*input->sample_max_count*(long long)sizeof(float);
    if(input->counts > 0 && input->sample_record_size 
        > (LLONG_MAX-input->sample_header_size)/input->counts){
      return -1;
    }
    long double estimated_bytes = (long double)input->sample_header_size
        +(long double)input->counts*(long double)input->sample_record_size;
    long double estimated_gib = estimated_bytes
        /(1024.0L*1024.0L*1024.0L);
    if(!isfinite(dream->max_sample_file_gb)
        || estimated_gib > (long double)dream->max_sample_file_gb){
      if(mpi->world_rank == 0){
        snprintf(message_buffer, sizeof(message_buffer),
            "Estimated sample file size is %.6Lf GiB, exceeding "
            "max_sample_file_gb=%.6g. Reduce sampling_generations, chains, "
            "pixels, or output parameters.\n", estimated_gib,
            dream->max_sample_file_gb);
        LOG_ERROR(ERR_LVL_ERROR, "SAMPLE_FILE_INIT", message_buffer);
      }
      return -1;
    }

    int header[SAMPLE_HEADER_INTS] = {0};
    memcpy(header, SAMPLE_MAGIC, SAMPLE_MAGIC_SIZE);
    header[1] = input->cache_header.nx;
    header[2] = input->cache_header.ny;
    header[3] = noutput;
    header[4] = input->sample_max_count;
    header[5] = 32;
    header[6] = params->nmodel;
    for(int i=0; i<noutput; i++) header[7+i] = indices[i];
    header[SAMPLE_META_OFFSET] = input->sol_box[0][0];
    header[SAMPLE_META_OFFSET+1] = input->sol_box[0][1];
    header[SAMPLE_META_OFFSET+2] = input->sol_box[1][0];
    header[SAMPLE_META_OFFSET+3] = input->sol_box[1][1];
    memcpy(header+SAMPLE_META_OFFSET+4, &input->config_hash, 
        sizeof(input->config_hash));
    header[SAMPLE_META_OFFSET+6] = params->nlines;
    header[SAMPLE_META_OFFSET+7] = params->nregions;
    header[SAMPLE_META_OFFSET+8] = params->magnetic_mode;

    int file_exists = 0;
    if(mpi->world_rank == 0) file_exists = FILE_EXIST(input->sample_path);
    if(MPI_Bcast(&file_exists, 1, MPI_INT, 0, mpi->world_comm)
        != MPI_SUCCESS) return -1;
    int status = MPI_File_open(mpi->world_comm, input->sample_path, 
        MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &input->sample_file);
    int local_open = status == MPI_SUCCESS ? 1 : 0;
    int all_open = 0;
    if(MPI_Allreduce(&local_open, &all_open, 1, MPI_INT, MPI_MIN,
        mpi->world_comm) != MPI_SUCCESS || !all_open){
      /* MPI_File_close is collective and therefore cannot safely close a
         handle that exists on only a subset of ranks.  Leave any partial
         handle to MPI_Finalize and keep the public state consistently
         closed on every rank. */
      input->sample_file_open = false;
      input->sample_file = MPI_FILE_NULL;
      return -1;
    }
    input->sample_file_open = true;
    int io_status = 0;
    MPI_Offset expected_size = (MPI_Offset)input->sample_header_size 
        +(MPI_Offset)input->counts*input->sample_record_size;
    if(!input->cache_reused || !file_exists){
      int truncate_status = MPI_File_set_size(input->sample_file, 0);
      int resize_status = MPI_File_set_size(input->sample_file, expected_size);
      if(truncate_status != MPI_SUCCESS || resize_status != MPI_SUCCESS){
        io_status = -1;
      }
      if(mpi->world_rank == 0){
        if(MPI_File_write_at(input->sample_file, 0, header, 
            SAMPLE_HEADER_INTS, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){
          io_status = -1;
        }
      }
    }else{
      int disk_header[SAMPLE_HEADER_INTS] = {0};
      if(mpi->world_rank == 0){
        if(MPI_File_read_at(input->sample_file, 0, disk_header, 
            SAMPLE_HEADER_INTS, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){
          disk_header[0] = 0;
        }
      }
      if(MPI_Bcast(disk_header, SAMPLE_HEADER_INTS, MPI_INT, 0,
          mpi->world_comm) != MPI_SUCCESS) io_status = -1;
      MPI_Offset disk_size = 0;
      if(MPI_File_get_size(input->sample_file, &disk_size) != MPI_SUCCESS){
        io_status = -1;
      }
      if(memcmp(disk_header, header, sizeof(header)) != 0 
          || disk_size != expected_size){
        if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_WARNING, 
            "SAMPLE_FILE_INIT", "sample file is inconsistent with this run.\n");
        io_status = -1;
      }
    }
    if(MPI_Allreduce(MPI_IN_PLACE, &io_status, 1, MPI_INT, MPI_MIN,
        mpi->world_comm) != MPI_SUCCESS) io_status = -1;
    if(io_status != 0){
      if(MPI_File_close(&input->sample_file) != MPI_SUCCESS) io_status = -1;
      input->sample_file_open = false;
      return -1;
    }
    if(MPI_Barrier(mpi->world_comm) != MPI_SUCCESS){
      MPI_File_close(&input->sample_file);
      input->sample_file_open = false;
      return -1;
    }

    return 0;

}

/*--------------------------------------------------------------------------------*/

int SAMPLE_RECORD_RESET(STRUCT_PROFILE_IO *input, STRUCT_DREAM *dream, 
    STRUCT_MPI *mpi, const STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Mark the current pixel sample record as incomplete before inversion.
      Input parameters:
        input, open sample file and record layout.
        dream, configured sample-output mode.
        mpi, MPI island communicator and local rank.
        subset, coordinates of the current pixel.
      Return:
        0 on success or when sample output is disabled; -1 on I/O failure.
    ######################################################################*/

    if(!input || !dream || !mpi || !subset) return -1;
    if(dream->sample_output == SAMPLE_OUTPUT_NONE) return 0;
    if(!input->sample_file_open) return -1;

    int io_status = 0;
    if(mpi->island_rank == 0){
      int record_header[2] = {0, 0};
      long long offset = SAMPLE_RECORD_OFFSET(input, subset);
      if(MPI_File_write_at(input->sample_file, (MPI_Offset)offset, 
          record_header, 2, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){
        io_status = -1;
      }
    }
    if(mpi->island_size > 1){
      MPI_Bcast(&io_status, 1, MPI_INT, 0, mpi->island_comm);
    }
    return io_status;
}

/*--------------------------------------------------------------------------------*/

int SAMPLE_WRITE_BLOCK(STRUCT_PROFILE_IO *input, STRUCT_DREAM *dream,
    STRUCT_PARA *params, STRUCT_MPI *mpi, const STRUCT_SUBSET *subset,
    int first_generation, int generation_count){

    /*######################################################################
      Purpose:
        Write one contiguous block of sampling generations before their
        circular-history slots are reused.
      Input parameters:
        input, open sample file and fixed-record layout.
        dream, sample-output controls and circular chain history.
        params, free-parameter mapping and output-selection metadata.
        mpi, MPI island layout.
        subset, coordinates of the current profile.
        first_generation, first sampling-generation index in the block.
        generation_count, number of consecutive generations to write.
      Return:
        0 on success or when sample output is disabled; -1 on invalid state,
        allocation failure, or MPI-I/O failure.
      Note:
        All island ranks call this function, but only the island leader writes
        because each generation has already been gathered across the island.
    ######################################################################*/

    if(!input || !dream || !params || !mpi || !subset
        || first_generation < 0 || generation_count < 1
        || first_generation > dream->sampling_generations-generation_count){
      return -1;
    }
    if(dream->sample_output == SAMPLE_OUTPUT_NONE) return 0;
    if(!input->sample_file_open || !dream->chains.data_ptr
        || mpi->total_chains < 1
        || generation_count > INT_MAX/mpi->total_chains) return -1;

    int block_samples = generation_count*mpi->total_chains;
    float *output = NULL;
    if(mpi->island_rank == 0){
      output = (float *)malloc((size_t)block_samples*sizeof(*output));
    }
    int io_status = mpi->island_rank != 0 || output ? 0 : -1;
    if(mpi->island_size > 1){
      if(MPI_Allreduce(MPI_IN_PLACE, &io_status, 1, MPI_INT, MPI_MIN,
          mpi->island_comm) != MPI_SUCCESS) io_status = -1;
    }
    if(io_status != 0){
      free(output);
      return -1;
    }

    if(mpi->island_rank == 0){
      double ***chains = TENSOR_DBL(dream->chains);
      long long record_offset = SAMPLE_RECORD_OFFSET(input, subset);
      int output_idx = 0;
      for(int param_idx=0; param_idx<params->npar; param_idx++){
        int model_idx = params->free_index[param_idx];
        if(!SAMPLE_PARAM_SELECTED(dream, params, model_idx)) continue;

        int sample_idx = 0;
        for(int generation=first_generation;
            generation<first_generation+generation_count; generation++){
          int slot = generation%dream->history_size;
          for(int chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
            output[sample_idx++] = (float)chains[slot][chain_idx][param_idx];
          }
        }
        long long data_offset = record_offset+2*(long long)sizeof(int)
            +(long long)output_idx*input->sample_max_count*sizeof(float)
            +(long long)first_generation*mpi->total_chains*sizeof(float);
        if(MPI_File_write_at(input->sample_file, (MPI_Offset)data_offset,
            output, block_samples, MPI_FLOAT, MPI_STATUS_IGNORE)
            != MPI_SUCCESS){
          io_status = -1;
          break;
        }
        output_idx++;
      }
    }

    free(output);
    if(mpi->island_size > 1){
      if(MPI_Bcast(&io_status, 1, MPI_INT, 0, mpi->island_comm)
          != MPI_SUCCESS) io_status = -1;
    }
    return io_status;
}

/*--------------------------------------------------------------------------------*/

int SAMPLE_FILE_CLOSE(STRUCT_PROFILE_IO *input){

    /*######################################################################
      Purpose:
        Close the sample file when it is open.
      Input parameters:
        input, sample-file handle and open-state flag.
      Output parameters:
        input, with the sample file marked closed.
      Return:
        0 on success or when already closed; -1 on close failure.
    ######################################################################*/

    if(!input || !input->sample_file_open) return 0;
    int status = MPI_File_close(&input->sample_file);
    input->sample_file_open = false;

    return status == MPI_SUCCESS ? 0 : -1;

}

/*--------------------------------------------------------------------------------*/

int SAMPLE_ANALY(STRUCT_MPI *mpi, STRUCT_PARA *params, STRUCT_DREAM *dream, 
    STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset){
    
    /*######################################################################
      Purpose:
        Summarize the retained DREAM window and commit the sample record.
      Input parameters:
        mpi, MPI island layout and chain distribution.
        params, model parameter metadata.
        dream, retained circular sample history.
        input, result buffers and parallel sample file.
        subset, coordinates of the current profile.
      Return:
        0 on success.
     ######################################################################*/
    
    if(!mpi || !params || !dream || !input || !subset) return -1;

    int param_idx, generation_idx, chain_idx, index, nsamples;
    bool output_enabled = dream->sample_output != SAMPLE_OUTPUT_NONE;
    if(output_enabled && !input->sample_file_open) return -1;
    int history_count = dream->retained_generations<dream->history_size 
        ? dream->retained_generations : dream->history_size;
    if(history_count < 1) return -1;
    int begin_gener = dream->retained_generations-history_count;
    int end_gener = dream->retained_generations-1;
    double ***chains = TENSOR_DBL(dream->chains);
    double **likelihood = MAT_DBL(dream->likelihood);
    double stddev[params->npar], mean[params->npar];
    if(mpi->verbose_level >= 3){
      if(Chains_STD(mpi, begin_gener, end_gener, params->npar, chains,
          dream->history_size, stddev, mean) != 0) return -1;
    }

    nsamples = history_count*mpi->total_chains;
    size_t sample_count = (size_t)nsamples;
    double *array = NULL;
    if(mpi->island_rank == 0){
      array = (double *)malloc(sample_count*sizeof(*array));
    }
    double error[params->npar][2];
    int quantile_idx = (int)(nsamples*0.15865525393145707);

    int allocation_status = mpi->island_rank != 0 || array ? 0 : -1;
    if(mpi->island_size > 1){
      MPI_Allreduce(MPI_IN_PLACE, &allocation_status, 1, MPI_INT, MPI_MIN, 
          mpi->island_comm);
    }
    if(allocation_status != 0){
      free(array);
      return -1;
    }

    long long record_offset = output_enabled
        ? SAMPLE_RECORD_OFFSET(input, subset) : 0;
    int io_status = 0;
    int sort_status = 0;

    if(mpi->island_rank == 0){
      int best_generation = begin_gener;
      int best_chain = 0;
      for(generation_idx=begin_gener; generation_idx<=end_gener;
          generation_idx++){
        int slot = generation_idx%dream->history_size;
        for(chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
          int best_slot = best_generation%dream->history_size;
          if(likelihood[slot][chain_idx] > likelihood[best_slot][best_chain]){
            best_generation = generation_idx;
            best_chain = chain_idx;
          }
        }
      }
      int best_slot = best_generation%dream->history_size;

      for(param_idx=0; param_idx<params->npar ;param_idx++){
        int model_idx = params->free_index[param_idx];
        double best_value = chains[best_slot][best_chain][param_idx];
        index = 0;
        for(generation_idx=begin_gener; generation_idx<=end_gener; 
            generation_idx++){
          for(chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
            array[index] = 
                chains[generation_idx%dream->history_size][chain_idx][param_idx];
            index++;
          }
        }

        bool periodic_phi = params->magnetic_mode == MAGNETIC_SPHERICAL
            && params->kind[model_idx] == PARAM_B2;
        if(periodic_phi){
          for(int isample=0; isample<nsamples; isample++){
            array[isample] = remainder(array[isample]-best_value, M_PI);
          }
        }
        if(SORT_VALUES(array, (size_t)nsamples) != 0){
          sort_status = -1;
          break;
        }
        double lower_quantile = array[quantile_idx];
        double upper_quantile = array[nsamples-quantile_idx-1];
        if(periodic_phi){
          error[param_idx][0] = fmax(0.0, -lower_quantile);
          error[param_idx][1] = fmax(0.0, upper_quantile);
        }else{
          error[param_idx][0] = fmax(0.0, best_value-lower_quantile);
          error[param_idx][1] = fmax(0.0, upper_quantile-best_value);
        }
      }

      if(sort_status == 0){
        for(param_idx=0; param_idx<params->nmodel; param_idx++){
          input->result_buffer[param_idx] = params->value_const[param_idx];
          input->error_buffer[param_idx*2] = 0.0;
          input->error_buffer[param_idx*2+1] = 0.0;
        }
        for(param_idx=0; param_idx<params->npar; param_idx++){
          int model_idx = params->free_index[param_idx];
          double best_value = chains[best_slot][best_chain][param_idx];
          if(params->magnetic_mode == MAGNETIC_SPHERICAL
              && params->kind[model_idx] == PARAM_B2){
            best_value = fmod(best_value, M_PI);
            if(best_value < 0.0) best_value += M_PI;
          }
          input->result_buffer[model_idx] = best_value;
          input->error_buffer[model_idx*2] = error[param_idx][0];
          input->error_buffer[model_idx*2+1] = error[param_idx][1];
        }
        if(params->magnetic_mode == MAGNETIC_SPHERICAL
            && params->nmodel > 2){
          input->result_buffer[2] = fmod(input->result_buffer[2], M_PI);
          if(input->result_buffer[2] < 0.0){
            input->result_buffer[2] += M_PI;
          }
        }
        input->result_buffer[params->nmodel] =
            likelihood[best_slot][best_chain];

        if(mpi->verbose_level >= 3){
          for(param_idx = 0; param_idx<params->npar; param_idx++){
            char name[64];
            MODEL_PARAMETER_NAME(params, params->free_index[param_idx],
                name, sizeof(name));
            snprintf(message_buffer, sizeof(message_buffer),
                "  %s statistics: best=%e, mean=%e, std=%e, "
                "68.27%% errors=[-%e, +%e]", name,
                input->result_buffer[params->free_index[param_idx]],
                mean[param_idx], stddev[param_idx],
                error[param_idx][0], error[param_idx][1]);
            LOG_WRITE(message_buffer, true, true);
          }
        }
      }
    }

    if(mpi->island_size > 1){
      MPI_Bcast(&sort_status, 1, MPI_INT, 0, mpi->island_comm);
    }
    if(sort_status != 0){
      free(array);
      return -1;
    }

    free(array);

    /* Commit the record only after every full-history sample block was
       written successfully during DREAM sampling. */
    if(output_enabled && mpi->island_rank == 0 && io_status == 0){
      int record_header[2] = {
          1, dream->sampling_generations*mpi->total_chains};
      if(MPI_File_write_at(input->sample_file, (MPI_Offset)record_offset, 
          record_header, 2, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){
        io_status = -1;
      }
    }

    if(mpi->island_size > 1){
      MPI_Bcast(&io_status, 1, MPI_INT, 0, mpi->island_comm);
    }
    
    return io_status;
}

/*--------------------------------------------------------------------------------*/
