#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "input_reader.h"
#include "logger.h"
#include "timer.h"
#include "me_solver.h"
#include "profile_io.h"
#include "profile_inversion.h"
#include "parallel_runtime.h"
#include "inversion_init.h"
#include "dream.h"
#include "sample_analysis.h"

#define TAG_WORK  1000
#define TAG_COORD 1001
#define TAG_DATA  1002
#define TAG_NOISE 1003
#define TAG_R_WORK  2000
#define TAG_R_COORD 2001
#define TAG_R_DATA  2002
#define TAG_R_ERROR 2003


#define Input_Path "./input.inv"

/*--------------------------------------------------------------------------------*/

static void REPORT_PROGRESS(const STRUCT_MPI *mpi,
    const STRUCT_PROFILE_IO *input, int completed_profiles,
    int *last_reported_row){

    /*######################################################################
      Purpose:
        Report approximate inversion progress in completed image rows.
      Input parameters:
        mpi, MPI rank and logging configuration.
        input, selected inversion geometry and total profile count.
        completed_profiles, number of cached or newly inverted profiles.
      Input/output parameters:
        last_reported_row, most recent row milestone already reported.
      Note:
        Only world rank 0 reports progress. A row count derived from the
        completed profile count is approximate when islands finish out of
        coordinate order.
    ######################################################################*/

    if(!mpi || !input || !last_reported_row || mpi->world_rank != 0
        || mpi->verbose_level < 1 || input->cache_header.nx < 1
        || input->cache_header.ny < 1 || input->counts < 1) return;

    int completed_rows = completed_profiles/input->cache_header.nx;
    if(completed_profiles >= input->counts){
      completed_rows = input->cache_header.ny;
    }

    while(*last_reported_row+5 <= completed_rows){
      *last_reported_row += 5;
      snprintf(message_buffer, sizeof(message_buffer),
          "Inversion progress: y = %d/%d.", *last_reported_row,
          input->cache_header.ny);
      LOG_WRITE(message_buffer, true, true);
    }

    if(completed_rows >= input->cache_header.ny
        && *last_reported_row < input->cache_header.ny){
      *last_reported_row = input->cache_header.ny;
      snprintf(message_buffer, sizeof(message_buffer),
          "Inversion progress: y = %d/%d.", *last_reported_row,
          input->cache_header.ny);
      LOG_WRITE(message_buffer, true, true);
    }
}

/*--------------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {

    /*######################################################################
      Purpose:
        Run the Stokes-profile inversion program with MPI island parallelism.
      Input parameters:
        argc, number of command-line arguments.
        argv, command-line arguments; argv[1] optionally selects the input file.
      Return:
        0 on success; 1 for initialization failure; 2 when one or more
        profiles fail during inversion.
     ######################################################################*/
  
    /*---------- Initialize the MPI Execution Environment ----------*/
    

    if(MPI_Init(&argc, &argv) != MPI_SUCCESS) return EXIT_FAILURE;
    MPI_Status status;

    setbuf(stderr, NULL);

    STRUCT_MPI mpi = {0};
    if(MPI_SETUP(&mpi) != 0){
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      return EXIT_FAILURE;
    }
    LOG_SET_CONSOLE(mpi.world_rank == 0);

    STRUCT_TIMER timer = {0};
    TIMER_START(&timer);

    const char *filename = (argc>=2) ? argv[1] : Input_Path;

    if(!FILE_EXIST(filename) && mpi.world_rank==0){
      LOG_ERROR(ERR_LVL_ERROR, "main", "input doesn't exist! \n");
    }

    STRUCT_PROFILE_IO input = {0};
    STRUCT_STK stokes = {0};
    STRUCT_PARA params = {0};
    STRUCT_SUBSET subset = {0};
    STRUCT_SUBSET result_subset = {0};
    STRUCT_DREAM dream = {0};
    int exit_code = 0;
    bool islands_initialized = false;

    int init_status = RDINPUT(filename, &input, &params, &dream, &stokes, &mpi);
    if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    if(init_status != 0){
      exit_code = 1;
      goto cleanup;
    }

    init_status = MPI_Init_Islands(&mpi, dream.nchains, dream.max_pairs);
    islands_initialized = true;
    if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    if(init_status != 0){
      exit_code = 1;
      goto cleanup;
    }
    if(mpi.verbose_level >= 2 && mpi.verbose_level < 4 
        && mpi.island_rank == 0 && mpi.world_rank != 0){
      char island_log[Max_Line_Length+32];
      snprintf(island_log, sizeof(island_log), "%s_island_%03d",
          input.log_path, mpi.island_id);
      if(LOG_INIT(island_log) != 0){
        LOG_ERROR(ERR_LVL_ERROR, "main", "Cannot open the island log file.");
      }
    }

    init_status = READ_WAVELENGTH(&input, &stokes, &mpi);
    if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    if(init_status != 0){
      exit_code = 1;
      goto cleanup;
    }


    init_status = INIT_INV(&input, &stokes, &params, &dream, &mpi, &subset);
    if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    if(init_status != 0){
      exit_code = 1;
      goto cleanup;
    }

    int stop_signal = -1, work_status = 1, cache_idx = 0;
    const int result_count = params.nmodel+1;
    const int error_count = params.nmodel*2;
    int local_failed_profiles = 0;
    int completed_profiles = 0, source_rank, has_work = 0;
    int last_reported_row = 0;
    double *data_ptr = NULL, *result_ptr = NULL;

    init_status = CACHE_INIT(&input, &mpi, &params, &dream, &stokes);
    if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    if(init_status == 0){
      init_status = SAMPLE_FILE_INIT(&input, &params, &dream, &mpi);
      if(MPI_SYNC_STATUS(MPI_COMM_WORLD, &init_status) != 0) init_status = -1;
    }
    if(init_status != 0){
      exit_code = 1;
      goto cleanup;
    }

    TIMER_REPORT(&timer, &mpi, "Initialization completed");
    if(mpi.world_rank == 0) LOG_WRITE("Inversion started.", 
        true, mpi.verbose_level>=1); 

    subset.processed = 0;
    if(mpi.nislands == 1){

      /* In single-island mode there is no dedicated master.  Island rank 0
         reads one task and writes its result.  With multiple ranks, it also
         broadcasts the task to the rest of the island. */
      while(1){
        if(mpi.island_rank == 0){
          while(subset.processed<input.counts 
              && input.cache[cache_idx]){
            PIXEL_ADVANCE(&input, &subset);
            cache_idx++;
            completed_profiles++;
            REPORT_PROGRESS(&mpi, &input, completed_profiles,
                &last_reported_row);
          }

          if(subset.processed>=input.counts){
            work_status = stop_signal;
          }else{
            result_subset.coord[0] = subset.coord[0];
            result_subset.coord[1] = subset.coord[1];
            READ_PROFILE(&input, &stokes, &subset);
            work_status = 1;
            cache_idx++;
          }
        }

        if(mpi.island_size > 1){
          MPI_REQUIRE_SUCCESS(MPI_Bcast(&work_status, 1, MPI_INT, 0,
              mpi.island_comm), "broadcast work status");
        }
        if(work_status < 0) break;

        if(mpi.island_size > 1){
          MPI_REQUIRE_SUCCESS(MPI_Bcast(result_subset.coord, 2, MPI_INT, 0,
              mpi.island_comm), "broadcast profile coordinates");
          MPI_REQUIRE_SUCCESS(MPI_Bcast(stokes.profile, stokes.nw*4, MPI_DOUBLE, 
              0, mpi.island_comm), "broadcast Stokes profile");
          if(stokes.noise_mode == NOISE_PER_PIXEL){
            MPI_REQUIRE_SUCCESS(MPI_Bcast(stokes.noise, stokes.nw*4, MPI_DOUBLE, 
                0, mpi.island_comm), "broadcast profile noise");
          }
        }

        int profile_status = INVERT_PROFILE(&mpi, &input, &stokes, &params,
            &dream, &result_subset);

        if(mpi.island_rank == 0){
          WRITE_RESULT(&input, &result_subset, profile_status >= 0);
          if(profile_status < 0) local_failed_profiles++;
          completed_profiles++;
          REPORT_PROGRESS(&mpi, &input, completed_profiles,
              &last_reported_row);
        }
      }
    }

    else if(mpi.is_master){

      data_ptr = stokes.profile;
      int workers_per_island = (mpi.world_size-1)/mpi.nislands;
      int extra_workers = (mpi.world_size-1)%mpi.nislands;
      int island_id = 0;

      for(cache_idx=0; island_id<mpi.nislands 
          && subset.processed<input.counts; cache_idx++){

        if(input.cache[cache_idx]){
          PIXEL_ADVANCE(&input, &subset);
          completed_profiles++;
          REPORT_PROGRESS(&mpi, &input, completed_profiles,
              &last_reported_row);
        }else{
          int leader = 1+island_id*workers_per_island 
              +(island_id<extra_workers ? island_id : extra_workers);

          MPI_REQUIRE_SUCCESS(MPI_Send(&work_status, 1, MPI_INT, leader, 
              TAG_WORK, MPI_COMM_WORLD), "send work status");
          MPI_REQUIRE_SUCCESS(MPI_Send(subset.coord, 2, MPI_INT, leader, 
              TAG_COORD, MPI_COMM_WORLD), "send profile coordinates");

          READ_PROFILE(&input, &stokes, &subset);

          MPI_REQUIRE_SUCCESS(MPI_Send(data_ptr, stokes.nw*4, MPI_DOUBLE, 
              leader, TAG_DATA, MPI_COMM_WORLD), "send Stokes profile");
          if(stokes.noise_mode == NOISE_PER_PIXEL){
            MPI_REQUIRE_SUCCESS(MPI_Send(stokes.noise, stokes.nw*4, MPI_DOUBLE,
                leader, TAG_NOISE, MPI_COMM_WORLD), "send profile noise");
          }
          island_id++;
        }
      }

      for(; island_id<mpi.nislands; island_id++){
        int leader = (1+island_id*workers_per_island 
            +(island_id<extra_workers ? island_id : extra_workers));
        MPI_REQUIRE_SUCCESS(MPI_Send(&stop_signal, 1, MPI_INT, leader, 
            TAG_WORK, MPI_COMM_WORLD), "send stop signal");
      }

      while(completed_profiles<input.counts){
        int profile_status = 0;
        MPI_RECV_EXACT(&profile_status, 1, MPI_INT, MPI_ANY_SOURCE,
            TAG_R_WORK, MPI_COMM_WORLD, &status);

        source_rank = status.MPI_SOURCE;
        if(!MPI_IS_ISLAND_LEADER(&mpi, source_rank)){
          LOG_ERROR(ERR_LVL_ERROR, "main",
              "received a result from a non-leader MPI rank.\n");
        }
        MPI_RECV_EXACT(result_subset.coord, 2, MPI_INT, source_rank,
            TAG_R_COORD, MPI_COMM_WORLD, &status);

        result_ptr = input.result_buffer;
          
        MPI_RECV_EXACT(result_ptr, result_count, MPI_DOUBLE, source_rank,
            TAG_R_DATA, MPI_COMM_WORLD, &status);

        MPI_RECV_EXACT(input.error_buffer, error_count, MPI_DOUBLE,
            source_rank, TAG_R_ERROR, MPI_COMM_WORLD, &status);

        if(profile_status >= 0
            && (!INVERSION_RESULT_FINITE(input.result_buffer, result_count)
            || !INVERSION_RESULT_FINITE(input.error_buffer, error_count))){
          profile_status = -1;
          memset(input.result_buffer, 0,
              (size_t)result_count*sizeof(*input.result_buffer));
          memset(input.error_buffer, 0,
              (size_t)error_count*sizeof(*input.error_buffer));
        }
        WRITE_RESULT(&input, &result_subset, profile_status >= 0);
        if(profile_status < 0) local_failed_profiles++;

        completed_profiles++;
        REPORT_PROGRESS(&mpi, &input, completed_profiles,
            &last_reported_row);
        has_work = 0;
        while(subset.processed < input.counts){
          if(input.cache[cache_idx]){
            cache_idx++;
            PIXEL_ADVANCE(&input, &subset);
            completed_profiles++;
            REPORT_PROGRESS(&mpi, &input, completed_profiles,
                &last_reported_row);
          }else{
            has_work = 1;
            break;
          }
        }

        if(has_work){
          cache_idx++;
          MPI_REQUIRE_SUCCESS(MPI_Send(&work_status, 1, MPI_INT, source_rank,
              TAG_WORK, MPI_COMM_WORLD), "send work status");
          MPI_REQUIRE_SUCCESS(MPI_Send(subset.coord, 2, MPI_INT, source_rank,
              TAG_COORD, MPI_COMM_WORLD), "send profile coordinates");

     
          READ_PROFILE(&input, &stokes, &subset);
            
          MPI_REQUIRE_SUCCESS(MPI_Send(data_ptr, stokes.nw*4, MPI_DOUBLE,
              source_rank, TAG_DATA, MPI_COMM_WORLD), "send Stokes profile");
          if(stokes.noise_mode == NOISE_PER_PIXEL){
            MPI_REQUIRE_SUCCESS(MPI_Send(stokes.noise, stokes.nw*4, MPI_DOUBLE,
                source_rank, TAG_NOISE, MPI_COMM_WORLD),
                "send profile noise");
          }
        }else{
          MPI_REQUIRE_SUCCESS(MPI_Send(&stop_signal, 1, MPI_INT, source_rank,
              TAG_WORK, MPI_COMM_WORLD), "send stop signal");
        }
      }
  
    }else{
      while(1){
        if(mpi.island_rank == 0){
          MPI_RECV_EXACT(&work_status, 1, MPI_INT, 0, TAG_WORK,
              MPI_COMM_WORLD, &status);
        }
        MPI_REQUIRE_SUCCESS(MPI_Bcast(&work_status, 1, MPI_INT, 0,
            mpi.island_comm), "broadcast work status");
        if(work_status < 0) break;

        if(mpi.island_rank == 0){
          MPI_RECV_EXACT(subset.coord, 2, MPI_INT, 0, TAG_COORD,
              MPI_COMM_WORLD, &status);
          MPI_RECV_EXACT(stokes.profile, stokes.nw*4, MPI_DOUBLE, 0,
              TAG_DATA, MPI_COMM_WORLD, &status);
          if(stokes.noise_mode == NOISE_PER_PIXEL){
            MPI_RECV_EXACT(stokes.noise, stokes.nw*4, MPI_DOUBLE, 0,
                TAG_NOISE, MPI_COMM_WORLD, &status);
          }
        }
        MPI_REQUIRE_SUCCESS(MPI_Bcast(subset.coord, 2, MPI_INT, 0,
            mpi.island_comm), "broadcast profile coordinates");
        MPI_REQUIRE_SUCCESS(MPI_Bcast(stokes.profile, stokes.nw*4, MPI_DOUBLE, 
            0, mpi.island_comm), "broadcast Stokes profile");
        if(stokes.noise_mode == NOISE_PER_PIXEL){
          MPI_REQUIRE_SUCCESS(MPI_Bcast(stokes.noise, stokes.nw*4, MPI_DOUBLE, 
              0, mpi.island_comm), "broadcast profile noise");
        }

        int profile_status = INVERT_PROFILE(&mpi, &input, &stokes, &params,
            &dream, &subset);
        if(mpi.island_rank == 0){
          MPI_REQUIRE_SUCCESS(MPI_Send(&profile_status, 1, MPI_INT, 0, 
              TAG_R_WORK, MPI_COMM_WORLD), "send profile status");
          MPI_REQUIRE_SUCCESS(MPI_Send(subset.coord, 2, MPI_INT, 0, 
              TAG_R_COORD, MPI_COMM_WORLD), "send result coordinates");
          MPI_REQUIRE_SUCCESS(MPI_Send(input.result_buffer, result_count, 
              MPI_DOUBLE, 0, TAG_R_DATA, MPI_COMM_WORLD), 
              "send inversion result");
          MPI_REQUIRE_SUCCESS(MPI_Send(input.error_buffer, error_count, 
              MPI_DOUBLE, 0, TAG_R_ERROR, MPI_COMM_WORLD), 
              "send parameter errors");
        }
      }
    }
    
    int total_failed_profiles = 0;
    MPI_REQUIRE_SUCCESS(MPI_Allreduce(&local_failed_profiles, 
        &total_failed_profiles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD),
        "reduce failed-profile count");
    if(total_failed_profiles > 0){
      exit_code = 2;
      if(mpi.world_rank == 0){
        char message[128];
        snprintf(message, sizeof(message),
            "%d profile inversion(s) failed; zero results were written "
            "without marking them complete in the cache.\n",
            total_failed_profiles);
        LOG_ERROR(ERR_LVL_WARNING, "main", message);
      }
    }

cleanup:
    if(SAMPLE_FILE_CLOSE(&input) != 0 && mpi.world_rank == 0){
      LOG_ERROR(ERR_LVL_WARNING, "main", "failed to close the sample file.\n");
      exit_code = 1;
    }

    FREE_TENSOR(dream.chains);
    FREE_MATRIX(dream.likelihood);
    Free_Dream(&dream);

    if(params.bounds){
      free(params.bounds[0]);
      free(params.bounds);
    }
    free(params.proposal_scale);
    free(params.type);
    free(params.narrower_guess);
    free(params.step);
    free(params.free_index);
    free(params.lines);
    free(params.regions);

    free(stokes.wavelength);
    free(stokes.profile);
    free(stokes.synthetic);
    free(stokes.inv_noise);
    free(stokes.noise);
    free(stokes.faddeeva.exp_a2n2);

    free(input.result_buffer);
    free(input.error_buffer);
    free(input.cache);
    TIMER_STOP(&timer, &mpi);
    if(CLOSE_FILES() != 0) exit_code = 1;

#ifdef USE_GSL
    if(mpi.rank_rng && mpi.rank_rng->gsl) gsl_rng_free(mpi.rank_rng->gsl);
#endif
    free(mpi.rank_rng);

    if(islands_initialized && mpi.island_comm != MPI_COMM_NULL 
        && mpi.island_comm != MPI_COMM_WORLD){
      if(MPI_Comm_free(&mpi.island_comm) != MPI_SUCCESS) exit_code = 1;
    }
    if(MPI_Finalize() != MPI_SUCCESS) exit_code = 1;

    return exit_code;
}

/*--------------------------------------------------------------------------------*/
