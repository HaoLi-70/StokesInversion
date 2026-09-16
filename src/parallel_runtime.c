
#include "parallel_runtime.h"
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        MPI world initialization, balanced island construction, per-island
        DREAM chain distribution, and rank-local random-state setup.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

int MPI_SYNC_STATUS(MPI_Comm comm, int *status){

    /*######################################################################
      Purpose:
        Synchronize a status value across an MPI communicator.
      Input parameters:
        comm, communicator whose ranks participate in the synchronization.
        status, local status value; negative errors take precedence.
      Output parameters:
        status, synchronized status shared by every participating rank.
      Return:
        0 on success; -1 for invalid input or an MPI failure.
    ######################################################################*/

    if(!status || comm == MPI_COMM_NULL) return -1;
    int min_status = 0, max_status = 0;
    int min_result = MPI_Allreduce(status, &min_status, 1, MPI_INT, MPI_MIN,
        comm);
    int max_result = MPI_Allreduce(status, &max_status, 1, MPI_INT, MPI_MAX,
        comm);
    if(min_result != MPI_SUCCESS || max_result != MPI_SUCCESS) return -1;
    *status = min_status < 0 ? min_status : max_status;

    return 0;
}

/*--------------------------------------------------------------------------------*/

void MPI_REQUIRE_SUCCESS(int mpi_status, const char *operation){

    /*######################################################################
      Purpose:
        Abort the MPI job after a failed required communication operation.
      Input parameters:
        mpi_status, MPI return code to validate.
        operation, short description used in the error message.
    ######################################################################*/

    if(mpi_status == MPI_SUCCESS) return;
    char message[256];
    snprintf(message, sizeof(message), "MPI operation failed: %s.\n",
        operation ? operation : "unknown");
    LOG_ERROR(ERR_LVL_ERROR, "MPI_REQUIRE_SUCCESS", message);
}

/*--------------------------------------------------------------------------------*/

void MPI_RECV_EXACT(void *buffer, int count, MPI_Datatype datatype,
    int source, int tag, MPI_Comm comm, MPI_Status *status){

    /*######################################################################
      Purpose:
        Receive an MPI message and require its element count to match exactly.
      Input parameters:
        buffer, destination message buffer.
        count, expected number of elements.
        datatype, MPI datatype of each element.
        source, source rank or MPI_ANY_SOURCE.
        tag, required message tag.
        comm, communicator used for the receive.
      Output parameters:
        buffer, received message data.
        status, MPI receive status.
    ######################################################################*/

    MPI_REQUIRE_SUCCESS(
        MPI_Recv(buffer, count, datatype, source, tag, comm, status),
        "MPI_Recv");
    int received = 0;
    MPI_REQUIRE_SUCCESS(MPI_Get_count(status, datatype, &received),
        "MPI_Get_count");
    if(received != count){
      LOG_ERROR(ERR_LVL_ERROR, "MPI_RECV_EXACT",
          "MPI message has an invalid size.\n");
    }
}

/*--------------------------------------------------------------------------------*/

bool MPI_IS_ISLAND_LEADER(const STRUCT_MPI *mpi, int world_rank){

    /*######################################################################
      Purpose:
        Determine whether a world rank is a configured island leader.
      Input parameters:
        mpi, balanced island layout.
        world_rank, world rank to inspect.
      Return:
        true for an island leader; otherwise false.
    ######################################################################*/

    if(!mpi || mpi->nislands < 2 || world_rank <= 0
        || world_rank >= mpi->world_size) return false;
    int workers = mpi->world_size-1;
    int workers_per_island = workers/mpi->nislands;
    int extra_workers = workers%mpi->nislands;
    for(int island=0; island<mpi->nislands; island++){
      int leader = 1+island*workers_per_island
          +(island < extra_workers ? island : extra_workers);
      if(world_rank == leader) return true;
    }

    return false;
}

/*--------------------------------------------------------------------------------*/

int MPI_SETUP(STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
          Initialize the basic MPI communicator, rank, and process count.
        Input parameters:
        mpi, MPI state to initialize.
        Output parameters:
          mpi, populated with the active communicator, rank, and process count.
    ######################################################################*/
  
    if(!mpi) return -1;

    mpi->world_comm = MPI_COMM_WORLD;
    mpi->island_comm = mpi->world_comm;
    if(MPI_Comm_rank(mpi->world_comm, &mpi->world_rank) != MPI_SUCCESS
        || MPI_Comm_size(mpi->world_comm, &mpi->world_size) != MPI_SUCCESS
        || mpi->world_size < 1) return -1;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int MPI_Init_Islands(STRUCT_MPI *mpi, int nchains, int max_pairs){

    /*######################################################################
      Purpose:
        Divide worker ranks into balanced MPI islands and assign DREAM chains.
      Input parameters:
        mpi, basic MPI state initialized by MPI_SETUP.
        nchains, requested minimum number of DREAM chains in each island.
        max_pairs, maximum number of chain pairs used by DREAM.
      Output parameters:
        mpi, populated island topology, chain range, and RNG states.
      Return:
        0 on success.
    ######################################################################*/
    
    if(!mpi || mpi->requested_island_size < 1 || nchains <= 10 || max_pairs < 1
        || mpi->world_size < 1) return -1;

    int worker_count = mpi->world_size - 1;
    int island_size = mpi->requested_island_size;

    if(mpi->world_size == 1) worker_count = 1;
    if(worker_count < 1) island_size = 1;

    if(island_size > worker_count && worker_count > 0){
      island_size = worker_count;
    }
    
    mpi->nislands = worker_count/island_size;

    /* A dedicated master is useful only when it dispatches work to multiple
       islands.  For a single island, rank 0 joins the island and performs I/O. */
    if(mpi->nislands == 1){
      mpi->is_master = false;
      worker_count = mpi->world_size;
      island_size = worker_count;
    }else{
      mpi->is_master = (mpi->world_rank == 0);
    }

    /* Distribute workers as evenly as possible among the islands.  The first
       extra_workers islands receive one additional worker. */
    int workers_per_island = worker_count/mpi->nislands;
    int extra_workers = worker_count%mpi->nislands;

    if(mpi->nislands == 1){
      mpi->island_id = 0;
    }else if(mpi->is_master){
      mpi->island_id = -1;
    }else{
      int worker_rank = mpi->world_rank-1;
      int workers_in_larger_islands = 
        (workers_per_island+1)*extra_workers;

      if(worker_rank < workers_in_larger_islands){
        mpi->island_id = worker_rank/(workers_per_island+1);
      }else{
        mpi->island_id = extra_workers 
          +(worker_rank-workers_in_larger_islands)/workers_per_island;
      }
    }

    int mpi_status = MPI_Comm_split(mpi->world_comm,
      mpi->is_master ? MPI_UNDEFINED : mpi->island_id, mpi->world_rank, 
      &mpi->island_comm);
    if(!mpi->is_master){
      if(mpi_status == MPI_SUCCESS){
        mpi_status = MPI_Comm_rank(mpi->island_comm, &mpi->island_rank);
      }
      if(mpi_status == MPI_SUCCESS){
        mpi_status = MPI_Comm_size(mpi->island_comm, &mpi->island_size);
      }
    }

    int local_ok = mpi_status == MPI_SUCCESS;
    int all_ok = 0;
    if(MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_LAND,
        mpi->world_comm) != MPI_SUCCESS || !all_ok) return -1;

    if(mpi->is_master){
      mpi->island_rank = -1;
      mpi->island_size = 0;
      mpi->rank_rng = NULL;
      mpi->chains_per_rank = 0;
      mpi->total_chains = 0;
      return 0;
    }

    mpi->chains_per_rank = (nchains+mpi->island_size-1)
      /mpi->island_size;

    mpi->chain_begin = mpi->chains_per_rank*mpi->island_rank;
    mpi->chain_end = mpi->chain_begin+mpi->chains_per_rank-1;
    mpi->total_chains = mpi->chains_per_rank*mpi->island_size;

    if(max_pairs < 1 || 2LL*max_pairs > (long long)mpi->total_chains-1){
      if(mpi->island_id == 0 && mpi->island_rank == 0){
        LOG_ERROR(ERR_LVL_WARNING, "MPI_Init_Islands", 
            "max_pair requires at least 2*max_pair+1 chains per island.\n");
      }
      return -1;
    }

    mpi->rank_rng = (STRUCT_RNGState *)calloc(1, sizeof(*mpi->rank_rng));
    if(!mpi->rank_rng) return -1;
    if(RNG_INIT(mpi->rank_rng, mpi->world_rank) != 0){
      free(mpi->rank_rng);
      mpi->rank_rng = NULL;
      return -1;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/
