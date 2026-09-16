
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>
#include <mpi.h>
#include "random.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_MPI{
   
    int world_rank, world_size;
    int island_id, island_rank, island_size, nislands;
    int requested_island_size;
    int verbose_level;
    bool is_master;
    MPI_Comm world_comm, island_comm;
    int chain_begin, chain_end, chains_per_rank, total_chains;

    STRUCT_RNGState *rank_rng;

}STRUCT_MPI;

/*--------------------------------------------------------------------------------*/

extern int MPI_SETUP(STRUCT_MPI *mpi);

extern int MPI_Init_Islands(STRUCT_MPI *mpi, int nchains, int max_pairs);

extern int MPI_SYNC_STATUS(MPI_Comm comm, int *status);

extern void MPI_REQUIRE_SUCCESS(int mpi_status, const char *operation);

extern void MPI_RECV_EXACT(void *buffer, int count, MPI_Datatype datatype,
    int source, int tag, MPI_Comm comm, MPI_Status *status);

extern bool MPI_IS_ISLAND_LEADER(const STRUCT_MPI *mpi, int world_rank);

/*--------------------------------------------------------------------------------*/
