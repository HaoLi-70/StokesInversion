
#ifndef MPI_CONTROL_h
#define MPI_CONTROL_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

/*--------------------------------------------------------------------------------*/

extern int Mpi_pid;
extern int Mpi_num_procs;

typedef struct Struct_MPI{
   
    int Rank, Num_procs;
    bool Debug;

    int indxb_gemc, indxe_gemc, num_gemc, tot_gemc;
    int indxb_dream, indxe_dream, num_dream, tot_dream;

    long *idum;
}STRUCT_MPI;

extern int CONTROL(void);

extern int ABORTED(void);

/*--------------------------------------------------------------------------------*/

#endif /* MPI_CONTROL_h */
