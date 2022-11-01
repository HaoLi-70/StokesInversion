
#ifndef MPI_INIT_h
#define MPI_INIT_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include "MPI_CONTROL.h"
#include "READ_INPUT_MCMC.h"

/*--------------------------------------------------------------------------------*/

extern int MPI_INIT(STRUCT_INPUT *Input, STRUCT_MPI *Mpi);

/*--------------------------------------------------------------------------------*/

#endif /* MPI_INIT_h */
