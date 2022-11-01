
#ifndef SAMPLE_ANALY_h
#define SAMPLE_ANALY_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>

#include "ALLOCATION.h"
#include "MPI_CONTROL.h"
#include "READ_INPUT_MCMC.h"
#include "SORT.h"
#include "DREAM.h"

/*--------------------------------------------------------------------------------*/

extern int SAMPLE_ANALY(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_PAR *Str_Par, STRUCT_DREAM *Dream, int Chain_bestfit);

/*--------------------------------------------------------------------------------*/

#endif /* SAMPLE_ANALY_h */
