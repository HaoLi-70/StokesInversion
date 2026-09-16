
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>

#include "parallel_runtime.h"

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Timer{
    double start_time;
    bool running;
}STRUCT_TIMER;

/*--------------------------------------------------------------------------------*/

extern void TIMER_START(STRUCT_TIMER *timer);

extern double TIMER_REPORT(const STRUCT_TIMER *timer, STRUCT_MPI *mpi, 
    const char *label);

extern double TIMER_STOP(STRUCT_TIMER *timer, STRUCT_MPI *mpi);

/*--------------------------------------------------------------------------------*/
