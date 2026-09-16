
#include "timer.h"
#include "logger.h"
#include <stdio.h>

/*--------------------------------------------------------------------------------*/

static double TIMER_NOW(void){
    /*######################################################################
      Purpose:
        Read a monotonic wall-clock timestamp for the active runtime.
      Return:
        Current wall-clock value in seconds.
      Note:
        The program uses MPI_Wtime because MPI is always enabled.
    ######################################################################*/

    return MPI_Wtime();
}

/*--------------------------------------------------------------------------------*/

void TIMER_START(STRUCT_TIMER *timer){

    /*######################################################################
      Purpose:
        Start or restart a wall-clock timer.
      Input parameters:
        timer, timer state to initialize.
      Output parameters:
        timer, populated with the current time and marked as running.
    ######################################################################*/

    if(!timer) return;
    timer->start_time = TIMER_NOW();
    timer->running = true;
}

/*--------------------------------------------------------------------------------*/

double TIMER_REPORT(const STRUCT_TIMER *timer, STRUCT_MPI *mpi, 
    const char *label){

    /*######################################################################
      Purpose:
        Report elapsed wall time without stopping the timer.
      Input parameters:
        timer, a timer previously initialized by TIMER_START.
        mpi, MPI communicator, rank, and logging configuration.
        label, description printed before the elapsed time; NULL or an empty
          string selects the default label.
      Return:
        Maximum elapsed wall time in seconds; -1 for invalid state.
      Note:
        Every MPI world rank must call this function collectively.
    ######################################################################*/

    if(!timer || !mpi || !timer->running) return -1.0;

    double elapsed = TIMER_NOW()-timer->start_time;
    if(elapsed < 0.0) elapsed = 0.0;
    double global_elapsed = 0.0;
    MPI_Allreduce(&elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX,
        mpi->world_comm);
    elapsed = global_elapsed;
    if(mpi->world_rank == 0 && mpi->verbose_level >= 1){
      const char *text = label && label[0] ? label : "Elapsed time";
      long hours = (long)(elapsed/3600.0);
      long minutes = (long)((elapsed-(double)hours*3600.0)/60.0);
      double seconds = elapsed-(double)hours*3600.0
          -(double)minutes*60.0;
      if(hours > 0){
        snprintf(message_buffer, sizeof(message_buffer),
            "%s: %ld h %ld min %.2f s", text, hours, minutes, seconds);
      }else if(minutes > 0){
        snprintf(message_buffer, sizeof(message_buffer),
            "%s: %ld min %.2f s", text, minutes, seconds);
      }else{
        snprintf(message_buffer, sizeof(message_buffer),
            "%s: %.2f s", text, seconds);
      }
      LOG_WRITE(message_buffer, true, true);
    }

    return elapsed;
}

/*--------------------------------------------------------------------------------*/

double TIMER_STOP(STRUCT_TIMER *timer, STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
        Report final elapsed wall time and stop the timer.
      Input parameters:
        timer, a timer previously initialized by TIMER_START.
        mpi, MPI communicator, rank, and logging configuration.
      Output parameters:
        timer, marked as no longer running after a valid call.
      Return:
        Maximum elapsed wall time in seconds; -1 for invalid state.
      Note:
        Every MPI world rank must call this function collectively.
    ######################################################################*/

    double elapsed = TIMER_REPORT(timer, mpi, "Execution time");
    if(elapsed >= 0.0) timer->running = false;

    return elapsed;

}

/*--------------------------------------------------------------------------------*/
