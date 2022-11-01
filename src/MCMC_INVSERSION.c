#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "GEMC.h"
#include "READ_INPUT_MCMC.h"
#include "TIME_PRINT.h"
#include "MPI_INIT.h"
#include "STR.h"
#include "FORWARD_ME.h"
#include "INIT_BOUNDS.h"
#include "DREAM.h"
#include "MATH_TOOL.h"
#include "SAMPLE_ANALY.h"


#define Path_Input "./input.dat"

/*--------------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
  
    /*---------- Initialize the MPI Execution Environment ----------*/
    
    MPI_Init(&argc,&argv);
    
    STRUCT_MPI *Mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &Mpi_pid);
    MPI_Comm_size(MPI_COMM_WORLD, &Mpi_num_procs);
    Mpi = (STRUCT_MPI *)malloc(sizeof(STRUCT_MPI));
    Mpi->Rank = Mpi_pid;
    Mpi->Num_procs = Mpi_num_procs;
    
    /*---------- Begin to calculate the Running Time ----------*/
    
    if(Mpi->Rank == 0) Time_Print();
    char filename[Max_Line_Length]; 
    if(argc > 1){
      String_Copy(filename, argv[1], strlen(argv[1]), true);
    }else{
      String_Copy(filename, Path_Input, strlen(Path_Input), true);
    }

    STRUCT_INPUT *Input = (STRUCT_INPUT *)malloc(sizeof(STRUCT_INPUT));
    STRUCT_OBSERVATION *Observation = (STRUCT_OBSERVATION *)\
        malloc(sizeof(STRUCT_OBSERVATION));
    STRUCT_PAR *Str_Par = (STRUCT_PAR *)malloc(sizeof(STRUCT_PAR));
    STRUCT_ATOM *Atom = (STRUCT_ATOM *)malloc(sizeof(STRUCT_ATOM));
    STRUCT_GEMC *Gemc = (STRUCT_GEMC *)malloc(sizeof(STRUCT_GEMC));
    STRUCT_DREAM *Dream = (STRUCT_DREAM *)malloc(sizeof(STRUCT_DREAM));

    //int (*Forward)(STRUCT_ATOM *, double *, STRUCT_OBSERVATION *) \
        = Forward_ME;

    int Chain_bestfit;

    RDINPUT(filename, Input, Mpi, Atom, Observation);

    MPI_INIT(Input, Mpi);

    Init_Bounds(Str_Par, Observation, Atom);

    if(Input->GEMC_Simulation){
      Chain_Init_Gemc(Mpi, Gemc, Atom, Observation, Str_Par);

      Chain_bestfit = GEMC(Mpi, Gemc, Atom, Observation, Str_Par);

      GEMC2DREAM(Input, Mpi, Gemc, Dream, Str_Par, Observation, Atom);
      FREE_MATRIX(Gemc->Chains, 0, 0, enum_dbl);  
      FREE_VECTOR(Gemc->Likelihood, 0, enum_dbl);  
      
      if(Mpi->Rank == 0) Time_Print();
    }else{
	    Chain_Init_Dream(Input, Mpi, Dream, Atom, Observation, \
          Str_Par);
    }

    Chain_bestfit = DREAM(Input, Mpi, Dream, Atom, Observation, \
        Str_Par);

    SAMPLE_ANALY(Input, Mpi, Str_Par, Dream, Chain_bestfit);

    if(Mpi->Rank == 0) Time_Print();

    FREE_MATRIX(Dream->Likelihood, 0, 0, enum_dbl);  
    FREE_TENSOR_DBL(Dream->Chains, 0, 0, 0);
    free(Dream);  
    free(Gemc);
    free(Input);
    FREE_MATRIX(Observation->Profile_Data, 0, 0, enum_dbl);  	
    FREE_MATRIX(Observation->Sigma, 0, 0, enum_dbl);  	
    FREE_MATRIX(Observation->Stokes, 0, 0, enum_dbl);  	
    FREE_VECTOR(Observation->Lambda, 0, enum_dbl);
    free(Observation);
    FREE_MATRIX(Str_Par->Par_Bounds, 0, 0, enum_dbl);  
    FREE_VECTOR(Str_Par->narrower_guess, 0, enum_dbl);
    FREE_VECTOR(Str_Par->type, 0, enum_dbl);
    FREE_VECTOR(Str_Par->Distribution, 0, enum_dbl);	
    free(Str_Par);
    free(Atom);
    free(Mpi->idum);
    free(Mpi);
    

    MPI_Finalize();
    return 0;

    /*----------                       END                      ----------*/
    
}

/*--------------------------------------------------------------------------------*/
