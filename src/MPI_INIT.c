
#include "MPI_INIT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:
        8 Sept. 2021.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern int MPI_INIT(STRUCT_INPUT *Input, STRUCT_MPI *Mpi){

    /*######################################################################
      Purpose:
        initialize the mpi structure.
      Record of revisions:
        30 Otc., 2022
      Input parameters:
        Mpi, a structure save the information for mpi.
      Output parameters:
        Mpi, the structure save the information for mpi.
     ######################################################################*/
    
    int i;
    Mpi->num_gemc = Input->Gemc_Nchains/Mpi->Num_procs;

    if ((Mpi->num_gemc+1)*Mpi->Num_procs < 250) Mpi->num_gemc++;

    Mpi->indxb_gemc = Mpi->num_gemc*Mpi->Rank;
    Mpi->indxe_gemc = Mpi->indxb_gemc+Mpi->num_gemc-1;
    Mpi->tot_gemc = Mpi->num_gemc*Mpi->Num_procs;

    Mpi->num_dream = Input->Dream_Nchains/Mpi->Num_procs;

    if ((Mpi->num_dream+1)*Mpi->Num_procs < 60) Mpi->num_dream++;

    Mpi->indxb_dream = Mpi->num_dream*Mpi->Rank;
    Mpi->indxe_dream = Mpi->indxb_dream+Mpi->num_dream-1;
    Mpi->tot_dream = Mpi->num_dream*Mpi->Num_procs;

    Mpi->idum=(long *)malloc(sizeof(long)*Mpi->Num_procs);

    if (Mpi->Rank == 0){
      srand((unsigned int)time(NULL));
      srand(rand());  
      for(i=0; i<Mpi->Num_procs; i++){  
        Mpi->idum[i] = -rand();
      }
    }
    MPI_Bcast(Mpi->idum, Mpi->Num_procs, MPI_LONG, 0, MPI_COMM_WORLD);

    return 0;
}

/*--------------------------------------------------------------------------------*/

