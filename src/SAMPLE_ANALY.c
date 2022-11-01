
#include "SAMPLE_ANALY.h"

/*--------------------------------------------------------------------------------*/

extern int SAMPLE_ANALY(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_PAR *Str_Par, STRUCT_DREAM *Dream, int Chain_bestfit){
    
    /*######################################################################
      Purpose:
        Calculate the standard deviations for each parameters
      Record of revisions::
        15 Jun 2018.
      Input parameters:
        Num_Chain, the total number of chains.
        Begin_Gener, the fist generation (include) used in the calculation.
        End_Gener, the last generation (not include) used in the calculation.
        Num_Par, the total number of parameters.
        Chains[][][], model parameters in all chains.
      Output parameters:
        STD[], return the standard deviation of each parameter.
        Mean[], return the mean value of each patameter.
      Method:
        Sigma=sqrt((sum(Xi^2)-N*Xbar*Xbar)/(N-1)).
     ######################################################################*/
    
    int Indx_Par, Indx_Gener, Indx_Chain, Indx, Nsamples;
    double STD[Str_Par->Num_Par], Mean[Str_Par->Num_Par];
    Chains_STD(Mpi, 0, Input->Num_Gener-1, Str_Par->Num_Par, \
        Dream->Chains, STD, Mean);

    Nsamples = Input->Num_Gener*Mpi->tot_dream;
    double *array = (double *)malloc(sizeof(double)*Nsamples);
    double sigma[Str_Par->Num_Par][2];
    int sigmaindx = (int)(Nsamples*0.12);

    FILE *fa;
    if(Input->Output_Samples && Mpi->Rank == 0) {
      fa = fopen("./Samples.bin", "wb+");
      fwrite(&Nsamples, sizeof(int), 1, fa);
      fwrite(&Str_Par->Num_Par, sizeof(int), 1, fa);
    }

    if(Mpi->Rank == 0 && Input->Output_Samples && Input->Sigma){
      for(Indx_Par=0; Indx_Par<Str_Par->Num_Par ;Indx_Par++){
        Indx = 0;
        for(Indx_Gener=0; Indx_Gener<Input->Num_Gener; Indx_Gener++){
          for(Indx_Chain=0; Indx_Chain<Mpi->tot_dream; Indx_Chain++){
            array[Indx] = Dream->Chains[Indx_Gener][Indx_Chain][Indx_Par];
            Indx++;
          }
        }
        if(Input->Output_Samples) fwrite(array, sizeof(double), Nsamples, fa);
        if(Input->Sigma) quick_sort(0, Nsamples-1, array);
        sigma[Indx_Par][0] = array[sigmaindx];
        sigma[Indx_Par][1] = array[Nsamples-sigmaindx-1];
        //if(Input->Output_Samples) fwrite(array, sizeof(double), Nsamples, fa);
      }
    } 

    if(Mpi->Rank == 0){
      if(Input->Sigma){
        for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
          fprintf(stderr, "%d bestfit %e, mean %e, std %e, 1 sigma %e %e \n", \
              Indx_Par, Dream->Chains[Input->Num_Gener-1][Chain_bestfit][Indx_Par], \
              Mean[Indx_Par], STD[Indx_Par], \
              sigma[Indx_Par][0], sigma[Indx_Par][1]);
        }
      }else{
        for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
          fprintf(stderr, "%d bestfit %e, mean %e, std %e \n", \
              Indx_Par, Dream->Chains[Input->Num_Gener-1][Chain_bestfit][Indx_Par], \
              Mean[Indx_Par], STD[Indx_Par]);
        }
      }
    }

    if(Input->Output_Samples && Mpi->Rank == 0) fclose(fa);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/
