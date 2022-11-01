
#include "GEMC.h"

/*--------------------------------------------------------------------------------*/

extern int GEMC(STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, STRUCT_ATOM *Atom, \
    STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par){

    /*######################################################################
      Purpose:
        GEMC simulation.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Gemc, a structure saved the GEMC samples and likelihood.
        Atom, a structure the atomic information
        Observation, a structure saved observed profile.
        Str_Par, a structure saved the parameter bounds.
      Output parameters:
        Gemc, a structure saved the GEMC samples and likelihood.
      Reference:
          Tregloan-Reed 2013 MNRAS
     ######################################################################*/

    if(Mpi->Rank ==0) fprintf(stderr, "----------GEMC SAMPLING ----------\n");
    Observation->weights_flag = false;

    int Indx_Chain=0, Indx_Par, Chain_bestfit=0, ITERATION=0, REIT=0, i;
    int count_sample=0, count_accept=0, tot_sample=0, tot_accept=0;
    double Bestliklihood, likelihood;
    double TMP_LIKELIHOOD, Varance = 0, Varance_new = 0.0;
    double Delta=0;

    double *Bestsample = (double *)VECTOR(0, Str_Par->Num_Par, enum_dbl, false);
    double *Sample = (double *)VECTOR(0, Str_Par->Num_Par, enum_dbl, false);

    Chain_bestfit = BESTFIT_CHAIN(Mpi, Gemc->Likelihood, Mpi->tot_gemc, \
        Gemc->Chains, Str_Par->Num_Par, &Bestliklihood, Bestsample);

    likelihood = Bestliklihood;

    Varance = Chain_Varance(Mpi, Gemc->Likelihood, Mpi->tot_gemc);

    for(REIT=0; REIT<=10; REIT++){
      for(ITERATION=0; ITERATION<=70; ITERATION++){
        for(Indx_Chain=Mpi->indxb_gemc; Indx_Chain<=Mpi->indxe_gemc; Indx_Chain++){
          if(Indx_Chain!=Chain_bestfit){
            for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){         
              Delta = Bestsample[Indx_Par] - Gemc->Chains[Indx_Chain][Indx_Par];   
              
              if(fabs(Delta)>Str_Par->Distribution[Indx_Par]*0.1){
                Sample[Indx_Par] = Gemc->Chains[Indx_Chain][Indx_Par] \
                    +Delta*Ran1(Mpi->idum+Mpi->Rank)*2.;
              }else{
                Sample[Indx_Par] = Gemc->Chains[Indx_Chain][Indx_Par] \
                    +Str_Par->Distribution[Indx_Par] \
                    *(Ran1(Mpi->idum+Mpi->Rank)-0.5)*0.01; 
              }  
            }
            Bounds_Enforce(Sample, Str_Par);
          
            count_sample++;
            TMP_LIKELIHOOD = Likelihood_Log(Atom, Sample, Observation);

            if(log(Ran1(Mpi->idum+Mpi->Rank)) \
                <(TMP_LIKELIHOOD-Gemc->Likelihood[Indx_Chain])){
              Gemc->Likelihood[Indx_Chain] = TMP_LIKELIHOOD;
              for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
                Gemc->Chains[Indx_Chain][Indx_Par] = Sample[Indx_Par];
              }
              count_accept++;               
            }

          }else{
            for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){         
              Sample[Indx_Par] = Gemc->Chains[Indx_Chain][Indx_Par] \
                  +Str_Par->Distribution[Indx_Par] \
                  *(Ran1(Mpi->idum+Mpi->Rank)-0.5)*0.01;      
            }
            Bounds_Enforce(Sample, Str_Par);
            count_sample++;
            TMP_LIKELIHOOD = Likelihood_Log(Atom, Sample, Observation);

            if(log(Ran1(Mpi->idum+Mpi->Rank)) \
                <TMP_LIKELIHOOD-Gemc->Likelihood[Indx_Chain]){
              Gemc->Likelihood[Indx_Chain] = TMP_LIKELIHOOD;
              for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
                Gemc->Chains[Indx_Chain][Indx_Par] = Sample[Indx_Par];
              }
              count_accept++;               
            }
          }
        }

        for(i=0; i<Mpi->Num_procs; i++){
          MPI_Bcast(Gemc->Chains[i*Mpi->num_gemc], Str_Par->Num_Par*Mpi->num_gemc, \
              MPI_DOUBLE, i, MPI_COMM_WORLD);
          MPI_Bcast(Gemc->Likelihood+i*Mpi->num_gemc, Mpi->num_gemc, MPI_DOUBLE, i, \
              MPI_COMM_WORLD );
        }

        Chain_bestfit = BESTFIT_CHAIN(Mpi, Gemc->Likelihood, Mpi->tot_gemc, \
            Gemc->Chains, Str_Par->Num_Par, &Bestliklihood, Bestsample);

        Varance_new = Chain_Varance(Mpi, Gemc->Likelihood, Mpi->tot_gemc);

        if(fabs((Varance_new-Varance)/Varance_new)<=1e-2) break;
        Varance = Varance_new;
      }

      if(fabs((Bestliklihood-likelihood)/likelihood)<=1e-2) break;

      likelihood = Bestliklihood;

      for(Indx_Chain=Mpi->indxb_gemc; Indx_Chain<=Mpi->indxe_gemc; Indx_Chain++){
        if(Indx_Chain!=Chain_bestfit){
          for(Indx_Par = 0; Indx_Par < Str_Par->Num_Par; Indx_Par++){
            if(Str_Par->narrower_guess[Indx_Par]){
              Gemc->Chains[Indx_Chain][Indx_Par] = \
                  (Ran1(Mpi->idum+Mpi->Rank)-0.5) \
                  *(Str_Par->Par_Bounds[Indx_Par][1] \
                  -Str_Par->Par_Bounds[Indx_Par][0])*0.3 \
                  +Str_Par->Par_Bounds[Indx_Par][0]*0.5 \
                  +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
            } else {
              Gemc->Chains[Indx_Chain][Indx_Par] = \
                  (Ran1(Mpi->idum+Mpi->Rank)-0.5) \
                  *(Str_Par->Par_Bounds[Indx_Par][1]  \
                  -Str_Par->Par_Bounds[Indx_Par][0]) \
                  +Str_Par->Par_Bounds[Indx_Par][0]*0.5 \
                  +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
            }
          }
          Gemc->Likelihood[Indx_Chain]=Likelihood_Log(Atom, \
              Gemc->Chains[Indx_Chain], Observation);
        }
      }

      for(i=0; i<Mpi->Num_procs; i++){
        MPI_Bcast(Gemc->Chains[i*Mpi->num_gemc], \
            Str_Par->Num_Par*Mpi->num_gemc, MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(Gemc->Likelihood+i*Mpi->num_gemc, Mpi->num_gemc, \
            MPI_DOUBLE, i, MPI_COMM_WORLD );
      }
      Chain_bestfit = BESTFIT_CHAIN(Mpi, Gemc->Likelihood, Mpi->tot_gemc, \
          Gemc->Chains, Str_Par->Num_Par, &Bestliklihood, Bestsample);
      Varance = Chain_Varance(Mpi, Gemc->Likelihood, Mpi->tot_gemc);
    }

    MPI_Allreduce(&count_accept, &tot_accept, 1, MPI_INT, MPI_SUM, \
        MPI_COMM_WORLD);
    MPI_Allreduce(&count_sample, &tot_sample, 1, MPI_INT, MPI_SUM, \
        MPI_COMM_WORLD);

    if(Mpi->Rank == 0){
      fprintf(stderr, "GEMC : total samples = %d accepted samples = %d rates = %.1f%% \n", \
        tot_sample, tot_accept, tot_accept*100./tot_sample);
      fprintf(stderr, "Best likelihood (in log scale) is %e \n", \
          Bestliklihood);
      fprintf(stderr, "Best fit parameters are \n");
      for(Indx_Par =0; Indx_Par<9; Indx_Par++){
        fprintf(stderr, "%d %e \n", Indx_Par, \
            Gemc->Chains[Chain_bestfit][Indx_Par]);
      }

    }

    Str_Par->Par_Bounds[2][0] = Gemc->Chains[Chain_bestfit][2]-3.1415926/2;
    Str_Par->Par_Bounds[2][1] = Gemc->Chains[Chain_bestfit][2]+3.1415926/2;

/*
    if(Mpi->Rank == 0){
      fprintf(stderr, "%d %d \n",count_accept,count_sample);
      fprintf(stderr, "%e \n", Gemc->Likelihood[Chain_bestfit]);
      for(Indx_Par =0; Indx_Par<9; Indx_Par++){
        fprintf(stderr, "%d %e \n", Indx_Par, Gemc->Chains[Chain_bestfit][Indx_Par]);
      }
    }
*/
    return Chain_bestfit;
}

/*--------------------------------------------------------------------------------*/

extern int Chain_Init_Gemc(STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, STRUCT_ATOM *Atom, \
    STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par){
    
    /*######################################################################
      Purpose:
        Initialize the model parameters in GEMC chains.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Gemc, a structure saved the GEMC samples and likelihood.
        Atom, a structure the atomic information
        Observation, a structure saved observed profile.
        Str_Par, a structure saved the parameter bounds.
      Output parameters:
        Gemc, a structure saved the GEMC samples and likelihood.
     ######################################################################*/
    
    int Indx_Rank, Indx_Par, Indx_Chain;
    
    Gemc->Chains = (double **)MATRIX(0, Mpi->tot_gemc-1, 0, \
        Str_Par->Num_Par-1, enum_dbl, true);
    Gemc->Likelihood = (double *)VECTOR(0, Mpi->tot_gemc-1, enum_dbl, false);

    for(Indx_Chain = Mpi->indxb_gemc; Indx_Chain <= Mpi->indxe_gemc; Indx_Chain++){
      for(Indx_Par = 0; Indx_Par < Str_Par->Num_Par; Indx_Par++){
        if(Str_Par->narrower_guess[Indx_Par]){
          Gemc->Chains[Indx_Chain][Indx_Par] = \
              (Ran1(Mpi->idum+Mpi->Rank)-0.5)*(Str_Par->Par_Bounds[Indx_Par][1] \
              -Str_Par->Par_Bounds[Indx_Par][0])*0.3 \
              +Str_Par->Par_Bounds[Indx_Par][0]*1.5 \
              +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
        } else {
            Gemc->Chains[Indx_Chain][Indx_Par] = \
              (Ran1(Mpi->idum+Mpi->Rank)-0.5)*(Str_Par->Par_Bounds[Indx_Par][1] \
                  -Str_Par->Par_Bounds[Indx_Par][0]) \
                  +Str_Par->Par_Bounds[Indx_Par][0]*1.5 \
                  +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
        }
      }
      Gemc->Likelihood[Indx_Chain]=Likelihood_Log(Atom, \
          Gemc->Chains[Indx_Chain], Observation);
    }

    for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
      MPI_Bcast(Gemc->Chains[Indx_Rank*Mpi->num_gemc], \
          Str_Par->Num_Par*Mpi->num_gemc, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
      MPI_Bcast(Gemc->Likelihood+Indx_Rank*Mpi->num_gemc, Mpi->num_gemc, \
          MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD );
    }

    return 0;    
}

/*--------------------------------------------------------------------------------*/

extern int BESTFIT_CHAIN(STRUCT_MPI *Mpi, double *Likelihood, int Num_Chains, \
    double **Chains, int Num_Par, double *Bestliklihood, double *Bestsample){

    /*######################################################################
      Purpose:
        find the bestfit chain.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Likelihood[], the likilihood in log scale.
        Num_Chains, the number of the chains.
        Chains[][], the model parameters in all the chains.
        Num_Par, the number of the model parameters.
      Output parameter,
        Bestliklihood, the bestfit likilihood in log scale.
        Bestsample, the bestfit sample parameters.
      Return:
        the index of the bestfit chain.
     ######################################################################*/

    int Chain_bestfit, Indx_Chain, Indx_Par;
    if(Mpi->Rank == 0){
      Chain_bestfit = 0;
      for(Indx_Chain=1; Indx_Chain<Num_Chains; Indx_Chain++){
        if(Likelihood[Indx_Chain]>Likelihood[Chain_bestfit]){
          Chain_bestfit = Indx_Chain;
        }
      }
    }

    MPI_Bcast(&Chain_bestfit, 1, MPI_INT, 0, MPI_COMM_WORLD);

    *Bestliklihood = Likelihood[Chain_bestfit];
    for(Indx_Par = 0; Indx_Par<Num_Par; Indx_Par++){
      Bestsample[Indx_Par] = Chains[Chain_bestfit][Indx_Par];
    }

    return Chain_bestfit;
}

/*--------------------------------------------------------------------------------*/

extern double Chain_Varance(STRUCT_MPI *Mpi, double *Likelihood, int Num_Chains){

    /*######################################################################
      Purpose:
        Compute the variance of the Likelihood.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Likelihood[], the likilihood in log scale.
        Num_Chains, the number of the chains.
      Return:
        the variance of the likilihood.
     ######################################################################*/

    double AVG_LIKELIHOOD = 0.0, SUM_SQUARE = 0.0, tmp = 0.0, Varance;
    int Indx_Chain;

    tmp = 0.0;
    for(Indx_Chain=Mpi->indxb_gemc; Indx_Chain<=Mpi->indxe_gemc; Indx_Chain++){
      tmp += Likelihood[Indx_Chain];
    }

    MPI_Allreduce(&tmp, &AVG_LIKELIHOOD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    AVG_LIKELIHOOD /= Mpi->tot_gemc;                

    tmp = 0.0;
    for(Indx_Chain=Mpi->indxb_gemc; Indx_Chain<=Mpi->indxe_gemc; Indx_Chain++){
      tmp += (Likelihood[Indx_Chain]-AVG_LIKELIHOOD) \
        *(Likelihood[Indx_Chain]-AVG_LIKELIHOOD);
    } 
            
    MPI_Allreduce(&tmp, &SUM_SQUARE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
    Varance = sqrt(SUM_SQUARE/(Num_Chains-1));

    return Varance;
}

/*--------------------------------------------------------------------------------*/

extern int Bounds_Enforce(double *Sample, STRUCT_PAR *Str_Par){
    
    /*######################################################################
      Purpose:
        Checks the bounds of the parameters.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Sample[], a vector which saves the sampled parameters.
        Str_Par, a structure saved the parameter bounds.
      Output parameters:
        Sample[], the checked parameters.
     ######################################################################*/
    
    int Indx_Par;
    double tmp;

    for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
      tmp = Str_Par->Par_Bounds[Indx_Par][1]-Str_Par->Par_Bounds[Indx_Par][0];
                
      if(tmp <= 0.){
        Sample[Indx_Par] = Str_Par->Par_Bounds[Indx_Par][0];
      }else if(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0] \
          || Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]){

        switch(Str_Par->type[Indx_Par]){
          case enum_fold://fold bounds
            if(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0]){
              do{
                Sample[Indx_Par] += tmp;
              }while(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0]);
            }else if(Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]){
              do{
                Sample[Indx_Par]-=tmp;
              }while(Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]);
            }
            break;
              
          case enum_reflect:///reflect bounds
            while(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0] \
              || Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]){
              if(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0]){              
                Sample[Indx_Par] = 2*Str_Par->Par_Bounds[Indx_Par][0] \
                    -Sample[Indx_Par];
              }
              if(Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]){              
                Sample[Indx_Par] = 2*Str_Par->Par_Bounds[Indx_Par][1] \
                    -Sample[Indx_Par];
              }              
            }
            break;
              
          case enum_set:///set to bound value
            if(Sample[Indx_Par]<Str_Par->Par_Bounds[Indx_Par][0]){              
              Sample[Indx_Par] = Str_Par->Par_Bounds[Indx_Par][0];
            }else if(Sample[Indx_Par]>Str_Par->Par_Bounds[Indx_Par][1]){
              Sample[Indx_Par] = Str_Par->Par_Bounds[Indx_Par][1];
            }
            break;
              
          default:    
            break;
        }
      }
    }
    
    return 0;   
}

/*--------------------------------------------------------------------------------*/
