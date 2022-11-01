
#include "DREAM.h"

/*--------------------------------------------------------------------------------*/

extern int DREAM(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, STRUCT_DREAM *Dream, \
    STRUCT_ATOM *Atom, STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par){

    /*######################################################################
      Purpose:
        DREAM simulation.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Input, a structure saved the input information.
        Mpi, a structure saved the Mpi information.
        Dream, a structure saved the DREAM samples and likelihood.
        Atom, a structure the atomic information
        Observation, a structure saved observed profile.
        Str_Par, a structure saved the parameter bounds.
      Output parameters:
        Dream, a structure saved the DREAM samples and likelihood.
      Returm:
        the index of the bestfit chains.
     ######################################################################*/

    if(Mpi->Rank ==0){ 
      fprintf(stderr, "----------DREAM SAMPLING----------\n");
      fprintf(stderr, "generation number = %d, chain number = %d \n", \
        Input->Num_Gener, Mpi->tot_dream);
    }

    STRUCT_CR *Cr = (STRUCT_CR *)malloc(sizeof(STRUCT_CR));
    Cr->Num_Cr = Input->Num_Cr;

    bool Converged = false;
    int Num_Pair, Indx_Cr;
    int count_sample=0, count_accept=0, tot_sample=0, tot_accept=0;
    int Indx_Chain=0, Indx_Gener, Indx_Par, Indx_Rank;
    int Burn_in, Chain_bestfit;
    double TMP_LIKELIHOOD, Bestliklihood;
    double Sample[Str_Par->Num_Par], Statis_R[Str_Par->Num_Par];
    double Bestsample[Str_Par->Num_Par];
    double STD[Str_Par->Num_Par], Mean[Str_Par->Num_Par];

    Init_Cr(Cr);

    for(Burn_in = 0; Burn_in<2; Burn_in++){
      count_sample=0;
      count_accept=0;
      for(Indx_Gener=1; Indx_Gener<Input->Num_Gener; Indx_Gener++){
        for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
            Indx_Chain++){
         
          Num_Pair = Sample_PairNum(Input->MaxNum_Pair, Mpi);
          Indx_Cr = Dream_Sample(Mpi, Cr, Str_Par, Num_Pair, \
              Dream->Chains[Indx_Gener-1], Indx_Chain, Indx_Gener, Sample);
          count_sample++;
          TMP_LIKELIHOOD = Likelihood_Log(Atom, Sample, Observation);

          if(log(Ran1(Mpi->idum+Mpi->Rank))< \
              (TMP_LIKELIHOOD-Dream->Likelihood[Indx_Gener-1][Indx_Chain])){    
            count_accept++;
            Dream->Likelihood[Indx_Gener][Indx_Chain] = TMP_LIKELIHOOD;
            for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
              Dream->Chains[Indx_Gener][Indx_Chain][Indx_Par] = Sample[Indx_Par];
            }
                  
          }else{
            Dream->Likelihood[Indx_Gener][Indx_Chain] = \
                Dream->Likelihood[Indx_Gener-1][Indx_Chain];
            for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
              Dream->Chains[Indx_Gener][Indx_Chain][Indx_Par] = \
                  Dream->Chains[Indx_Gener-1][Indx_Chain][Indx_Par];
            }
          }

          if(Cr->Num_Cr>1 && Burn_in == 0 && Input->Cr_Update && Indx_Gener<1500){
            Cr_distance(Dream->Chains, Mpi->tot_dream, Str_Par->Num_Par, \
                Indx_Chain, Indx_Gener, Indx_Cr, Cr);
          }
        }

        for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
          MPI_Bcast(Dream->Chains[Indx_Gener][Indx_Rank*Mpi->num_dream], \
              Str_Par->Num_Par*Mpi->num_dream, MPI_DOUBLE, Indx_Rank, \
              MPI_COMM_WORLD);

          MPI_Bcast(Dream->Likelihood[Indx_Gener]+Indx_Rank*Mpi->num_dream, \
              Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD );
        }

        if(Indx_Gener%500 == 0){
          if(Mpi->Rank == 0) fprintf(stderr, "generation = %d \n", Indx_Gener);
        }

        if(Indx_Gener%30 == 0 && Burn_in == 0){
        
          PSRF(Mpi, Dream->Chains, Indx_Gener, Str_Par->Num_Par, Statis_R);

          if(!Converged){
            Converged = true;   
            for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){   
              if(Statis_R[Indx_Par]>=1.2){
                Converged = false;          
                break;
              }          
            }
          }

          if(Input->Debug || (Converged)){
            Chain_bestfit = BESTFIT_CHAIN(Mpi, Dream->Likelihood[Indx_Gener], \
              Mpi->tot_dream, Dream->Chains[Indx_Gener], Str_Par->Num_Par, \
              &Bestliklihood, Bestsample);
            if(Mpi->Rank == 0) fprintf(stderr, \
              "Bestfit chains: generation = %d chain index = %d likelihood = %e\n", \
              Indx_Gener, Chain_bestfit, Bestliklihood);
          }

          if(Converged){
            if(Mpi->Rank == 0){
              fprintf(stderr, "burn in perid finished! \n Best fit results : \n");
              for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par;Indx_Par++){
                fprintf(stderr, " %e ", \
                    Dream->Chains[Indx_Gener-1][Chain_bestfit][Indx_Par]);
              }
              fprintf(stderr, "\n log likelihood = %e \n ", \
                  Dream->Likelihood[Indx_Gener-1][Chain_bestfit]);
            }
            break;
          }

          if(Cr->Num_Cr > 1 && Input->Cr_Update &&Indx_Gener<1500) Cr_Prob(Cr, Mpi);
          if(Indx_Gener < 2000)Rm_Outlierchain(Mpi, Dream, Str_Par->Num_Par, \
              Indx_Gener);
	
        }
      }

      if(Converged && Burn_in == 0){
        for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
            Indx_Chain++){
          for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par;Indx_Par++){
            Dream->Chains[0][Indx_Chain][Indx_Par] = \
                Dream->Chains[Indx_Gener-1][Indx_Chain][Indx_Par];
          }
          Dream->Likelihood[0][Indx_Chain] = \
              Dream->Likelihood[Indx_Gener-1][Indx_Chain];
        }

        for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
          MPI_Bcast(Dream->Chains[0][Indx_Rank*Mpi->num_dream], \
              Str_Par->Num_Par*Mpi->num_dream, MPI_DOUBLE, Indx_Rank, \
              MPI_COMM_WORLD);
          MPI_Bcast(Dream->Likelihood[0]+Indx_Rank*Mpi->num_dream, \
              Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
        }
      
      }else if( Burn_in == 1){
        Chain_bestfit = BESTFIT_CHAIN(Mpi, Dream->Likelihood[Indx_Gener-1], \
            Mpi->tot_dream, Dream->Chains[Indx_Gener-1], \
            Str_Par->Num_Par, &Bestliklihood, Bestsample);
        if(Mpi->Rank == 0) fprintf(stderr, \
            "Bestfit chains: generation = %d chain index = %d likelihood = %e\n", \
            Indx_Gener, Chain_bestfit, Bestliklihood);
        
	    }else{
        for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
            Indx_Chain++){
          for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par;Indx_Par++){
            Dream->Chains[0][Indx_Chain][Indx_Par] = \
                Dream->Chains[Indx_Gener-1][Indx_Chain][Indx_Par];
          }
          Dream->Likelihood[0][Indx_Chain] = \
              Dream->Likelihood[Indx_Gener-1][Indx_Chain];
        }

        for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
          MPI_Bcast(Dream->Chains[0][Indx_Rank*Mpi->num_dream], \
              Str_Par->Num_Par*Mpi->num_dream, MPI_DOUBLE, \
              Indx_Rank, MPI_COMM_WORLD);
          MPI_Bcast(Dream->Likelihood[0]+Indx_Rank*Mpi->num_dream, \
              Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
        }
      }

      MPI_Allreduce(&count_sample, &tot_sample, 1, MPI_INT, \
          MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&count_accept, &tot_accept, 1, MPI_INT, \
          MPI_SUM, MPI_COMM_WORLD);

      tot_accept += Mpi->tot_dream;
      tot_sample += Mpi->tot_dream;
	    if(Mpi->Rank == 0) {
        fprintf(stderr, "accepted samples: %d \n", tot_accept);
        fprintf(stderr, "total samples: %d \n", tot_sample);
        fprintf(stderr, "acceptance rated: %.1f %%\n", \
            tot_accept*100./tot_sample);
      }
      if(Input->Distribution_Update && Burn_in == 0){
        Chains_STD(Mpi, Indx_Gener*3/4, Indx_Gener-1, \
            Str_Par->Num_Par, Dream->Chains, STD, Mean);
        if(Mpi->Rank == 0){
          for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
            fprintf(stderr, "%d %e %e %e %e\n", Indx_Par, \
                Dream->Chains[Indx_Gener-1][Chain_bestfit][Indx_Par], \
                Mean[Indx_Par], STD[Indx_Par], STD[Indx_Par] \
                /(Str_Par->Par_Bounds[Indx_Par][1] \
                -Str_Par->Par_Bounds[Indx_Par][0]));
            Str_Par->Distribution[Indx_Par] = STD[Indx_Par];
          }
        }
      }      
    }

    FREE_VECTOR(Cr->Cr, 1, enum_dbl);
    FREE_VECTOR(Cr->Prob, 1, enum_dbl);
    FREE_VECTOR(Cr->Delta, 1, enum_dbl);
    FREE_VECTOR(Cr->Delta_tot, 1, enum_dbl);
    FREE_VECTOR(Cr->Delta_sum, 1, enum_dbl);
    FREE_VECTOR(Cr->counts, 1, enum_int);
    FREE_VECTOR(Cr->counts_tot, 1, enum_int);
    FREE_VECTOR(Cr->counts_sum, 1, enum_int);
    free(Cr);

    return Chain_bestfit;
}

/*--------------------------------------------------------------------------------*/

extern int Chain_Init_Dream(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_DREAM *Dream, STRUCT_ATOM *Atom, \
    STRUCT_OBSERVATION *Observation, STRUCT_PAR *Str_Par){
    
    /*######################################################################
      Purpose:
        Initialize the model parameters in DREAM chains.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Input, a structure saved the input information.
        Mpi, a structure saved the Mpi information.
        Dream, a structure saved the DREAM samples and likelihood.
        Atom, a structure the atomic information
        Observation, a structure saved observed profile.
        Str_Par, a structure saved the parameter bounds.
      Output parameters:
        Dream, a structure saved the DREAM samples and likelihood.
     ######################################################################*/

    int Indx_Rank, Indx_Par, Indx_Chain;

    Dream->Chains = (double ***)TENSOR_DBL(0, Input->Num_Gener-1, 0, \
        Mpi->tot_dream-1, 0, Str_Par->Num_Par-1, true);
    Dream->Likelihood = (double **)MATRIX(0, Input->Num_Gener-1, 0, \
        Mpi->tot_dream-1, enum_dbl, false);

    for(Indx_Chain = Mpi->indxb_dream; Indx_Chain <= Mpi->indxe_dream; \
        Indx_Chain++){
      for(Indx_Par = 0; Indx_Par < Str_Par->Num_Par; Indx_Par++){
        if(Str_Par->narrower_guess[Indx_Par]){
          Dream->Chains[0][Indx_Chain][Indx_Par] = \
              (Ran1(Mpi->idum+Mpi->Rank)-0.5) \
              *(Str_Par->Par_Bounds[Indx_Par][1] \
              -Str_Par->Par_Bounds[Indx_Par][0])*0.3 \
              +Str_Par->Par_Bounds[Indx_Par][0]*1.5 \
              +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
        }else{
           Dream->Chains[0][Indx_Chain][Indx_Par] = \
              (Ran1(Mpi->idum+Mpi->Rank)-0.5) \
              *(Str_Par->Par_Bounds[Indx_Par][1] \
              -Str_Par->Par_Bounds[Indx_Par][0]) \
              +Str_Par->Par_Bounds[Indx_Par][0]*1.5\
                +Str_Par->Par_Bounds[Indx_Par][1]*0.5;
        }
      }
      Dream->Likelihood[0][Indx_Chain] = Likelihood_Log(Atom, \
          Dream->Chains[0][Indx_Chain], Observation);
    }

    for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
      MPI_Bcast(Dream->Chains[0][Indx_Rank*Mpi->num_dream], \
          Str_Par->Num_Par*Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
      MPI_Bcast(Dream->Likelihood[0]+Indx_Rank*Mpi->num_dream, \
          Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD );
    }

    return 0;    
}

/*--------------------------------------------------------------------------------*/

extern int GEMC2DREAM(STRUCT_INPUT *Input, STRUCT_MPI *Mpi, STRUCT_GEMC *Gemc, \
    STRUCT_DREAM *Dream, STRUCT_PAR *Str_Par, \
    STRUCT_OBSERVATION *Observation, STRUCT_ATOM *Atom){

    /*######################################################################
      Purpose:
        Initialized the DREAM chains with the results from GEMC.
      Record of revisions:
        30 Otc. 2022.
      Input parameters:
        Input, a structure saved the input information.
        Mpi, a structure saved the Mpi information.
        Gemc, a structure saved the GEMC samples and likelihood.
        Dream, a structure saved the DREAM samples and likelihood.
        Str_Par, a structure saved the model parameters, bounds.
        Observation, a structure saved the lambda.
        Atom, a structure saved the Landu facor.
      Output parameters:
        Sample, the sampled model parameters.
     ######################################################################*/

    int indx[Mpi->tot_gemc];
        
    if(Mpi->Rank == 0) qsort_index(0, Mpi->tot_gemc-1, \
        Gemc->Likelihood, indx);

    MPI_Bcast(indx, Mpi->tot_gemc, MPI_INT, 0, MPI_COMM_WORLD);
    Dream->Chains = (double ***)TENSOR_DBL(0, Input->Num_Gener-1, 0, \
      Mpi->tot_dream-1, 0, Str_Par->Num_Par-1, true);
    Dream->Likelihood = (double **)MATRIX(0, Input->Num_Gener-1, 0, \
      Mpi->tot_dream-1, enum_dbl, false);

    int Indx_Rank, Indx_Par, Indx_Chain, i;
    
    for(Indx_Chain=0, i = Mpi->tot_gemc-Mpi->tot_dream; \
      Indx_Chain<Mpi->tot_dream; Indx_Chain++, i++){
      for(Indx_Par = 0; Indx_Par<Str_Par->Num_Par;Indx_Par++){
        Dream->Chains[0][Indx_Chain][Indx_Par] = \
            Gemc->Chains[indx[i]][Indx_Par];
      }
      Dream->Likelihood[0][Indx_Chain] = Gemc->Likelihood[indx[i]];
    }

    if(Observation->weights_flag){
      Observation->weights_flag = false;
      for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
          Indx_Chain++){
        Dream->Likelihood[0][Indx_Chain]=Likelihood_Log(Atom, \
            Dream->Chains[0][Indx_Chain], Observation);
      }
      for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
        MPI_Bcast(Dream->Likelihood[0]+Indx_Rank*Mpi->num_dream, \
            Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD );
      }
    }
      
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Dream_Sample(STRUCT_MPI *Mpi, STRUCT_CR *Cr, STRUCT_PAR *Str_Par, \
    int Num_Pair, double **Chains, int Indx_Chain, int indx_Gener, \
    double *Sample){
    
    /*######################################################################
      Purpose:
        Sample the model parameters.
      Record of revisions:
        30 Otc. 2022.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Cr, a structure saved the Cr information.
        Str_Par, a structure saved the model parameters, bounds.
        Num_Pair, the number of pairs of crossover chains.
        Chains[][], the sampled model parameters in the specified.
        Indx_Chain, index of the chains.
      Output parameters:
        Sample, the sampled model parameters.
     ######################################################################*/
    
    int *Jump_dim = VECTOR(0, Str_Par->Num_Par-1, enum_int, false);
    double *Diff = VECTOR(0, Str_Par->Num_Par-1, enum_dbl, false);

    int Indx_Cr = Sample_Cr(Cr, Mpi);

    int Num_Jump = Dream_Dim(Mpi, Cr, Str_Par, Indx_Cr, Jump_dim);

    Dream_Diff(Mpi, Chains, Indx_Chain, Num_Pair, Num_Jump, Jump_dim, Diff);

    double *Factor_E = VECTOR(0, Num_Jump-1, enum_dbl, false);
    double *Factor_Epsilon = VECTOR(0, Num_Jump-1, enum_dbl, false);

    int i, Indx_Par;
    double gamma; 

    for(i=0; i<Num_Jump; i++){    
      Factor_E[i]=(Ran1(Mpi->idum+Mpi->Rank)-0.5)*0.1;
      Factor_Epsilon[i]=GASDEV(Mpi->idum+Mpi->Rank) \
          *Str_Par->Distribution[Jump_dim[i]]*1e-5;
    }
      
    for(Indx_Par=0; Indx_Par<Str_Par->Num_Par; Indx_Par++){
      Sample[Indx_Par] = Chains[Indx_Chain][Indx_Par];
    }
    
    if(indx_Gener%5==0){
      gamma = 1;
    }else if(Str_Par->Num_Par > Num_Jump){
      //gamma = 2.38/sqrt(2.*Num_Pair*(Str_Par->Num_Par-Num_Jump));
      gamma = 2.38/sqrt(2.*Num_Pair*(Str_Par->Num_Par));
    }else{
      gamma = 2.38/sqrt(2.*Num_Pair*(Str_Par->Num_Par));
    }
    
    //Num_Jump
    for(i=0; i<Num_Jump; i++){
      Sample[Jump_dim[i]] += Diff[i]*(1.0+Factor_E[i])*gamma \
        +Factor_Epsilon[i]; 
    }

    Bounds_Enforce(Sample, Str_Par);
    
    FREE_VECTOR(Jump_dim, 0, enum_int);
    FREE_VECTOR(Diff, 0, enum_dbl);
    FREE_VECTOR(Factor_E, 0, enum_dbl);
    FREE_VECTOR(Factor_Epsilon, 0, enum_dbl);

    return Indx_Cr;
}

/*--------------------------------------------------------------------------------*/

extern int Dream_Dim(STRUCT_MPI *Mpi, STRUCT_CR *Cr, STRUCT_PAR *Str_Par, \
    int Indx_Cr, int *Jump_dim){
    
    /*######################################################################
      Purpose:
        Choose the dimensions in which a jump is to be made.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Cr, a structure saved the Cr information.
        Str_Par, a structure saved the model parameters, bounds.
        Indx_Cr, index of the choosed Cr.
      Output parameters:
        Jump_dim, a vector which saves the indexs of jumping parameters.
     ######################################################################*/
    
    int i=0, Num_Jump=0;
    
    for(i=0; i<Str_Par->Num_Par; i++){
      Jump_dim[i] = -1;  
    }
    
    for(i=0; i<Str_Par->Num_Par; i++){
      if(Ran1(Mpi->idum+Mpi->Rank)<=1.0-Cr->Cr[Indx_Cr]){
        Jump_dim[Num_Jump] = i;
        Num_Jump++;
      }  
    }
    
    if(Num_Jump==0){
        Jump_dim[0]=(int)(Ran1(Mpi->idum+Mpi->Rank)*Str_Par->Num_Par);
        Num_Jump = 1;
    }
    
    return Num_Jump;
    
}

/*--------------------------------------------------------------------------------*/

extern int Dream_Diff(STRUCT_MPI *Mpi, double **Chains, int Indx_Chain, \
    int Num_Pair, int Num_Jump, int *Jump_dim, double *Diff){
    
    /*######################################################################
      Purpose:
        Calculate the pairs differences used to sample new candidates.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Chains[][], the sampled model parameters in the specified 
          generation.
        Indx_Chain, the index of the current chain.
        Num_Pair, the number of pairs of crossover chains.
        Num_Jump, the number of dimensions in which a jump will be made.
        Jump_dim, the dimensions in which a jump is to be made.
      Output parameters:
        Diff[], a vector which saves pairs differences.
     ######################################################################*/
        
    int Pairs[2*Num_Pair];
    int i, j;
    int tmp1, tmp2;

    for(i=0; i<Num_Pair; i++){
      do{      
        tmp1 = (int)(Ran1(Mpi->idum+Mpi->Rank)*Mpi->tot_dream);    
        tmp2 = (int)(Ran1(Mpi->idum+Mpi->Rank)*Mpi->tot_dream);
      }while(tmp1==tmp2||tmp1==Indx_Chain||tmp2==Indx_Chain);
        
      Pairs[i*2] = tmp1;      
      Pairs[i*2+1] = tmp2;
    }

    for(i=0; i<Num_Jump; i++){   
      Diff[i] = 0;
  //    if (Mpi->Rank ==0) fprintf(stderr, "%d ", Jump_dim[i]);
      for(j=0; j<Num_Pair; j++){
        Diff[i] += Chains[Pairs[j*2]][Jump_dim[i]] \
            -Chains[Pairs[j*2+1]][Jump_dim[i]];
     //   if (Mpi->Rank ==0) fprintf(stderr, "%e %e %e ", Chains[Pairs[j*2]][Jump_dim[i]], Chains[Pairs[j*2+1]][Jump_dim[i]],Diff[i]);
      }     
     // if (Mpi->Rank ==0)fprintf(stderr, "\n ");
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Cr_distance(double ***Chains, int Num_Chain, int Num_Par, \
    int Indx_Chain, int Indx_Gener, int Indx_Cr, STRUCT_CR *Cr){
    
    /*######################################################################
      Purpose:
        Compute the squared normalized jumping distance.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Chains[][][], the sampled model parameters.
        Num_Chain, the number of chains
        Num_Par, the number of model parameters.
        Indx_Chain, the index of current chain.
        Indx_Gener, the index of current generation.
        Indx_Cr, the index of CR value choosed.
        Cr, a structure saved the Cr information.
      Output parameters:
        Cr, a structure saved the Cr information.
     ######################################################################*/
    
    int Indx_Par;
    double STD[Num_Par], Mean[Num_Par];

    Chains_STD_Single(Num_Chain, Indx_Gener-1, Indx_Gener-1, Num_Par, Chains, \
        STD, Mean);
    
    for(Indx_Par=0; Indx_Par<Num_Par; Indx_Par++){ 
      Cr->Delta[Indx_Cr] += pow((Chains[Indx_Gener][Indx_Chain][Num_Par] \
        -Chains[Indx_Gener-1][Indx_Chain][Num_Par])/STD[Indx_Par],2);
    }
    
    Cr->counts[Indx_Cr]++;

    return 0 ; 
}

/*--------------------------------------------------------------------------------*/

extern int Cr_Prob(STRUCT_CR *Cr, STRUCT_MPI *Mpi){
    
    /*######################################################################
      Purpose:
        Update the probability of each individual CR values.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Cr, a structure saved the Cr information.
        Mpi, a structure saved the Mpi information.
      Output parameters:
        Cr, a structure saved the Cr information.

     ######################################################################*/
    
    int i;
    double Tot_Prob=0, Tot_Dis=0;

    MPI_Allreduce(Cr->Delta+1, Cr->Delta_sum+1, Cr->Num_Cr, MPI_DOUBLE, \
        MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(Cr->counts+1, Cr->counts_sum+1, Cr->Num_Cr, MPI_INT, \
        MPI_SUM, MPI_COMM_WORLD);

    if(Mpi->Rank == 0){
      for(i=1; i<=Cr->Num_Cr; i++){
        Cr->Delta_tot[i] += Cr->Delta_sum[i];
        Cr->counts_tot[i] += Cr->counts_sum[i];
        Tot_Dis += Cr->Delta_tot[i];
      }

      for(i=1; i<=Cr->Num_Cr; i++){
        Cr->Prob[i] = Cr->Delta_tot[i]/Cr->counts_tot[i]/Tot_Dis;
        Tot_Prob += Cr->Prob[i];
      }
      for(i=1; i<=Cr->Num_Cr; i++){        
        Cr->Prob[i] /= Tot_Prob;
      }
    }

    MPI_Bcast(Cr->counts_tot+1, Cr->Num_Cr, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(Cr->Delta_tot+1, Cr->Num_Cr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Cr->Prob+1, Cr->Num_Cr, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(i=1; i<=Cr->Num_Cr; i++){
      Cr->Delta[i] = 0;
      Cr->counts[i] = 0;
    }
   
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Init_Cr(STRUCT_CR *Cr){
    
    /*######################################################################
      Purpose:
        Initialize the CR (crossover probability) values.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        .
      Output parameters:
        Cr, a structure saved the Cr information.
     ######################################################################*/
    
    Cr->Cr=VECTOR(1, Cr->Num_Cr, enum_dbl, false);
    Cr->Prob=VECTOR(1, Cr->Num_Cr, enum_dbl, false);
    Cr->Delta=VECTOR(1, Cr->Num_Cr, enum_dbl, false);
    Cr->counts=VECTOR(1, Cr->Num_Cr, enum_int, false);
    Cr->Delta_sum=VECTOR(1, Cr->Num_Cr, enum_dbl, false);
    Cr->counts_sum=VECTOR(1, Cr->Num_Cr, enum_int, false);
    Cr->Delta_tot=VECTOR(1, Cr->Num_Cr, enum_dbl, false);
    Cr->counts_tot=VECTOR(1, Cr->Num_Cr, enum_int, false);

    int i;

    for(i=1; i<=Cr->Num_Cr; i++){
      Cr->Cr[i] = ((double)(i))/Cr->Num_Cr;
      Cr->Delta[i] = 0;
      Cr->Delta_tot[i] = 1;
      Cr->Prob[i] = 1.0/Cr->Num_Cr;
      Cr->counts[i] = 0;
    }
    
    return 0 ;  
}

/*--------------------------------------------------------------------------------*/

extern int Sample_Cr(STRUCT_CR *Cr, STRUCT_MPI *Mpi){
    
    /*######################################################################
      Purpose:
        Sample a CR index from the total number of Cr values according to 
          the probability of each individual CR values.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Cr, a structure saved the Cr information.
        Mpi, a structure saved the Mpi information.
      Output parameters:
        .
      Return:
        the CR index
     ######################################################################*/
    
    int i;
    double tmp = Ran1(Mpi->idum+Mpi->Rank);
    
    for(i=1; i<=Cr->Num_Cr; i++){
      tmp -= Cr->Prob[i];
      if(tmp<0) return i;
    }
    
    return Cr->Num_Cr; 
}

/*--------------------------------------------------------------------------------*/

extern int Sample_PairNum(int MaxNum_Pair, STRUCT_MPI *Mpi){
    
    /*######################################################################
      Purpose:
        Sample the number of chain pairs used to generate the jump.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        MaxNum_Pair, the max number of chain pairs used to generate the jump.
        Mpi, a structure saved the Mpi information.
      Output parameters:
        .
      Return:
        the number of chain pairs used to generate the jump.
     ######################################################################*/
    
    int i, Num_Pair=1;
    double tmp = Ran1(Mpi->idum+Mpi->Rank);
    
    for(i=1; i<=MaxNum_Pair; i++){
      if(tmp<1.0/MaxNum_Pair){
        Num_Pair=i;
        break;
      }else{
        tmp-=1.0/MaxNum_Pair;
      }
    }
    
    return Num_Pair; 
}

/*--------------------------------------------------------------------------------*/

extern int Rm_Outlierchain(STRUCT_MPI *Mpi, STRUCT_DREAM *Dream, int Num_Par, \
    int Gener){

    /*######################################################################
      Purpose:
        Remove the outlier chains.
      Record of revisions:
        30 Otc. 2022
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Dream, a structure saved the samples and likelihood.
        Num_Par, the number of the model parameters.
        Gener, the current generation.
      Output parameters:
        Dream, a structure saved the samples and likelihood.
      Return:
        the number of outlier chains.
      Reference:
        Vrugt, J. A., et al. 2009. International Journal of Nonlinear 
          Science & Numerical Simulation, 10(3), 273-290
     ######################################################################*/
    
    int Begin = Gener/2;
    int Length = Gener-Begin+1;
    double AVG_likelihood[Mpi->tot_dream];
    int indx[Mpi->tot_dream];

    int Indx_Chain, Indx_Gener, Indx_Par, Indx_Rank;

    for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; Indx_Chain++){
	    AVG_likelihood[Indx_Chain] = 0;
      for(Indx_Gener=Begin; Indx_Gener<=Gener; Indx_Gener++){
        AVG_likelihood[Indx_Chain] += Dream->Likelihood[Indx_Gener][Indx_Chain];
      }
      AVG_likelihood[Indx_Chain] /= Length;
    }

    for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
      MPI_Bcast(AVG_likelihood+Indx_Rank*Mpi->num_dream, Mpi->num_dream, \
        MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
    }

    if(Mpi->Rank == 0) qsort_index(0, Mpi->tot_dream-1, AVG_likelihood, indx);
    MPI_Bcast(indx, Mpi->tot_dream, MPI_INT, 0, MPI_COMM_WORLD);
    
    double Qr = AVG_likelihood[indx[Mpi->tot_dream*3/4]] \
        -AVG_likelihood[indx[Mpi->tot_dream/4]];
    double Omega = AVG_likelihood[indx[Mpi->tot_dream/4]]-Qr*2;
    
    int Num_Outlierchain=0;


    for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; Indx_Chain++){
      if(AVG_likelihood[Indx_Chain]<Omega){
        Num_Outlierchain++;
               
        for(Indx_Par=0; Indx_Par<Num_Par; Indx_Par++){
          Dream->Chains[Gener][Indx_Chain][Indx_Par] = \
              Dream->Chains[Gener][indx[Mpi->tot_dream-1]][Indx_Par];
        }
            
        for(Indx_Gener=0; Indx_Gener<=Gener; Indx_Gener++){
          Dream->Likelihood[Indx_Gener][Indx_Chain] = \
              Dream->Likelihood[Indx_Gener][indx[Mpi->tot_dream-1]];
        }     
      }        
    }

    for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
      MPI_Bcast(Dream->Chains[Gener][Indx_Rank*Mpi->num_dream], \
          Num_Par*Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
      for(Indx_Gener=0; Indx_Gener<=Gener; Indx_Gener++){
        MPI_Bcast(Dream->Likelihood[Indx_Gener]+Indx_Rank*Mpi->num_dream, \
            Mpi->num_dream, MPI_DOUBLE, Indx_Rank, MPI_COMM_WORLD);
      }
    }
    
    if(Mpi->Rank == 0 && Num_Outlierchain >0) \
        fprintf(stderr,"Generation: %d, Outlier chains = %d \n", Gener, \
            Num_Outlierchain);
 
    return Num_Outlierchain;    
}

/*--------------------------------------------------------------------------------*/
