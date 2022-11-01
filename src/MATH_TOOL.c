
#include "MATH_TOOL.h"

/*--------------------------------------------------------------------------------*/

extern int PSRF(STRUCT_MPI *Mpi, double ***Chains, int Gener, int Num_Par, \
    double *Statis_R){
    
    /*######################################################################
      Purpose:
        Computes the potential scale reducton factor (PSRF) used to check 
          multichain MCMC convergence (omit the correction factor due to 
          D tends to be large at convergence).
      Record of revisions:
        12 Mar 2018.
      Input parameters:
        Mpi, a structure saved the Mpi information.
        Gener, the current generation.
        Num_Par, the number of the model parameters.     
      Output parameters:
        Statis_R[], an array saved PSRF for each parameter.
      Method:
        Gelman Rubin R statistics.
      References:
        Brooks, S. and Gelman, A. (1998). General methods for monitoring 
          convergence of iterative simulations. Journal of Computational 
          and Graphical Statistics, 7(4), 434-55.
        Gelman, A. and Rubin, D. B. (1992). Inference from iterative 
          simulation using multiple sequences. Statistical Science, 7, 
          457-72.
     ######################################################################*/
    
    int begin = Gener/2;
    int Length = Gener-begin+1;
    if(Length > 1500){
      Length = 1500;
      begin = Gener+1-Length;
    }
    int Indx_Par, Indx_Chain, Indx_Gener;
    
    //double *AVG_Chain = (double *)malloc(Num_Chain*sizeof(double));
    double AVG_Chain[Mpi->tot_dream];
    double tmp, Average, Par_B, Par_W, Par_Sigma;

    for(Indx_Par=0; Indx_Par<Num_Par; Indx_Par++){      
      tmp = 0;
      for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; Indx_Chain++){
        AVG_Chain[Indx_Chain]=0;
        for(Indx_Gener=begin; Indx_Gener<=Gener; Indx_Gener++){
          AVG_Chain[Indx_Chain] += Chains[Indx_Gener][Indx_Chain][Indx_Par];
        }
        AVG_Chain[Indx_Chain] /= Length;

        tmp += AVG_Chain[Indx_Chain];
      }
      /*
      for(Indx_Rank=0; Indx_Rank<Mpi->Num_procs; Indx_Rank++){
        MPI_Bcast(AVG_Chain+Indx_Rank*Mpi->num_dream, Mpi->num_dream, MPI_DOUBLE, \
          Indx_Rank, MPI_COMM_WORLD );
      }*/
      MPI_Allreduce(&tmp, &Average, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      Average /= Mpi->tot_dream;
        
      tmp = 0;
      for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; Indx_Chain++){
        tmp += (AVG_Chain[Indx_Chain]-Average)*(AVG_Chain[Indx_Chain]-Average);
      }

      MPI_Allreduce(&tmp, &Par_B, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      Par_B /= (Mpi->tot_dream-1.);

      tmp = 0;
      for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; Indx_Chain++){
        for(Indx_Gener=begin; Indx_Gener<=Gener; Indx_Gener++){
          tmp += (Chains[Indx_Gener][Indx_Chain][Indx_Par]-AVG_Chain[Indx_Chain])
              *(Chains[Indx_Gener][Indx_Chain][Indx_Par]-AVG_Chain[Indx_Chain]);

        }
      }
      MPI_Allreduce(&tmp, &Par_W, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      Par_W = Par_W/(Length-1.0)/Mpi->tot_dream;
      Par_Sigma = (Length-1.)/Length*Par_W+Par_B;
      Statis_R[Indx_Par] = sqrt((Mpi->tot_dream+1.0)/Mpi->tot_dream \
          *Par_Sigma/Par_W-(Length-1.0)/Length/Mpi->tot_dream);

    }
    
    return 0; 
}

/*--------------------------------------------------------------------------------*/

extern int Chains_STD(STRUCT_MPI *Mpi, int Begin_Gener, int End_Gener, \
    int Num_Par, double ***Chains, double *STD, double *Mean){
    
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
    
    int Indx_Par, Indx_Gener, Indx_Chain;
    
    double AVG, SUM, tot_AVG, tot_SUM;
    for(Indx_Par=0; Indx_Par<Num_Par; Indx_Par++){
      AVG=0.;  
      SUM=0.;
        
      for(Indx_Gener=Begin_Gener; Indx_Gener<=End_Gener; Indx_Gener++){      
        for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
            Indx_Chain++){
          AVG += Chains[Indx_Gener][Indx_Chain][Indx_Par];
        }
      }
      MPI_Allreduce(&AVG, &tot_AVG, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      Mean[Indx_Par] = tot_AVG/(End_Gener-Begin_Gener+1)/Mpi->tot_dream;
      for(Indx_Gener=Begin_Gener; Indx_Gener<=End_Gener; Indx_Gener++){      
        for(Indx_Chain=Mpi->indxb_dream; Indx_Chain<=Mpi->indxe_dream; \
            Indx_Chain++){      
          SUM += (Chains[Indx_Gener][Indx_Chain][Indx_Par]-Mean[Indx_Par]) \
            *(Chains[Indx_Gener][Indx_Chain][Indx_Par]-Mean[Indx_Par]);
        }
      }
      MPI_Allreduce(&SUM, &tot_SUM, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      STD[Indx_Par] = sqrt(tot_SUM/((End_Gener-Begin_Gener+1) \
          *Mpi->tot_dream-1));  
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Chains_STD_Single(int Num_Chain, int Begin_Gener, int End_Gener, \
    int Num_Par, double ***Chains, double *STD, double *Mean){
    
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
    
    int Indx_Par, Indx_Gener, Indx_Chain;
    double AVG, SUM;

    for(Indx_Par=0; Indx_Par<Num_Par; Indx_Par++){
      AVG = 0.;  
      SUM = 0.;
        
      for(Indx_Gener=Begin_Gener; Indx_Gener<=End_Gener; Indx_Gener++){      
        for(Indx_Chain=0; Indx_Chain<Num_Chain; Indx_Chain++){
          AVG += Chains[Indx_Gener][Indx_Chain][Indx_Par];
        }
      }

      Mean[Indx_Par] = AVG/(End_Gener-Begin_Gener+1)/Num_Chain;
      for(Indx_Gener=Begin_Gener; Indx_Gener<=End_Gener; Indx_Gener++){      
        for(Indx_Chain=0; Indx_Chain<Num_Chain; Indx_Chain++){      
          SUM += (Chains[Indx_Gener][Indx_Chain][Indx_Par]-Mean[Indx_Par]) \
            *(Chains[Indx_Gener][Indx_Chain][Indx_Par]-Mean[Indx_Par]);
        }
      }
      STD[Indx_Par] = sqrt(SUM/((End_Gener-Begin_Gener+1)*Num_Chain-1));  
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/
