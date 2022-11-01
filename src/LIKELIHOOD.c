#include "LIKELIHOOD.h"

/*--------------------------------------------------------------------------------*/

double Likelihood_Log(STRUCT_ATOM *Atom, double *Par, STRUCT_OBSERVATION *Observation){
    
    /*######################################################################
      Purpose:
        Compute the log-likelihood of Stokes profiles.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        Atom, a structure saved the Landu facor.
        Par, an array saved the parameters.
        Observation, a structure saved the lambda.
      Return:
        return the log-likelihood.
     ######################################################################*/
    
    int i, j;
    double Likelihood = 0;
    
    Forward_ME(Atom, Par, Observation);

    double weights[4];
    if (Observation->weights_flag){
      weights[0] = 1;
      weights[1] = 3;
      weights[2] = 3;
      weights[3] = 2;
      weights[0] = 1;
      weights[1] = 1;
      weights[2] = 1;
      weights[3] = 1;
    } else {
      weights[0] = 1;
      weights[1] = 1;
      weights[2] = 1;
      weights[3] = 1;
    }

    
    for (i=0; i<Observation->Num_Lam; i++) {
      for (j=0; j<4; j++){
        Likelihood += (Observation->Profile_Data[j][i]-Observation->Stokes[j][i]) \
          *(Observation->Profile_Data[j][i]-Observation->Stokes[j][i]) \
          /Observation->Sigma[j][i]/Observation->Sigma[j][i]*weights[j]*weights[j];
      }     
    }
    
    return -0.5*Likelihood;
}

/*--------------------------------------------------------------------------------*/

