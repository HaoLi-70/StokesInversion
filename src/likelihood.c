#include "likelihood.h"
#include <math.h>

/*--------------------------------------------------------------------------------*/

double Likelihood_Log(double *model_values, STRUCT_STK *stokes, 
  STRUCT_PARA *params, int *numerical_error){
    
    /*######################################################################
      Purpose:
        Compute the log-likelihood of Stokes profiles.
      Input parameters:
        model_values, proposed model parameter values.
        stokes, observed profile, synthesis workspace, and inverse noise.
        params, spectral-line data and model configuration.
      Output parameters:
        numerical_error, set to -1 for invalid state or a non-finite value.
      Return:
        Gaussian log-likelihood, equal to -0.5 times chi-squared.
     ######################################################################*/
    
    if(numerical_error) *numerical_error = 0;
    if(!model_values || !stokes || !params || !numerical_error){
      if(numerical_error) *numerical_error = -1;
      return -INFINITY;
    }
    /* Model layout, fixed values, profile, and inverse-noise weights are
       validated once by Model_Layout_Init(), Profile_Prepare(), and
       Noise_Init().  Do not repeat those invariant checks for every MCMC
       proposal. */
    double full_model[MAX_MODEL_PARAMS];
    int ifree = 0;
    for(int imodel=0; imodel<params->nmodel; imodel++){
      if(params->inv[imodel]){
        full_model[imodel] = model_values[ifree++];
      }else{
        full_model[imodel] = params->value_const[imodel];
      }
    }

    if(Milne_Eddington(full_model, stokes, params) != 0){
      *numerical_error = -1;
      return -INFINITY;
    }

    double chisq = 0.0;
    double correction = 0.0;
    size_t sample_count = (size_t)stokes->nw*4U;

    const double *synthetic = stokes->synthetic;
    const double *profile = stokes->profile;
    const double *ivnoi = stokes->inv_noise;

    for(size_t index=0; index<sample_count; index++){
      if(!isfinite(synthetic[index])){
        *numerical_error = -1;
        return -INFINITY;
      }
      double diff = synthetic[index]-profile[index];
      double term = diff*diff*ivnoi[index];
      if(!isfinite(term)){
        *numerical_error = -1;
        return -INFINITY;
      }

      /* Kahan summation retains small contributions when Stokes weights
         differ substantially. */
      double corrected_term = term-correction;
      double updated_chisq = chisq+corrected_term;
      if(!isfinite(updated_chisq)){
        *numerical_error = -1;
        return -INFINITY;
      }
      correction = (updated_chisq-chisq)-corrected_term;
      chisq = updated_chisq;
    }
    
    return -0.5*chisq;
}

/*--------------------------------------------------------------------------------*/
