#include "LIKELIHOOD.h"

double Likelihood_Log(double *Parameter_Sample, int Num_Wav, double *Wavelength_Data, double **Profile_Data, TRANSITION_INFORMATION Atom){
    
    /***********************************************************************************
     Purpose:
     Calculate the log-likelihood of Stokes profiles.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     Parameter_Sample[9], parameters sampled.
     Num_Wav, number of profile points.
     Wavelength_Data[Num_Wav], wavelength of each profile points.
     Profile_Data[Num_Wav][4], observed Stokes profiles.
     Atom, atomic information.
     Output parameters:
     Likelihood, return the log-likelihood.
     ***********************************************************************************/
    
    int i, k;
    
    double Likelihood=0;
    
    double Stokes[4], Weight[4];
    
    Weight[0]=1;
    
    Weight[1]=9;
    
    Weight[2]=9;
    
    Weight[3]=4;
    
    for (i=0; i<Num_Wav; i++) {
        
        Stokes_profile_Normal(Parameter_Sample[0], Parameter_Sample[1], Parameter_Sample[2], Parameter_Sample[3], Parameter_Sample[4], Parameter_Sample[5], Parameter_Sample[6], Parameter_Sample[7], Parameter_Sample[8], Wavelength_Data[i], Atom, Stokes);
        
        for (k=0; k<4; k++)
            
            Likelihood += (Profile_Data[i][k]-Stokes[k])*(Profile_Data[i][k]-Stokes[k])/Profile_Data[i][0]/Profile_Data[i][0]*Weight[k];
            
    }
    
    return -Likelihood/(Num_Wav*4-9.)/2*1e8;//25;
    
}
