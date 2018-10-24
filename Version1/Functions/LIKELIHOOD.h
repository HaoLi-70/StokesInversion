#ifndef LIKELIHOOD_h

#define LIKELIHOOD_h


#include <stdio.h>

#include "ANGULAR_MOMENTUM.h"

#include "ME_PROFILE.h"


double Likelihood_Log(double *Parameter_Sample, int Num_Wav, double *Wavelength_Data, double **Profile_Data, TRANSITION_INFORMATION Atom);


#endif /* LIKELIHOOD_h */
