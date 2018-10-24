#ifndef INVERSION_h

#define INVERSION_h


#include "INIT.h"

#include "MCMC.h"

#include "LIKELIHOOD.h"

#include "PATH.h"


extern int Inversion(double **Profile_Data, double *Wavelength_Data, int Num_Wav, int Flag_Correlation, int Flag_Histdata, double **Result, double *Accept_Rates);


#endif /* INVERSION_h */
