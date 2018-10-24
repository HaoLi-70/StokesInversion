#ifndef INIT_h

#define INIT_h


#include <stdio.h>

#include <math.h>

#include <string.h>

#include "PATH.h"

#include "ALLOCATION.h"

#include "ANGULAR_MOMENTUM.h"


void Init_Wav(int *Wav_Num);

void Init_Atom(TRANSITION_INFORMATION *Atom);

void Init_file(int *Wav_begin, int *Wav_end, int *Pixel_begin, int *Pixel_end, int *Num_file);

void filelist(int Num_file, char **Name_fits, char **Res_fits);

void Init_Wavelength(double *Wavelength_Data, int Num_Wav);

void Init_Parameter(int *Num_Parameter, int *Num_Chain_Gemc, int *Num_Chain_Dream, int *Begin_Generation, int *Num_Generation, int *Num_Cr, int *Pair_Number, int *Indx_Bounds_Gemc, int *Indx_Bounds_Dream, _Bool *Cr_adapt_flag);

void Init_Bounds(double **Profile_Data, double *Wavelength_Data, int Num_Wav, int Num_Parameter, float Lambda, double **Parameter_Bounds);


#endif /* INIT_h */
