#ifndef MCMC_h

#define MCMC_h


#include <stdio.h>

#include <math.h>

#include <time.h>

#include "ANGULAR_MOMENTUM.h"

#include "ALLOCATION.h"

#include "RANDOM_NUMBER.h"

#include "LIKELIHOOD.h"


double Dream(int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, int Num_Generation, int Num_Chain, int Num_Cr, int Pair_Number, double **Parameter_Bounds, int Indx_Bounds, TRANSITION_INFORMATION ATOM, _Bool Cr_adapt_flag, double **RESULT_DREAM, double **LIKELIHOOD_DREAM, double *Best_fit, int *Convergence);

int Dream_Sample(double **RESULT_DREAM, int Num_Chain, int Num_Parameter, int Num_pair, int Indx_Chain, int Indx_Generation, int Num_Cr, double *Cr, double *Cr_Prob, double **Parameter_Bounds, int Indx_Bounds, long *idum, double *Dream_Sample);

int Sample_PairNum(int Pair_Number, long *idum);

void Ini_Cr(int Num_Cr, double *Cr, double *Cr_Prob, double *Jump_Dis, int *Cr_ups);

int Sample_Cr(int Num_Cr, double *Cr_Prob, long *idum);

void Cr_dis_Update(double **RESULT_DREAM, int Num_Chain, int Num_Parameter, int Indx_Chain, int Indx_Generation, int Indx_Cr, int *Cr_L, double *Cr_Delta);

void Cr_Pro_Update(int Num_Cr, double *Cr_Delta, int *Cr_L, double *Cr_Prob);

void Dream_Diff(double **RESULT_DREAM, int Indx_Generation, int Num_Chain, int Num_Parameter, int Num_Pair, int Num_Jump, int *Jump_dim, long *idum, double *Diff);

void Dream_Dim(int Num_Parameter, double *Cr, int Indx_Cr, long *idum, int *Num_Jump, int *Jump_dim);

void Bounds_Enforce(double *Parameter_Sample, int Num_Parameter, double **Parameter_Bounds, int Index_case);

void PSRF(double **RESULT_MCMC, int Indx_Generation, int Num_Chain, int Num_Parameter, double *R);

int Rm_Outlierchain(int Num_Chain, int Indx_Generation, int Num_parameter, double **RESULT_MCMC, double **likelihood);

void HPSORT(long n, double *ra);

void MCMC_STD(int Begin_Generation, int End_Generation, int Num_Chain, int Num_Parameter, double **RESULT_MCMC, double *STD, double *Mean);

void Chain_Ini(int Num_Chain, int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, double **Parameter_Bounds, TRANSITION_INFORMATION ATOM, double **RESULT_MCMC,  double *likelihood);

int GEMC(int Num_Chain_Gemc, int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, double **Parameter_Bounds, TRANSITION_INFORMATION ATOM, double **GEMC_RESULT, double *LIKELIHOOD, int Indx_Bounds);


#endif /* MCMC_h */
