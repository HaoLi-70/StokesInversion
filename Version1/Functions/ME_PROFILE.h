#ifndef ME_PROFILE_h

#define ME_PROFILE_h


#include <stdio.h>

#include "FADDEEVA.h"

#include "ALLOCATION.h"

#include "ANGULAR_MOMENTUM.h"


void Stokes_profile_Normal(double Bmag, double Theta, double Chi, double Velocity, double Doppler_width, double Source_function, double Coeffi_Abs, double Beta, double a, double Lambda, TRANSITION_INFORMATION Atom, double *Stokes);

void Stokes_profile(double Bmag, double Theta, double Chi, double Velocity, double Doppler_width, double Source_function, double Coeffi_Abs, double Beta, double a, double Lambda, TRANSITION_INFORMATION Atom, double *Stokes);


#endif /* ME_PROFILE_h */
