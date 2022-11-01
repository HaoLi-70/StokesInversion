
#ifndef LIKELIHOOD_h
#define LIKELIHOOD_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include "ALLOCATION.h"
#include "FORWARD_ME.h"
#include "READ_INPUT_MCMC.h"

/*--------------------------------------------------------------------------------*/

double Likelihood_Log(STRUCT_ATOM *Atom, double *Par, \
    STRUCT_OBSERVATION *Observation);

/*--------------------------------------------------------------------------------*/

#endif /* LIKELIHOOD_h */
