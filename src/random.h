
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>
#include <stdint.h>

#ifdef USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

/*--------------------------------------------------------------------------------*/

#define NTAB 32

typedef struct Struct_RNG_State{
#ifdef USE_GSL
    gsl_rng *gsl;
#else
    int32_t idum;
    int32_t iy;
    int32_t iv[NTAB];

    bool iset;
    double gset; 
#endif
}STRUCT_RNGState;

/*--------------------------------------------------------------------------------*/

extern int RNG_INIT(STRUCT_RNGState *state, int rank);

extern double RNG_UNIFORM(STRUCT_RNGState *state);

extern double RNG_GAUSS(STRUCT_RNGState *state);

/*--------------------------------------------------------------------------------*/
