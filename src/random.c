
#include "random.h"

#include <math.h>
#include <stddef.h>
#include <time.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Independent MPI-rank random states with uniform and Gaussian
        generators backed by GSL or the Numerical Recipes implementation.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

static uint64_t RNG_MIX(uint64_t value){

    /*######################################################################
      Purpose:
        Mix a 64-bit seed so nearby input values produce unrelated states.
      Input parameters:
        value, raw seed material.
      Return:
        Mixed 64-bit seed value.
      Note:
        SplitMix64 finalizer.
    ######################################################################*/

    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value^(value >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
    value = (value^(value >> 27))*UINT64_C(0x94d049bb133111eb);

    return value^(value >> 31);
}

/*--------------------------------------------------------------------------------*/

int RNG_INIT(STRUCT_RNGState *state, int rank){

    /*######################################################################
      Purpose:
        Initialize a rank-specific random-number-generator state.
      Input parameters:
        state, random-number-generator state to initialize.
        rank, MPI rank used to diversify the seed.
      Return:
        0 on success; -1 for a null state or GSL allocation failure.
      Note:
        Numerical Recipes in C, 2nd edition, or GSL.
    ######################################################################*/

    if(!state) return -1;

    struct timespec now = {0, 0};
    if(timespec_get(&now, TIME_UTC) != TIME_UTC){
      now.tv_sec = time(NULL);
    }
    uint64_t seed = (uint64_t)now.tv_sec;
    seed ^= (uint64_t)now.tv_nsec << 32;
    seed ^= (uint64_t)(uint32_t)rank << 1;
    seed = RNG_MIX(seed);

#ifdef USE_GSL
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    state->gsl = gsl_rng_alloc(T);
    if(!state->gsl) return -1;
    gsl_rng_set(state->gsl, (unsigned long)seed);
#else
    uint32_t initial = (uint32_t)(seed%UINT32_C(2147483646))+1U;
    state->idum = -(int32_t)initial;
    state->iy = 0;
    state->iset = false;
    state->gset = 0.0;
#endif

    return 0;
}

/*--------------------------------------------------------------------------------*/

double RNG_UNIFORM(STRUCT_RNGState *state){

    /*######################################################################
      Purpose:
        Generate a uniformly distributed random number in the interval 
            [0, 1) using either GSL or the method described 
            in Numerical Recipes.
      Input parameters:
        state, initialized random-number-generator state.
      Return:
        A uniformly distributed value in [0, 1), or NAN for an invalid state.
      Note:
        Numerical Recipes in C 2nd edition. 
    ######################################################################*/

#ifdef USE_GSL
    if(!state || !state->gsl) return NAN;

    return gsl_rng_uniform(state->gsl); 
#else

    if(!state) return NAN;

    static const int IA = 16807;
    static const int IM = 2147483647;
    static const double AM = (1.0/2147483647);
    static const int IQ = 127773;
    static const int IR = 2836;
    static const int NDIV = (1+(2147483647-1)/NTAB);
    static const double RNMX = (1.0-1.2e-7); 

    int j;
    int32_t k;
    double temp;

    if(state->idum<=0 || !state->iy){
      if (-(state->idum) < 1)
        state->idum = 1;
      else
        state->idum = -(state->idum);

      for(j=NTAB+7; j>= 0; j--){
        k = state->idum/IQ;
        state->idum = IA*(state->idum-k*IQ)-IR*k;
        if (state->idum<0)
          state->idum += IM;
        if (j<NTAB)
          state->iv[j] = state->idum;
      }
      state->iy = state->iv[0];
    }

    k = state->idum/IQ;
    state->idum = IA*(state->idum-k*IQ)-IR*k;
    if(state->idum<0)
      state->idum += IM;

    j = state->iy/NDIV;
    state->iy = state->iv[j];
    state->iv[j] = state->idum;

    temp = AM*state->iy;

    return (temp>RNMX) ? RNMX : temp;
#endif
}

/*--------------------------------------------------------------------------------*/

double RNG_GAUSS(STRUCT_RNGState *state){

    /*######################################################################
      Purpose:
        Generate a Gaussian-distributed random number with unit width using 
            either GSL or the method described in Numerical Recipes.
      Input parameters:
        state, initialized random-number-generator state.
      Return:
        A normally distributed value with mean 0 and standard deviation 1,
        or NAN for an invalid state.
      Note:
        Numerical Recipes in C 2nd edition / GSL. 
    ######################################################################*/

#ifdef USE_GSL
    if(!state || !state->gsl) return NAN;
    
    return gsl_ran_gaussian(state->gsl, 1.0); 
#else
    if(!state) return NAN;
    double v1,v2,rsq,fac;
    if(!state->iset){
      do{
        v1 = 2.0*RNG_UNIFORM(state)-1.0;
        v2 = 2.0*RNG_UNIFORM(state)-1.0;
        rsq = v1*v1+v2*v2;
      }while(rsq>= 1.0 || rsq==0.0);
      fac = sqrt(-2.0*log(rsq)/rsq);
      state->gset = v1*fac;
      
      state->iset = true;
      return v2*fac;
    }else{
      state->iset = false;

      return state->gset;
    }
#endif
}

/*--------------------------------------------------------------------------------*/
