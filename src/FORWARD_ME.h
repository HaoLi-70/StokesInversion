
#ifndef FORWARD_ME_h
#define FORWARD_ME_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>

#include "ALLOCATION.h"
#include "FADDEEVA.h"
#include "INIT_ATOM.h"

/*--------------------------------------------------------------------------------*/

enum bounds_type {enum_fold, enum_reflect, enum_set};

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Parameter{
    int Num_Par;
    double **Par_Bounds;
    double *Distribution;
    enum bounds_type *type;
    bool *narrower_guess;
}STRUCT_PAR;


typedef struct Struct_Observation{
    int Num_Lam;
    double *Lambda;
    double **Profile_Data;
    double **Sigma;
    double **Stokes;
    bool weights_flag;
}STRUCT_OBSERVATION;

/*--------------------------------------------------------------------------------*/

extern int Forward_ME(STRUCT_ATOM *Atom, double *Par, \
    STRUCT_OBSERVATION *Observation);

/*--------------------------------------------------------------------------------*/

#endif /* Forward_ME_h */
