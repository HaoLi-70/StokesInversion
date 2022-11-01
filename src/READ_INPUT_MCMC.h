
#ifndef READ_INPUT_FBLINV_h
#define READ_INPUT_FBLINV_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>

#include "ALLOCATION.h"
#include "STR.h"
#include "ERROR.h"
#include "MPI_CONTROL.h"
#include "FORWARD_ME.h"
#include "INIT_ATOM.h"

/*--------------------------------------------------------------------------------*/

enum keywordtype {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Keywords{
    char keyword[Key_Length];
    char line[Max_Line_Length];
    bool Set, Required;
}STRUCT_KEYS;

typedef struct Struct_Input{
    bool Verbose;
    char Path2Profile[Max_Line_Length];
    bool GEMC_Simulation, Cr_Update, Distribution_Update, Sigma, Output_Samples, \
        Debug;
    int Gemc_Nchains, Dream_Nchains;
    int Num_Gener, MaxNum_Pair;
    int Num_Cr;
}STRUCT_INPUT;

/*--------------------------------------------------------------------------------*/

extern void Keywords_Conversion(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input, \
    STRUCT_MPI *Mpi, STRUCT_ATOM *Atom);

extern void RDINPUT(char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_ATOM *Atom, STRUCT_OBSERVATION *Observation);

extern int RDOBSERVATION(char Filename[], STRUCT_OBSERVATION *OBSERVATION);

/*--------------------------------------------------------------------------------*/

#endif /* READ_INPUT_MCMC_h */
