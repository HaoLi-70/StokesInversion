
#ifndef FADDEEVA_H
#define FADDEEVA_H

/*--------------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "PARAMETER.h"

/*--------------------------------------------------------------------------------*/

extern int Faddeeva(double Nu, double y, double *H, double *L, int digits);

extern int Faddeeva916(double Nu, double y, double *H, double *L, int digits);

/*--------------------------------------------------------------------------------*/

#endif /* FADDEEVA_H */
