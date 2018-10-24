#ifndef FADDEEVA_H

#define FADDEEVA_H


#include <math.h>


extern double DAWSON(double x);

extern double Fadd_erfcx(double a);

extern void Faddeeva_Function_V1(double Nu, double y, double *H, double *L);

extern void Faddeeva_Function(double Nu, double y, double *H, double *L);


#endif /* FADDEEVA_H */
