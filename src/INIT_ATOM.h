
#ifndef INIT_ATOM_H
#define INIT_ATOM_H

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Atom{
  double Lambda;
  double Geffect;
  double BShift;
  double Ju, Lu, Su, Gu;
  double Jl, Ll, Sl, Gl;
}STRUCT_ATOM;

/*--------------------------------------------------------------------------------*/

extern double Geffect(double Gu, double Gl, double Ju, double Jl);

extern double Gfactor(double J, double L, double S);

extern int Init_Atom(STRUCT_ATOM *Atom);

/*--------------------------------------------------------------------------------*/

#endif /* INIT_ATOM_H */
