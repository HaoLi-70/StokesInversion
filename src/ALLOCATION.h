
#ifndef ALLOCATION_H
#define ALLOCATION_H

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <string.h>
#include <stdbool.h>

/*--------------------------------------------------------------------------------*/

enum data_type {enum_int, enum_flt, enum_dbl, enum_cplx, enum_char};

/*--------------------------------------------------------------------------------*/

extern int nrerror(char error_text[]);

extern void *VECTOR(long nl, long nh, enum data_type type, bool Init);

extern int FREE_VECTOR(void *v, long nl, enum data_type type);

extern void *MATRIX(long nrl, long nrh, long ncl, long nch, enum data_type type, \
    bool Init);

extern int FREE_MATRIX(void *m, long nrl, long ncl, enum data_type type);

/*--------------------------------------------------------------------------------*/

extern float **MATRIX_TRI_FLT(long nh, bool Init);

extern int FREE_MATRIX_TRI_FLT(float **m);

extern complex double **MATRIX_TRI_CPLX(long n1, bool Init);

extern int FREE_MATRIX_TRI_CPLX(complex double **Rho);

extern complex double **MATRIX_RHO_CPLX(long nh, bool Init);

extern int FREE_MATRIX_RHO_CPLX(complex double **Rho);

/*--------------------------------------------------------------------------------*/

extern double ***TENSOR_DBL(long Ni0, long Ni1, long Nj0, long Nj1, long Nk0, \
    long Nk1, bool Init);

extern int FREE_TENSOR_DBL(double ***C, long Ni0, long Nj0, long Nk0);

extern complex double ***TENSOR_TRI_CPLX(long n1, long n2, bool Init);

extern int FREE_TENSOR_TRI_CPLX(complex double ***Rho);

extern double ***TENSOR_RHO_DBL(long n1, long n2, bool Init);

extern int FREE_TENSOR_RHO_DBL(double ***Rho);

extern complex double ***TENSOR_RHO_CPLX(long n1, long n2, bool Init);

extern int FREE_TENSOR_RHO_CPLX(complex double ***Rho);

/*--------------------------------------------------------------------------------*/

#endif /* ALLOCATION_H */
