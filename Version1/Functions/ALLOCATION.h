#ifndef ALLOCATION_H

#define ALLOCATION_H


#include <stdio.h>

#include <stdlib.h>

#include <stddef.h>

#include <complex.h>


extern void nrerror(char error_text[]);

extern int *VECTOR_INT(long nl, long nh);

extern void FREE_VECTOR_INT(int *v, long nl);

extern float *VECTOR_FLOAT(long nl, long nh);

extern void FREE_VECTOR_FLOAT(float *v, long nl);

extern double *VECTOR_DOUBLE(long nl, long nh);

extern void FREE_VECTOR_DOUBLE(double *v, long nl);

extern complex double *VECTOR_COMPLEX(long nl, long nh);

extern void FREE_VECTOR_COMPLEX(complex double *v, long nl);

extern char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch);

extern void FREE_MATRIX_CHAR(char **m, long nrl, long ncl);

extern int **MATRIX_INT(long nrl, long nrh, long ncl, long nch);

extern void FREE_MATRIX_INT(int **m, long nrl, long ncl);

extern float **MATRIX_FLOAT(long nrl, long nrh, long ncl, long nch);

extern void FREE_MATRIX_FLOAT(float **m, long nrl, long ncl);

extern double **MATRIX_DOUBLE(long nrl, long nrh, long ncl, long nch);

extern void FREE_MATRIX_DOUBLE(double **m, long nrl, long ncl);

extern complex double **MATRIX_COMPLEX(long nrl, long nrh, long ncl, long nch);

extern void FREE_MATRIX_COMPLEX(complex double **m, long nrl, long ncl);

extern float **MATRIX_DIAGONAL(long nh);

extern void FREE_MATRIX_DIAGONAL(float **m);

extern complex double **MATRIX_RHO(long nh);

extern void FREE_MATRIX_RHO(complex double **Rho);


#endif /* ALLOCATION_H */
