
#pragma once

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Faddeeva{

    // number of figures and the corresponding index, max number of iteration.
    int precision_digits, precision_idx, max_terms;

    // the smallest value, log(min_value), sqrt(log(min_value)), a, and a*a;
    double min_value, log_min_value, sqrt_log_min_value, a, a_squared;
    
    // precomputed exp(-a^2n^2)
    double *exp_a2n2;

}STRUCT_FADDEEVA;

/*--------------------------------------------------------------------------------*/

extern int Faddeeva_init(STRUCT_FADDEEVA *faddeeva);

extern int Faddeeva(double nu, double y, double *h, double *l, 
    STRUCT_FADDEEVA *faddeeva);

/*--------------------------------------------------------------------------------*/
