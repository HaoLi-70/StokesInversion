
#ifndef SORT_h
#define SORT_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include "ALLOCATION.h"

/*--------------------------------------------------------------------------------*/

/*
#define M_SWAP(x,y) do{  \
    (void) (&x==&y);     \
    typeof(x) _x;        \
    _x=(x);              \
    x=(y);               \
    y=_x;                \
}while(0);
*/

#define M_SWAP(x,y) ({    \
    (void) (&x==&y);     \
    typeof(x) _x;        \
    _x=(x);              \
    x=(y);               \
    y=_x;                \
})

//#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

extern void HPSORT(long n, double *ra);

extern void qsort_index(int m, int n, double *arr, int *indx);

extern void quick_sort(int m, int n, double *array);

/*--------------------------------------------------------------------------------*/

#endif /* SORT_h */
