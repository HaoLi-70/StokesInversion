
#include "sorting.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

typedef struct{
    double value;
    int index;
}SORT_ITEM;

/*--------------------------------------------------------------------------------*/

static int COMPARE_VALUES(const void *lhs, const void *rhs){

    /*######################################################################
      Purpose:
        Compare two double values in ascending order for qsort.
      Input parameters:
        lhs, pointer to the left value.
        rhs, pointer to the right value.
      Return:
        -1 if lhs is smaller, 1 if lhs is larger, otherwise 0.
     ######################################################################*/ 

    const double left = *(const double *)lhs;
    const double right = *(const double *)rhs;

    if(left < right) return -1;
    if(left > right) return 1;
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int COMPARE_ITEMS(const void *lhs, const void *rhs){

    /*######################################################################
      Purpose:
        Compare two indexed values in ascending order for qsort. Use the
          original index to obtain a deterministic order for equal values.
      Input parameters:
        lhs, pointer to the left SORT_ITEM.
        rhs, pointer to the right SORT_ITEM.
      Return:
        -1 if lhs precedes rhs, 1 if rhs precedes lhs, otherwise 0.
     ######################################################################*/

    const SORT_ITEM *left = (const SORT_ITEM *)lhs;
    const SORT_ITEM *right = (const SORT_ITEM *)rhs;

    if(left->value < right->value) return -1;
    if(left->value > right->value) return 1;
    if(left->index < right->index) return -1;
    if(left->index > right->index) return 1;
    return 0;
}

/*--------------------------------------------------------------------------------*/

int SORT_VALUES(double *values, size_t count){

    /*######################################################################
      Purpose:
        Sort a zero-based array of double values in ascending order.
      Input parameters:
        count, number of values in the array.
        values, the array to sort.
      Output parameters:
        values, replaced by its sorted rearrangement.
      Return:
        0 on success; -1 for an invalid array.
      Note:
        Input values must not contain NaN.
     ######################################################################*/

    if(count == 0) return 0;
    if(!values || count > SIZE_MAX/sizeof(*values)) return -1;
    if(count > 1) qsort(values, count, sizeof(*values), COMPARE_VALUES);
    return 0;
}

/*--------------------------------------------------------------------------------*/

int SORT_INDICES(const double *values, size_t count, int *indices){

    /*######################################################################
      Purpose:
        Generate ascending-order indices without changing the input values.
      Input parameters:
        count, number of values in the arrays.
        values, input values to index.
      Output parameters:
        indices, zero-based indices such that values[indices[i]] is ordered.
      Return:
        0 on success; -1 for invalid input or an allocation failure.
      Note:
        Equal values retain their original index order. Input values must not
          contain NaN.
     ######################################################################*/

    if(count == 0) return 0;
    if(!values || !indices || count > (size_t)INT_MAX 
        || count > SIZE_MAX/sizeof(SORT_ITEM)) return -1;
    SORT_ITEM *items = (SORT_ITEM *)malloc(count*sizeof(*items));
    if(!items) return -1;

    for(size_t i=0; i<count; i++){
      items[i].value = values[i];
      items[i].index = (int)i;
    }

    if(count > 1) qsort(items, count, sizeof(*items), COMPARE_ITEMS);

    for(size_t i=0; i<count; i++){
      indices[i] = items[i].index;
    }

    free(items);
    return 0;
}

/*--------------------------------------------------------------------------------*/
