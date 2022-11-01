
#include "SORT.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:
        30 Oct. 2022.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern void HPSORT(long n, double *ra){
    
    /*######################################################################
      Purpose:
        sorts an array ra[1..n] into ascending numerical order using the
            Heapsort algorithm.
       Record of revisions:
        8 Sept. 2021.
       Input parameters:
        n, the number of elements the array.
        ra[], the array input.
      Output parameters:
        ra, replaced by its sorted rearrangement
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

   
    
    if(n < 2) return;
    long i,ir,j,l;
    double rra; 
    l = (n >> 1)+1;
    ir = n;

    for(;;){
      if(l>1){
        rra = ra[--l];
      }else{
        rra = ra[ir];
        ra[ir] = ra[1];
        if(--ir == 1){
          ra[1]=rra;
          break;
        }
      }
        
      i = l;
      j = l+l;
        
      while(j <= ir){
        if(j < ir && ra[j] < ra[j+1]) j++;
        if(rra < ra[j]){
          ra[i] = ra[j];
          i = j;
          j <<= 1;
        }else break;
      }
      ra[i] = rra;
    }
    return;
}

/*--------------------------------------------------------------------------------*/

extern void qsort_index(int m, int n, double *arr, int *indx){
    
    /*######################################################################
      Purpose:
        sorts an array ra[m..n] into ascending numerical order using the
            quick sort algorithm.
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        m, n, then indexes of elements in the array.
        arr, the array input.
      Output parameters:
        indx, arr[indx[j]] is in ascending order.
      reference:
        numerical recipes in C 2ed.
            Indexes an array arr[m..n], i.e., outputs the array indx[m..n]
                such that arr[indx[j]] is in ascending order The input
                quantities arr are not changed.
     ######################################################################*/

    if(n - m < 2) return;
    
    const int  M = 7;
    const int NSTACK = 100;
    int i, indxt, ir = n, j, k, l = m;
    int jstack=0;
    float a;

    int *istack = (int *)VECTOR(1, NSTACK, enum_int, false);

    for(j=l;j<=n;j++) indx[j] = j;
    
    for(;;){
      if(ir-l < M){
        for(j=l+1;j<=ir;j++){
          indxt = indx[j];
          a = arr[indxt];
          for (i=j-1;i>=l;i--) {
            if(arr[indx[i]] <= a) break;
            indx[i+1] = indx[i];
          }
          indx[i+1] = indxt;
        }
        if(jstack == 0) break;
        ir = istack[jstack--];
        l = istack[jstack--];
      }else{
        k = (l+ir) >> 1;
        M_SWAP(indx[k],indx[l+1]);
        if(arr[indx[l]] > arr[indx[ir]]) M_SWAP(indx[l],indx[ir]);
        if(arr[indx[l+1]] > arr[indx[ir]]) M_SWAP(indx[l+1],indx[ir]);
        if(arr[indx[l]] > arr[indx[l+1]]) M_SWAP(indx[l],indx[l+1]);
        i = l+1;
        j = ir;
        indxt = indx[l+1];
        a = arr[indxt];
        for(;;){
          do i++; while(arr[indx[i]] < a);
          do j--; while(arr[indx[j]] > a);
          if(j < i)break;
          M_SWAP(indx[i],indx[j]);
        }
        indx[l+1] = indx[j];
        indx[j] = indxt;
        jstack += 2;
        if(jstack > NSTACK) nrerror("NSTACK too small in indexx.");
        if(ir-i+1 >= j-l){
          istack[jstack] = ir;
          istack[jstack-1] = i;
          ir = j-1;
        }else{
          istack[jstack] = j-1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
    FREE_VECTOR(istack, 1, enum_int);
  
    return;
}
    
/*--------------------------------------------------------------------------------*/

extern void quick_sort(int m, int n, double *array){
    
    /*######################################################################
      Purpose:
        sorts an array ra[m..n] into ascending numerical order using the
            quick sort algorithm.
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        m, n, then indexes of elements in the array.
        array, the array input.
      Output parameters:
        array, replaced by its sorted rearrangement
      reference:
        numerical recipes in C 2ed.
            Sorts an array arr[1..n] into ascending numerical order using
                the Quicksort algorithm. arr is replaced on output by its
                sorted rearrangement.
     ######################################################################*/

    if(n - m < 2) return;

    const int M = 7;
    const int NSTACK = 100;
    int i, ir = n, j, k, l = m;
    int jstack = 0;
    double a;
    int *istack = (int *)VECTOR(1, NSTACK, enum_int, false);

    for(;;){ 
      if(ir-l < M){
        for(j=l+1;j<=ir;j++){
          a = array[j];
          for(i=j-1;i>=l;i--){
            if(array[i] <= a) break;
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if(jstack == 0) break;
        ir = istack[jstack--];
        l = istack[jstack--];
      }else{
        k = (l+ir) >> 1;
        M_SWAP(array[k],array[l+1]);
        if(array[l] > array[ir]) M_SWAP(array[l],array[ir]);
        if(array[l+1] > array[ir]) M_SWAP(array[l+1],array[ir]);
        if(array[l] > array[l+1]) M_SWAP(array[l],array[l+1]);
        i = l+1;
        j = ir;
        a = array[l+1];
        for(;;){
          do i++; while(array[i] < a);
          do j--; while(array[j] > a);
          if(j < i) break;
          M_SWAP(array[i],array[j]);
        }
        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if(jstack > NSTACK) nrerror("NSTACK too small in sort.");
        if(ir-i+1 >= j-l){
          istack[jstack] = ir;
          istack[jstack-1] = i;
          ir=j-1;
        }else{
          istack[jstack] = j-1;
          istack[jstack-1] = l;
          l=i;
        }
      }
    }
    FREE_VECTOR(istack, 1, enum_int);
    return;;
}

/*--------------------------------------------------------------------------------*/


