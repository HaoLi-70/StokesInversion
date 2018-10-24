#include "ALLOCATION.h"


extern void nrerror(char error_text[]){
    
    /***********************************************************************************
     Purpose:
     Numerical Recipes standard error handler.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     error_text[], the error text.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    fprintf(stderr,"Numerical Recipes run-time error...\n");
    
    fprintf(stderr,"%s\n",error_text);
    
    fprintf(stderr,"...now exiting to system...\n");
    
    exit(1);
    
}


extern int *VECTOR_INT(long nl, long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate an int vector with subscript range v[nl..nh].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    int *v;
    
    v=(int *)malloc(((nh-nl+1)*sizeof(int)));
    
    int i;
    
    for (i=0; i<=nh-nl; i++)
        v[i]=0;
    
    if (!v)
        
        nrerror("allocation failure in VECTOR_INT()");
    
    return v-nl;
    
}


extern void FREE_VECTOR_INT(int *v, long nl){
    
    /***********************************************************************************
     Purpose:
     free an int vector allocated with VECTOR_INT().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((v+nl));
    
    return;
    
}


extern float *VECTOR_FLOAT(long nl, long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate a float vector with subscript range v[nl..nh].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    float *v;
    
    v=(float *)malloc(((nh-nl+1)*sizeof(float)));
    
    int i;
    
    for (i=0; i<=nh-nl; i++)
        v[i]=0;
    
    if (!v)
        
        nrerror("allocation failure in VECTOR_FLOAT()");
    
    return v-nl;
    
}


extern void FREE_VECTOR_FLOAT(float *v, long nl){
    
    /***********************************************************************************
     Purpose:
     free a float vector allocated with VECTOR_FLOAT().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    free((v+nl));
    
    return;
    
}


extern double *VECTOR_DOUBLE(long nl, long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate a double vector with subscript range v[nl..nh].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    double *v;
    
    v=(double *)malloc(((nh-nl+1)*sizeof(double)));
    
    int i;
    
    for (i=0; i<=nh-nl; i++)
        
        v[i]=0;

    if (!v)
        
        nrerror("allocation failure in VECTOR_DOUBLE()");
    
    return v-nl;
    
}


extern void FREE_VECTOR_DOUBLE(double *v, long nl){
    
    /***********************************************************************************
     Purpose:
     free a double vector allocated with VECTOR_DOUBLE().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    free((v+nl));
    
    return;
    
}


extern complex double *VECTOR_COMPLEX(long nl, long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate a complex double vector with subscript range v[nl..nh].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    complex double *v;
    
    v=(complex double *)malloc((nh-nl+1)*sizeof(complex double));
    
    int i;
    
    for (i=0; i<=nh-nl; i++)
        
        v[i]=0;
    
    if (!v)
        
        nrerror("allocation failure in VECTOR_COMPLEX()");
    
    return v-nl;
    
}


extern void FREE_VECTOR_COMPLEX(complex double *v, long nl){
    
    /***********************************************************************************
     Purpose:
     free a complex double vector allocated with VECTOR_COMPLEX().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((v+nl));
    
    return;
    
}


extern char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch){
    
    /***********************************************************************************
     Purpose:
     Allocate an char matrix with subscript range m[nrl..nrh][ncl..nch].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     30 Sept. 2018, Hao Li
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    
    char **m;
    
    m=(char **) malloc(((nrow)*sizeof(char*)));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_CHAR()");
    
    m -= nrl;
    
    m[nrl]=(char *) malloc(((nrow*ncol)*sizeof(char)));
    
    for (i=0; i<(nrow*ncol); i++)
       
        m[nrl][i]=0;
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_CHAR()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++)
        
        m[i]=m[i-1]+ncol;
    
    return m;
    
}


extern void FREE_MATRIX_CHAR(char **m, long nrl, long ncl){
    
    /***********************************************************************************
     Purpose:
     free a char matrix allocated with MATRIX_CHAR().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
    
}


extern int **MATRIX_INT(long nrl, long nrh, long ncl, long nch){
    
    /***********************************************************************************
     Purpose:
     Allocate an int matrix with subscript range m[nrl..nrh][ncl..nch].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    
    int **m;
    
    m=(int **) malloc(((nrow)*sizeof(int*)));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_INT()");
    
    m -= nrl;
    
    m[nrl]=(int *) malloc(((nrow*ncol)*sizeof(int)));
    
    for (i=0; i<(nrow*ncol); i++)
        
        m[nrl][i]=0;
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_INT()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++)
        
        m[i]=m[i-1]+ncol;
    
    return m;
    
}


extern void FREE_MATRIX_INT(int **m, long nrl, long ncl){
    
    /***********************************************************************************
     Purpose:
     free an int matrix allocated with MATRIX_INT().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
    
}


extern float **MATRIX_FLOAT(long nrl, long nrh, long ncl, long nch){
    
    /***********************************************************************************
     Purpose:
     Allocate a float matrix with subscript range m[nrl..nrh][ncl..nch].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    
    float **m;
    
    /* allocate pointers to rows */
    m=(float **) malloc(((nrow)*sizeof(float*)));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_FLOAT()");
    
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(float *) malloc(((nrow*ncol)*sizeof(float)));
    
    for (i=0; i<(nrow*ncol); i++)
        
        m[nrl][i]=0;
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_FLOAT()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++)
        
        m[i]=m[i-1]+ncol;
    
    /* return pointer to array of pointers to rows */
    
    return m;
    
}


extern void FREE_MATRIX_FLOAT(float **m, long nrl, long ncl){
    
    /***********************************************************************************
     Purpose:
     free a float matrix allocated with MATRIX_FLOAT().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
    
}


extern double **MATRIX_DOUBLE(long nrl, long nrh, long ncl, long nch){
    
    /***********************************************************************************
     Purpose:
     Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    
    double **m;
    
    /* allocate pointers to rows */
    m=(double **) malloc(((nrow)*sizeof(double*)));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_DOUBLE()");
    
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc(((nrow*ncol)*sizeof(double)));
    
    for (i=0; i<(nrow*ncol); i++)
        
        m[nrl][i]=0;
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_DOUBLE()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++)
        
        m[i]=m[i-1]+ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
    
}


extern void FREE_MATRIX_DOUBLE(double **m, long nrl, long ncl){
    
    /***********************************************************************************
     Purpose:
     free a double matrix allocated with MATRIX_DOUBLE().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
    
}


extern complex double **MATRIX_COMPLEX(long nrl, long nrh, long ncl, long nch){
    
    /***********************************************************************************
     Purpose:
     Allocate a complex double matrix with subscript range m[nrl..nrh][ncl..nch].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    
    complex double **m;
    
    /* allocate pointers to rows */
    m=(complex double **) malloc(((nrow)*sizeof(complex double*)));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_COMPLEX()");
    
    m -= nrl;
    
    /* allocate rows and set pointers to them */
    m[nrl]=(complex double *) malloc(((nrow*ncol)*sizeof(complex double)));
    
    for (i=0; i<(nrow*nrow); i++)
        
        m[nrl][i]=0;
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_COMPLEX()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1;i<=nrh;i++)
        
        m[i]=m[i-1]+ncol;
    
    /* return pointer to array of pointers to rows */
    return m;
    
}


extern void FREE_MATRIX_COMPLEX(complex double **m, long nrl, long ncl){
    
    /***********************************************************************************
     Purpose:
     free a complex double matrix allocated with MATRIX_COMPLEX().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
    
}


extern float **MATRIX_DIAGONAL(long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate a float diagonal matrix with subscript range m[0..nh][0..(nh)].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nh, the dimension of the diagonal matrix..
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    long nrow=nh+1,i;
    
    float **m;
    
    /* allocate pointers to rows */
    m=(float **) malloc(nrow*sizeof(float*));
    
    if (!m)
        
        nrerror("allocation failure 1 in MATRIX_DIAGONAL()");
    
    /* allocate rows and set pointers to them */
    m[0]=(float *) malloc((nh+2)*nrow/2*sizeof(float));
    
    for (i=0; i<(nh+2)*nrow/2; i++)
        
        m[0][i]=0;
    
    if (!m[0])
        
        nrerror("allocation failure 2 in MATRIX_DIAGONAL()");
    
    for(i=1;i<=nh;i++)
        
        m[i]=m[i-1]+i;
    
    /* return pointer to array of pointers to rows */
    return m;
    
}


extern void FREE_MATRIX_DIAGONAL(float **m){
    
    /***********************************************************************************
     Purpose:
     free an float diagonal matrix allocated with MATRIX_DIAGONAL().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **m, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    free((m[0]));
    
    free((m));
    
    return;
    
}


extern complex double **MATRIX_RHO(long nh){
    
    /***********************************************************************************
     Purpose:
     Allocate an complex double rho matrix with subscript range m[0..nh][0..(nh*2+1)].
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     nh, the dimension of the rho matrix.
     Output parameters:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    long nrow=nh+1,i;
    
    complex double **R;
    
    /* allocate pointers to rows */
    R=(complex double **)malloc(nrow*sizeof(complex double *));
    
    if (!R)
        
        nrerror("allocation failure 1 in MATRIX_RHO()");
    
    R[0]=(complex double *)malloc(nrow*nrow*sizeof(complex double));
    
    for (i=0; i<(nrow*nrow); i++)
        
        R[0][i]=0;
    
    if (!R[0])
        
        nrerror("allocation failure 2 in MATRIX_RHO()");
    
    for (i=1; i<=nh; i++)
        
        R[i]=R[i-1]+2*i;
    
    /* return pointer to array of pointers to rows */
    return R;
    
}


extern void FREE_MATRIX_RHO(complex double **Rho){
    
    /***********************************************************************************
     Purpose:
     free a complex rho matrix allocated by MATRIX_RHO().
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     **Rho, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/

    free((Rho[0]));
    
    free((Rho));
    
    return;
    
}
