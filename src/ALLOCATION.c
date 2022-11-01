
#include "ALLOCATION.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        30 Oct. 2022.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

static int *VECTOR_INT(long nl, long nh, bool Init);

static int FREE_VECTOR_INT(int *v, long nl);

static float *VECTOR_FLT(long nl, long nh, bool Init);

static int FREE_VECTOR_FLT(float *v, long nl);

static double *VECTOR_DBL(long nl, long nh, bool Init);

static int FREE_VECTOR_DBL(double *v, long nl);

static complex double *VECTOR_CPLX(long nl, long nh, bool Init);

static int FREE_VECTOR_CPLX(complex double *v, long nl);

/*--------------------------------------------------------------------------------*/

static char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch, bool Init);

static int FREE_MATRIX_CHAR(char **m, long nrl, long ncl);

static int **MATRIX_INT(long nrl, long nrh, long ncl, long nch, bool Init);

static int FREE_MATRIX_INT(int **m, long nrl, long ncl);

static float **MATRIX_FLT(long nrl, long nrh, long ncl, long nch, bool Init);

static int FREE_MATRIX_FLT(float **m, long nrl, long ncl);

static double **MATRIX_DBL(long nrl, long nrh, long ncl, long nch, bool Init);

static int FREE_MATRIX_DBL(double **m, long nrl, long ncl);

static complex double **MATRIX_CPLX(long nrl, long nrh, long ncl, long nch, \
                                       bool Init);

static int FREE_MATRIX_CPLX(complex double **m, long nrl, long ncl);

/*--------------------------------------------------------------------------------*/

extern int nrerror(char error_text[]){

    /*######################################################################
      Purpose:
        Numerical Recipes standard error handler.
      Record of revisions:
        21 Mar. 2018.
      Input parameters:
        error_text[], the error text.
      Return:
        .
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
    
    return 0 ;
}

/*--------------------------------------------------------------------------------*/

static int *VECTOR_INT(long nl, long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a int vector v[i] with subscript range nl<=i<=nh.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nl, nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long byte = (nh-nl+1)*sizeof(int);
    int *v;

    v = (int *)malloc(byte);
  
    if(!v) nrerror("allocation failure in VECTOR_INT()");

    if(Init) memset(v,0,byte);
    
    return v-nl;
}

/*--------------------------------------------------------------------------------*/

static int FREE_VECTOR_INT(int *v, long nl){
    
    /*######################################################################
      Purpose:
        free an int vector allocated with VECTOR_INT().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        *v, the pointer.
        nl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((v+nl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static float *VECTOR_FLT(long nl, long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a float vector v[i] with subscript range nl<=i<=nh.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nl, nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long byte = (nh-nl+1)*sizeof(float);
    float *v;

    v = (float *)malloc(byte);

    if(!v) nrerror("allocation failure in VECTOR_FLT()");
    
    if(Init) memset(v, 0, byte);
    
    return v-nl;
}

/*--------------------------------------------------------------------------------*/

static int FREE_VECTOR_FLT(float *v, long nl){
    
    /*######################################################################
      Purpose:
        free a float vector allocated with VECTOR_FLT().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        *v, the pointer.
        nl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((v+nl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static double *VECTOR_DBL(long nl, long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a double vector v[i] with subscript range nl<=i<=nh.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nl, nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long byte = (nh-nl+1)*sizeof(double);
    double *v;

    v = (double *)malloc(byte);

    if(!v) nrerror("allocation failure in VECTOR_DBL()");
    
    if(Init) memset(v, 0, byte);
    
    return v-nl;
}

/*--------------------------------------------------------------------------------*/

static int FREE_VECTOR_DBL(double *v, long nl){
    
    /*######################################################################
      Purpose:
        free a double vector allocated with VECTOR_DBL().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        *v, the pointer.
        nl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((v+nl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static complex double *VECTOR_CPLX(long nl, long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a complex vector v[i] with subscript range nl<=i<=nh.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nl, nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long byte = (nh-nl+1)*sizeof(complex double);
    complex double *v;

    v = (complex double *)malloc(byte);
    
    if(!v) nrerror("allocation failure in VECTOR_CPLX()");
    
    if(Init) memset(v, 0, byte);
    
    return v-nl;
}

/*--------------------------------------------------------------------------------*/

static int FREE_VECTOR_CPLX(complex double *v, long nl){
    
    /*######################################################################
      Purpose:
        free a complex vector allocated with VECTOR_CPLX().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        *v, the pointer.
        nl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((v+nl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a char matrix m[i][j] with subscript range
            nrl<=i<=nrh, ncl<=j<=ncl.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    long byte = (nrow*ncol)*sizeof(char);
    char **m;

    m = (char **)malloc(((nrow)*sizeof(char*)));
    
    if(!m) nrerror("allocation failure 1 in MATRIX_CHAR()");
    
    m -= nrl;

    m[nrl] = (char *)malloc(byte);
    
    if(!m[nrl]) nrerror("allocation failure 2 in MATRIX_CHAR()");

    if(Init) memset(m[nrl], 0, byte);
    
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

static int FREE_MATRIX_CHAR(char **m, long nrl, long ncl){
    
    /*######################################################################
      Purpose:
        free a char matrix allocated with MATRIX_CHAR().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
        nrl, ncl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free(m[nrl]+ncl);
    free(m+nrl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int **MATRIX_INT(long nrl, long nrh, long ncl, long nch, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a int matrix m[i][j] with subscript range
            nrl<=i<=nrh, ncl<=j<=ncl.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    long byte = (nrow*ncol)*sizeof(int);
    int **m;

    m = (int **)malloc(((nrow)*sizeof(int*)));
    
    if(!m) nrerror("allocation failure 1 in MATRIX_INT()");
    
    m -= nrl;

    m[nrl] = (int *)malloc(byte);
    
    if(!m[nrl]) nrerror("allocation failure 2 in MATRIX_INT()");
    
    if(Init) memset(m[nrl], 0, byte);
    
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

static int FREE_MATRIX_INT(int **m, long nrl, long ncl){
    
    /*######################################################################
      Purpose:
        free an int matrix allocated with MATRIX_INT().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
        nrl, ncl, the first index of subscript range.
      Return:
        .
     ######################################################################*/


    free(m[nrl]+ncl);
    free(m+nrl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static float **MATRIX_FLT(long nrl, long nrh, long ncl, long nch, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a float matrix m[i][j] with subscript range
            nrl<=i<=nrh, ncl<=j<=ncl.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    long byte = (nrow*ncol)*sizeof(float);
    float **m;

    m = (float **)malloc(((nrow)*sizeof(float*)));
    
    if(!m) nrerror("allocation failure 1 in MATRIX_FLT()");
    
    m -= nrl;
    
    m[nrl] = (float *)malloc(byte);
    
    if(!m[nrl]) nrerror("allocation failure 2 in MATRIX_FLT()");

    if(Init) memset(m[nrl], 0, byte);
    
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

static int FREE_MATRIX_FLT(float **m, long nrl, long ncl){
    
    /*######################################################################
      Purpose:
        free a float matrix allocated with MATRIX_FLT().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
        nrl, ncl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free(m[nrl]+ncl);
    free(m+nrl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static double **MATRIX_DBL(long nrl, long nrh, long ncl, long nch, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a double matrix m[i][j] with subscript range
            nrl<=i<=nrh, ncl<=j<=nch..
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    long byte = (nrow*ncol)*sizeof(double);
    double **m;
    m = (double **)malloc(((nrow)*sizeof(double*)));
    
    if(!m) nrerror("allocation failure 1 in MATRIX_DBL()");
    
    m -= nrl;
    
    m[nrl] = (double *)malloc(byte);

    if(!m[nrl]) nrerror("allocation failure 2 in MATRIX_DBL()");
    
    if(Init) memset(m[nrl], 0, byte);
    
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

static int FREE_MATRIX_DBL(double **m, long nrl, long ncl){
    
    /*######################################################################
      Purpose:
        free a double matrix allocated with MATRIX_DBL().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
        nrl, ncl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((m[nrl]+ncl));
    free((m+nrl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static complex double **MATRIX_CPLX(long nrl, long nrh, long ncl, long nch, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a complex matrix m[i][j] with subscript range
             nrl<=i<=nrh, ncl<=j<=nch.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    long byte = (nrow*ncol)*sizeof(complex double);
    complex double **m;

    m = (complex double **)malloc(((nrow)*sizeof(complex double*)));

    if(!m) nrerror("allocation failure 1 in MATRIX_CPLX()");
    
    m -= nrl;
    
    m[nrl] = (complex double *)malloc(byte);

    if(!m[nrl]) nrerror("allocation failure 2 in MATRIX_CPLX()");

    if(Init) memset(m[nrl], 0, byte);
    
    m[nrl] -= ncl;

    for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

static int FREE_MATRIX_CPLX(complex double **m, long nrl, long ncl){
    
    /*######################################################################
      Purpose:
        free a complex matrix allocated with MATRIX_COMPLEX().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
        nrl, ncl, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free((m[nrl]+ncl));
    free((m+nrl));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern void *VECTOR(long nl, long nh, enum data_type type, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a vector v[i] with subscript range nl<=i<=nh.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        nl, nh, the subscript range.
        type, data type of the matrix (enum_int, enum_flt, enum_dbl,
            enum_cplx, enum_char).
        Init, intialize the vector with 0 or not.
      Return:
        return the pointer.
     ######################################################################*/

    void *v;
    int *p1;
    float *p2;
    double *p3;
    complex double *p4;
    
    switch (type){
      case enum_int:
        p1 = VECTOR_INT(nl, nh, Init);
        v = p1;
        break;
      case enum_flt:
        p2 = VECTOR_FLT(nl, nh, Init);
        v = p2;
        break;
      case enum_dbl:
        p3 = VECTOR_DBL(nl, nh, Init);
        v = p3;
        break;
      case enum_cplx:
        p4 = VECTOR_CPLX(nl, nh, Init);
        v = p4;
        break;

      default:
        nrerror("data type error in VECTOR()");
        v = NULL;
        break;
    }

    return v;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_VECTOR(void *v, long nl, enum data_type type){
    
    /*######################################################################
      Purpose:
        free an int vector allocated with VECTOR().
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        *v, the pointer.
        nl, the first index of subscript range.
        type, data type of the matrix (enum_int, enum_flt, enum_dbl,
            enum_cplx, enum_char).
      Return:
        .
     ######################################################################*/
    
    switch (type){
      case enum_int:
        FREE_VECTOR_INT((int *)v, nl);
        break;
      case enum_flt:
        FREE_VECTOR_FLT((float *)v, nl);
        break;
      case enum_dbl:
        FREE_VECTOR_DBL((double *)v, nl);
        break;
      case enum_cplx:
        FREE_VECTOR_CPLX((double complex *)v, nl);
        break;

      default:
        nrerror("data type error in FREE_VECTOR()");
        break;
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern void *MATRIX(long nrl, long nrh, long ncl, long nch, enum data_type type, \
    bool Init){
    
    /*######################################################################
      Purpose:
        allocate a matrix m[i][j] with subscript range nrl<=i<=nrh,
            ncl<=j<=nch.
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        nrl, nrh, ncl, nch, the subscript range.
        type, data type of the matrix (enum_int, enum_flt, enum_dbl,
            enum_cplx, enum_char).
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
     ######################################################################*/

    void *m;
    int **p1;
    float **p2;
    double **p3;
    complex double **p4;
    char **p5;

    switch (type){
      case enum_int:
        //size = sizeof(int);
        //psize = sizeof(int *);
        p1 = MATRIX_INT(nrl, nrh, ncl, nch, Init);
        m = p1;
        break;
      case enum_flt:
        //size = sizeof(float);
        //psize = sizeof(float *);
        p2 = MATRIX_FLT(nrl, nrh, ncl, nch, Init);
        m = p2;
        break;
      case enum_dbl:
        //size = sizeof(double);
        //psize = sizeof(double *);
        p3 = MATRIX_DBL(nrl, nrh, ncl, nch, Init);
        m = p3;
        break;
      case enum_cplx:
        //size = sizeof(complex double);
        //psize = sizeof(complex double *);
        p4 = MATRIX_CPLX(nrl, nrh, ncl, nch, Init);
        m = p4;
        break;
      case enum_char:
        //size = sizeof(char);
        //psize = sizeof(char *);
        p5 = MATRIX_CHAR(nrl, nrh, ncl, nch, Init);
        m = p5;
        break;
      default:    
        nrerror("data type error in MATRIX()");
        m = NULL;
        break;
    }
 
    return m;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_MATRIX(void *m, long nrl, long ncl, enum data_type type){
    
    /*######################################################################
      Purpose:
        free a matrix allocated with MATRIX().
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        m, the pointer.
        nrl, ncl, the first index of subscript range.
        type, data type of the matrix (enum_int, enum_flt, enum_dbl,
            enum_cplx, enum_char).
      Return:
        .
     ######################################################################*/

    switch(type){
      case enum_int:
        FREE_MATRIX_INT((int **)m, nrl, ncl);
        break;
      case enum_flt:
        FREE_MATRIX_FLT((float **)m, nrl, ncl);
        break;
      case enum_dbl:
        FREE_MATRIX_DBL((double **)m, nrl, ncl);
        break;
      case enum_cplx:
        FREE_MATRIX_CPLX((double complex **)m, nrl, ncl);
        break;
      case enum_char:
        FREE_MATRIX_CHAR((char **)m, nrl, ncl);
        break;

      default:
        nrerror("data type error in FREE_MATRIX()");
        break;
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern float **MATRIX_TRI_FLT(long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a float matrix m[i][j] with subscript range
            0<=i<=n1, 0<=j<=i.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nrow = nh+1,i;
    long byte = (nrow+1)*nrow/2*sizeof(float);
    float **m;

    m = (float **)malloc(nrow*sizeof(float*));
    
    if(!m) nrerror("allocation failure 1 in MATRIX_TRI_FLT()");
    
    m[0] = (float *)malloc(byte);

    if(!m[0]) nrerror("allocation failure 2 in MATRIX_TRI_FLT()");
    
    if(Init) memset(m[0], 0, byte);
   
    for(i=1; i<=nh; i++) m[i] = m[i-1]+i;
    
    return m;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_MATRIX_TRI_FLT(float **m){
    
    /*######################################################################
      Purpose:
        free an float matrix allocated with MATRIX_TRI_FLT().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **m, the pointer.
      Return:
        .
     ######################################################################*/

    free(m[0]);
    free(m);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern complex double **MATRIX_TRI_CPLX(long n1, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a complex matrix m[i][j], with subscript range
            0<=i<=n1, 0<=j<=i.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        n1, the subscript range..
        Init, initialization flag, if true, initialized the matrix with 0.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nrow = n1+1,i;
    long byte = (int)((nrow+1)*nrow/2)*sizeof(complex double);
    complex double **R;

    R = (complex double **)malloc(nrow*sizeof(complex double *));
    
    if(!R) nrerror("allocation failure 1 in MATRIX_TRI_CPLX()");
    
    R[0] = (complex double *)malloc(byte);
    
    if(!R[0]) nrerror("allocation failure 2 in MATRIX_TRI_CPLX()");

    if(Init) memset(R[0], 0, byte);
    
    for(i=1; i<=n1; i++) R[i] = R[i-1]+i;
    
    return R;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_MATRIX_TRI_CPLX(complex double **Rho){
    
    /*######################################################################
      Purpose:
        free a complex matrix allocated by MATRIX_TRI_CPLX().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **Rho, the pointer.
      Return:
        .
     ######################################################################*/

    free((Rho[0]));
    free((Rho));
    
    return 0; 
}

/*--------------------------------------------------------------------------------*/

extern complex double **MATRIX_RHO_CPLX(long nh, bool Init){
    
    /*######################################################################
      Purpose:
        allocate an complex matrix R[i][j] with subscript range
            0<=i<=nh, -i<=j<=i.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        nh, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nrow = nh+1,i;
    complex double **R;
    long byte = nrow*nrow*sizeof(complex double);

    R = (complex double **)malloc(nrow*sizeof(complex double *));
    
    if(!R) nrerror("allocation failure 1 in MATRIX_RHO_CPLX()");
    
    R[0] = (complex double *)malloc(byte);
    
    if(!R[0]) nrerror("allocation failure 2 in MATRIX_RHO_CPLX()");

    if(Init) memset(R[0], 0, byte);
    
    for(i=1; i<=nh; i++) R[i] = R[i-1]+2*i;
    
    return R;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_MATRIX_RHO_CPLX(complex double **Rho){
    
    /*######################################################################
      Purpose:
        free a complex tensor allocated by MATRIX_RHO().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **Rho, the pointer.
      Return:
        .
     ######################################################################*/

    free((Rho[0]));
    free((Rho));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern double ***TENSOR_DBL(long Ni0, long Ni1, long Nj0, long Nj1, long Nk0, \
                            long Nk1, bool Init){
    
    /*######################################################################
      Purpose:
        allocate an double tensor C[i][j][k] with subscript range
            Ni0<=i<=Ni1, Nj0<=j<=Nj1, Nk0<=k<=Nk1.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        Ni0, Ni1, Nj0, Nj1, Nk0, Nk1, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long Ni = Ni1-Ni0+1, Nj = Nj1-Nj0+1, Nk = Nk1-Nk0+1, i, j;
    double ***C;
    
    long byte1 = Ni*sizeof(double **);
    long byte2 = (Ni*Nj)*sizeof(double *);
    long byte3 = (Ni*Nj*Nk)*sizeof(double);

    C = (double ***)malloc(byte1);

    if(!C) nrerror("allocation failure 1 in TENSOR_DBL()");

    C -= Ni0;
    
    C[Ni0] = (double **)malloc(byte2);

    if(!C[Ni0]) nrerror("allocation failure 2 in TENSOR_DBL()");

    C[Ni0] -= Nj0;

    for(i=Ni0+1; i<=Ni1; i++) C[i] = C[i-1]+Nj;
    
    C[Ni0][Nj0] = (double *)malloc(byte3);

    if(!C[Ni0][Nj0]) nrerror("allocation failure 3 in TENSOR_DBL()");
    
    if(Init) memset(C[Ni0][Nj0], 0, byte3);
  
    C[Ni0][Nj0] -= Nk0;

    for(i=Ni0+1; i<=Ni1; i++) C[i][Nj0] = C[i-1][Nj0]+Nj*Nk;
    
    for(i=Ni0; i<=Ni1; i++){
      for(j=Nj0+1; j<=Nj1; j++){
        C[i][j] = C[i][j-1]+Nk;
      }
    }
    
    return C;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_TENSOR_DBL(double ***T, long Ni0, long Nj0, long Nk0){
    
    /*######################################################################
      Purpose:
        free a tensor allocated with TENSOR_DBL().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        **T, the pointer.
        Ni0, Nj0, Nk0, the first index of subscript range.
      Return:
        .
     ######################################################################*/

    free(T[Ni0][Nj0]+Nk0);
    free(T[Ni0]+Nj0);
    free((T+Ni0));
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern complex double ***TENSOR_TRI_CPLX(long n1, long n2, bool Init){
    
    /*######################################################################
      Purpose:
        allocate a complex tensor with subscript range m[i][j][k]
            0<=i<=n1, 0<=j<=nh, 0<=k<=j.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        n1, n2, the subscript range.
        Init, intialize the matrix with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nx = n1+1, ny = n2+1, i, j;
    long byte1 = nx*ny*sizeof(complex double *);
    long byte2 = (int)((nx*(ny+1)*ny)/2)*sizeof(complex double);
    complex double ***m;
    
    m = (complex double ***)malloc(nx*sizeof(complex double **));
    
    if(!m) nrerror("allocation failure 1 in TENSOR_TRI_CPLX()");
    
    m[0] = (complex double **)malloc(byte1);
    
    if(!m[0]) nrerror("allocation failure 2 in TENSOR_TRI_CPLX()");
    
    for(i=1; i<nx; i++) m[i] = m[i-1]+ny;
        
    m[0][0] = (complex double *)malloc(byte2);
    
    if(!m[0][0]) nrerror("allocation failure 3 in TENSOR_TRI_CPLX()");
    
    if(Init) memset(m[0][0], 0, byte2);
    
    for(i=1; i<nx; i++)  m[i][0] = m[i-1][0]+(int)(((ny+1)*ny)/2);
    
    for(i=0; i<nx; i++){
      for(j=1; j <ny; j++){
        m[i][j] = m[i][j-1]+j;
      }
    }
    
    return m;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_TENSOR_TRI_CPLX(complex double ***Rho){
    
    /*######################################################################
      Purpose:
        free a complex tensor allocated by TENSOR_TRI_CPLX().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        ***Rho, the pointer.
      Return:
        .
     ######################################################################*/

    free(Rho[0][0]);
    free(Rho[0]);
    free(Rho);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern double ***TENSOR_RHO_DBL(long n1, long n2, bool Init){
    
    /*######################################################################
      Purpose:
        allocate an double tensor m[i][j][k] with subscript
            range 0<=i<=n1, 0<=j<=n2, -j<=k<=j.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        n1, n2, the subscript range.
        Init, intialize the tensor with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nx = n1+1, ny = n2+1, i, j;
    long byte1 = nx*ny*sizeof(double *);
    long byte2 = nx*ny*ny*sizeof(double);
    double ***R;

    R = (double ***)malloc(nx*sizeof(double **));
    
    if(!R) nrerror("allocation failure 1 in TENSOR_RHO_DBL()");
    
    R[0] = (double **)malloc(byte1);

    if(!R[0]) nrerror("allocation failure 2 in TENSOR_RHO_DBL()");

    for(i=1; i<nx; i++) R[i] = R[i-1]+ny;

    R[0][0] = (double *)malloc(byte2);
    
    if(!R[0][0]) nrerror("allocation failure 3 in TENSOR_RHO_DBL()");
    
    if(Init) memset(R[0][0], 0, byte2);

    for(i=1; i<nx; i++) R[i][0] = R[i-1][0]+ny*ny;
    
    for(i=0; i<nx; i++){
      for(j=1; j <ny; j++){
        R[i][j] = R[i][j-1]+2*j;
      }
    }
    
    return R;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_TENSOR_RHO_DBL(double ***Rho){
    
    /*######################################################################
      Purpose:
        free a double tensor allocated by TENSOR_RHO_DBL().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        ***Rho, the pointer.
      Return:
        .
     ######################################################################*/

    free(Rho[0][0]);
    free(Rho[0]);
    free(Rho);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern complex double ***TENSOR_RHO_CPLX(long n1, long n2, bool Init){
    
    /*######################################################################
      Purpose:
        allocate an complex tensor m[i][j][k] with subscript
            range 0<=i<=n1, 0<=j<=n2, -j<=k<=j.
      Record of revisions:
        30 Oct. 2022.
      Input parameters:
        n1, n2, the subscript range
        Init, intialize the tensor with 0 or not.
      Return:
        return the pointer.
      Reference:
        numerical recipes in C 2ed.
     ######################################################################*/

    long nx = n1+1, ny = n2+1, i, j;
    long byte1 = nx*ny*sizeof(complex double *);
    long byte2 = nx*ny*ny*sizeof(complex double);
    complex double ***R;

    R = (complex double ***)malloc(nx*sizeof(complex double **));
    
    if(!R) nrerror("allocation failure 1 in TENSOR_RHO_CPLX()");

    R[0] = (complex double **)malloc(byte1);
    
    if(!R[0]) nrerror("allocation failure 2 in TENSOR_RHO_CPLX()");

    for(i=1; i<nx; i++) R[i] = R[i-1]+ny;

    R[0][0] = (complex double *)malloc(byte2);
    
    if(!R[0][0]) nrerror("allocation failure 3 in TENSOR_RHO_CPLX()");

    if(Init) memset(R[0][0], 0, byte2);

    for(i=1; i<nx; i++) R[i][0] = R[i-1][0]+ny*ny;
    
    for(i=0; i<nx; i++){
      for(j=1; j <ny; j++){
        R[i][j] = R[i][j-1]+2*j;
      }
    }
    
    return R;
}

/*--------------------------------------------------------------------------------*/

extern int FREE_TENSOR_RHO_CPLX(complex double ***Rho){
    
    /*######################################################################
      Purpose:
        free a complex tensor allocated by TENSOR_RHO_CPLX().
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        ***Rho, the pointer.
      Return:
        .
     ######################################################################*/

    free(Rho[0][0]);
    free(Rho[0]);
    free(Rho);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/
