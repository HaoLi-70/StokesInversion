
#include "allocation.h"
#include "logger.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Typed vector, matrix, and tensor allocation descriptors with checked
        index ranges, contiguous storage, and matching release operations.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

static inline size_t DTYPE_SIZE(DATA_TYPE t){

    /*######################################################################
      Purpose:
        Return the storage size of a supported data type.
      Input parameters:
        t, requested data type.
      Return:
        Size of one value, or 0 for an unsupported type.
    ######################################################################*/

    switch(t){
      case enum_int:  
        return sizeof(int);
      case enum_flt:  
        return sizeof(float);
      case enum_dbl:  
        return sizeof(double);
      case enum_cplx: 
        return sizeof(complex double);
      case enum_char: 
        return sizeof(char);
      default:
        LOG_ERROR(ERR_LVL_ERROR,"DTYPE_SIZE","unknown type");
        return 0;
    }
}

/*--------------------------------------------------------------------------------*/

STRUCT_VECTOR VECTOR(long il, long ih, DATA_TYPE type, bool initialize){

    /*######################################################################
      Purpose:
        Allocate a vector with the requested inclusive index range.
      Input parameters:
        il, ih, lower and upper vector indices.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Vector allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_VECTOR v;
    v.il = il; v.ih = ih;
    v.type = type;

    long n = ih-il+1;
    v.ni = n;
    size_t es = DTYPE_SIZE(type);

    if(initialize){
      v.data = calloc(n, es);
    }else{
      v.data = malloc(n*es);
    }

    if(!v.data){
      LOG_ERROR(ERR_LVL_ERROR,"VECTOR","allocation failure");
    }

    v.data_ptr = (char *)(v.data)-il*es;

    return v;
}

/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX(long il, long ih, long jl, long jh, 
    DATA_TYPE type, bool initialize){

    /*######################################################################
      Purpose:
        Allocate a matrix m[i][j] with subscript range il<=i<=ih, 
            jl<=j<=jh.
      Input parameters:
        il, ih, jl, jh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Matrix allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_MATRIX m;
    m.il = il; m.ih = ih;
    m.jl = jl; m.jh = jh;
    m.type = type;

    long nrows = ih-il+1;
    long ncols = jh-jl+1;
    m.ni = nrows;
    m.nj = ncols;

    size_t es = DTYPE_SIZE(type);

    m.data_ptr = malloc(nrows*sizeof(void *));
    if(!m.data_ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX", "Error: allocation failure 1\n");
    }
    m.data_ptr -= il;

    if(initialize){
      m.data_ptr[il] = calloc(nrows*ncols,es);
    }else{
      m.data_ptr[il] = malloc(nrows*ncols*es);
    }

    if(!m.data_ptr[il]){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX", "Error: allocation failure 2\n");
    }

    m.data = m.data_ptr[il];
    m.data_ptr[il] = (char *)m.data-jl*es;

    for(long i=il+1; i<=ih; i++){
      m.data_ptr[i] = (char *)m.data_ptr[i-1]+ncols*es;
    }

    return m;
}


/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX_TRI(long ih, DATA_TYPE type, bool initialize){

    /*######################################################################
      Purpose:
        Allocate a lower-triangular matrix with range 0<=j<=i<=ih.
      Input parameters:
        ih, maximum row index.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Matrix allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_MATRIX m;
    m.il = 0; m.ih = ih;
    m.type = type;

    long nrow = ih+1;
    long total = nrow*(nrow+1)/2; 

    size_t es = DTYPE_SIZE(type);

    m.data_ptr = malloc(nrow*sizeof(void *));
    if(!m.data_ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_TRI", "Error: allocation failure 1\n");
    }

    if(initialize){
      m.data_ptr[0] = calloc(total,es);
    }else{
      m.data_ptr[0] = malloc(total*es);
    }

    if(!m.data_ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_TRI", "Error: allocation failure 2\n");
    }
  
    m.data = m.data_ptr[0];

    for(long i = 1; i<=ih; i++){
      m.data_ptr[i] = (char *)m.data_ptr[i-1]+i*es;
    }

    return m;
}

/*--------------------------------------------------------------------------------*/

STRUCT_MATRIX MATRIX_RHO(long ih, DATA_TYPE type, bool initialize){
    
    /*######################################################################
      Purpose:
        Allocate a matrix with range 0<=i<=ih and -i<=j<=i.
      Input parameters:
        ih, maximum row index.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Matrix allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_MATRIX m;
    m.il = 0; m.ih = ih;
    m.type = type;

    long nrow = ih+1;
    long total = nrow*nrow;
    size_t es = DTYPE_SIZE(type);

    m.data_ptr = malloc(nrow*sizeof(void *));
    if(!m.data_ptr){
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_RHO", "Error: allocation failure 1\n");
    }

    if(initialize){
      m.data_ptr[0] = calloc(total, es);
    }else{
      m.data_ptr[0] = malloc(total*es);
    }

    if(!m.data_ptr[0]){ 
      LOG_ERROR(ERR_LVL_ERROR, "MATRIX_RHO", "Error: allocation failure 2\n");
    }

    m.data = m.data_ptr[0];

    for(long i=1; i<=ih; i++){
      m.data_ptr[i] = (char *)m.data_ptr[i-1]+2*i*es;
    }

    return m;
}

/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR(long il, long ih, long jl, long jh, long kl, 
    long kh, DATA_TYPE type, bool initialize){
   
    /*######################################################################
      Purpose:
        Allocate a tensor t[i][j][k] with the requested inclusive ranges:
            il<=i<=ih, jl<=j<=jh, kl<=k<=kh.
      Input parameters:
        il, ih, jl, jh, kl, kh, the subscript ranges.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Tensor allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_TENSOR t;
    t.il = il; t.ih = ih;
    t.jl = jl; t.jh = jh;
    t.kl = kl; t.kh = kh;
    t.type = type;

    long in = ih-il+1;
    long jn = jh-jl+1;
    long kn = kh-kl+1;
    
    t.ni = in;
    t.nj = jn;
    t.nk = kn;

    size_t es = DTYPE_SIZE(type);

    t.data_ptr = malloc(in*sizeof(void **));
    if(!t.data_ptr){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 1\n");
    }
    t.data_ptr -= il;

    t.data_ptr[il] = malloc(in*jn*sizeof(void *));
    if(!t.data_ptr[il]){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 2\n");
    }
    t.data_ptr[il] -= jl;
    for(long i=il+1; i<=ih; i++){
      t.data_ptr[i] = t.data_ptr[i-1]+jn;
    }

    if(initialize){
      t.data_ptr[il][jl] = calloc(in*jn*kn, es);
    }else{
      t.data_ptr[il][jl] = malloc(in*jn*kn*es);
    }

    if(!t.data_ptr[il][jl]){ 
      LOG_ERROR(ERR_LVL_ERROR, "TENSOR", "Error: allocation failure 3\n");
    }

    t.data = t.data_ptr[il][jl];
    t.data_ptr[il][jl] = (char*)t.data-kl*es;

    for(long i=il+1; i<=ih; i++){
      t.data_ptr[i][jl] = (char*)t.data_ptr[i-1][jl]+jn*kn*es;
    }

    for(long i=il; i<=ih; i++){
      for(long jj =jl+1; jj<=jh; jj++){
        t.data_ptr[i][jj] = (char*)t.data_ptr[i][jj-1]+kn*es;
      }
    }

    return t;
}

/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR_TRI(long ih, long jh, DATA_TYPE type, bool initialize){

    /*######################################################################
      Purpose:
        Allocate a triangular tensor with ranges
            0<=i<=ih, 0<=j<=jh, 0<=k<=j.
      Input parameters:
        ih, maximum first-dimension index.
        jh, maximum second-dimension index.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Tensor allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_TENSOR t;
    t.il = t.jl = t.kl = 0;                         
    t.ih = ih;
    t.jh = t.kh = jh;
    t.type = type;

    long nx = ih+1;
    long ny = jh+1;
    t.ni = nx;
    t.nj = ny;

    size_t es = DTYPE_SIZE(type);

    t.data_ptr = malloc(nx*sizeof(void **));
    if(!t.data_ptr){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", 
          "Error: allocation failure 1");
    }

    t.data_ptr[0] = malloc(nx*ny*sizeof(void *));
    if(!t.data_ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", 
          "Error: allocation failure 2");
    }

    for(long i=1; i<nx; i++){
      t.data_ptr[i] = t.data_ptr[i-1]+ny;
    }

    long tri = ny*(ny+1)/2;
    long total = nx*tri;

    if(initialize){
      t.data_ptr[0][0] = calloc(total, es);
    }else{
      t.data_ptr[0][0] = malloc(total*es);
    }

    if(!t.data_ptr[0][0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_TRI", 
          "Error: allocation failure 3");
    }
    
    t.data = t.data_ptr[0][0];

    for(long i=1; i<nx; i++){
      t.data_ptr[i][0] = (char *)t.data_ptr[i-1][0]+tri*es;
    }

    for(long i=0; i<nx; i++){
      for(long jj=1; jj<ny; jj++){
        t.data_ptr[i][jj] = (char*)t.data_ptr[i][jj-1]+jj*es;
      }
    }

    return t;
}

/*--------------------------------------------------------------------------------*/

STRUCT_TENSOR TENSOR_RHO(long ih,long jh, DATA_TYPE type,bool initialize){

    /*######################################################################
      Purpose:
        Allocate a tensor with ranges
            0<=i<=ih, 0<=j<=jh, -j<=k<=j.
      Input parameters:
        ih, maximum first-dimension index.
        jh, maximum second-dimension index.
        type, data type (enum_int, enum_flt, enum_dbl, enum_cplx, 
            enum_char).
        initialize, whether to zero-initialize the allocation.
      Return:
        Tensor allocation descriptor.
      Note:
        Numerical Recipes in C 2nd edition.
    ######################################################################*/

    STRUCT_TENSOR t;

    t.il = t.jl = t.kl = 0;                         
    t.ih = ih;
    t.jh = t.kh = jh;
    t.type = type;

    long nx = ih+1;
    long ny = jh+1;

    size_t es = DTYPE_SIZE(type);

    t.data_ptr = malloc(nx*sizeof(void **));
    if(!t.data_ptr){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", 
          "Error: allocation failure 1");
    }

    t.data_ptr[0] = malloc(nx*ny*sizeof(void *));
    if(!t.data_ptr[0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", 
          "Error: allocation failure 2");
    }

    for(long i=1; i<nx; i++) t.data_ptr[i] = t.data_ptr[i-1]+ny;

    long tri = ny*ny;
    long total = nx*tri;

    if(initialize){
      t.data_ptr[0][0] = calloc(total, es);
    }else{
      t.data_ptr[0][0] = malloc(total*es);
    }

    if(!t.data_ptr[0][0]){
      LOG_ERROR(ERR_LVL_ERROR,"TENSOR_RHO", 
          "Error: allocation failure 3");
    }

    t.data = t.data_ptr[0][0];

    for(long i=1; i<nx; i++){
      t.data_ptr[i][0] = (char*)t.data_ptr[i-1][0]+tri*es;
    }

    for(long i=0; i<nx; i++){
      for(long jj=1; jj<ny; jj++){
        t.data_ptr[i][jj] = (char*)t.data_ptr[i][jj-1]+(2*jj)*es;
      }
    }

    return t;
}

/*--------------------------------------------------------------------------------*/
