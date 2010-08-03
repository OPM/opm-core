#ifndef BLAS_LAPACK_H_INCLUDED
#define BLAS_LAPACK_H_INCLUDED

/* C <- a1*op(A)*op(B) + a2*C  where  op(X) in {X, X.'} */
void dgemm_(const char *transA  , const char *transB   ,
            const MAT_SIZE_T*  m, const MAT_SIZE_T* n  , const MAT_SIZE_T* k  ,
            const double*     a1, const double*     A  , const MAT_SIZE_T* ldA,
            const double*      B, const MAT_SIZE_T* ldB,
            const double*     a2,       double*     C  , const MAT_SIZE_T* ldC);


/* C <- a1*A*A' + a2*C   *or*   C <- a1*A'*A + a2*C */
void dsyrk_(const char       *uplo, const char       *trans,
            const MAT_SIZE_T *n   , const MAT_SIZE_T *k    ,
            const double     *a1  , const double     *A    , const MAT_SIZE_T *ldA,
            const double     *a2  ,       double     *C    , const MAT_SIZE_T *ldC);


void dgeqrf_(const MAT_SIZE_T *m    , const MAT_SIZE_T *n   ,
                   double     *A    , const MAT_SIZE_T *ld  ,
                   double     *tau  ,       double     *work,
             const MAT_SIZE_T *lwork,       MAT_SIZE_T *info);


void dorgqr_(const MAT_SIZE_T *m   , const MAT_SIZE_T *n    , const MAT_SIZE_T *k  ,
                   double     *A   , const MAT_SIZE_T *ld   , const double     *tau,
                   double     *work, const MAT_SIZE_T *lwork,       MAT_SIZE_T *info);


/* y <- a1*op(A)*x + a2*y */
void dgemv_(const char       *trans,
            const MAT_SIZE_T *m    , const MAT_SIZE_T *n,
            const double     *a1   , const double     *A, const MAT_SIZE_T *ldA ,
                                     const double     *x, const MAT_SIZE_T *incX,
            const double     *a2   ,       double     *y, const MAT_SIZE_T *incY);


/* y <- a*x + y */
void daxpy_(const MAT_SIZE_T *n, const double *a,
            const double *x, const MAT_SIZE_T *incx,
                  double *y, const MAT_SIZE_T *incy);

/* s <- x' * y */
double ddot_(const MAT_SIZE_T *n, const double *x, const MAT_SIZE_T *incx,
                                  const double *y, const MAT_SIZE_T *incy);


#endif /* BLAS_LAPACK_H_INCLUDED */
