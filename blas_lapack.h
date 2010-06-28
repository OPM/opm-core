// C <- a1*op(A)*op(B) + a2*C  where  op(X) in {X, X.'}
void dgemm_(char *transA        , char *transB         ,
            const MAT_SIZE_T*  m, const MAT_SIZE_T* n  , const MAT_SIZE_T* k  ,
            const double*     a1, const double*     A  , const MAT_SIZE_T* ldA,
            const double*      B, const MAT_SIZE_T* ldB,
            const double*     a2,       double*     C  , const MAT_SIZE_T* ldC);

// C <- a1*A*A' + a2*C   *or*   C <- a1*A'*A + a2*C
void dsyrk_(char             *uplo, char             *trans,
            const MAT_SIZE_T *n   , const MAT_SIZE_T *k    ,
            const double     *a1  , const double     *A    , const MAT_SIZE_T *ldA,
            const double     *a2  ,       double     *C    , const MAT_SIZE_T *ldC);


// B <- a*op(A)*B  *or*  B <- a*B*op(A)  where op(X) \in {X, X.', X'}
void dtrmm_(char *, char *, char *, char *,
            const MAT_SIZE_T *m, const MAT_SIZE_T* n  ,
            const double     *a,
            const double     *A, const MAT_SIZE_T* ldA,
                  double     *B, const MAT_SIZE_T* ldB);


void dgeqrf_(const MAT_SIZE_T *m    , const MAT_SIZE_T *n   ,
                   double     *A    , const MAT_SIZE_T *ld  ,
                   double     *tau  ,       double     *work,
             const MAT_SIZE_T *lwork,       MAT_SIZE_T *info);


void dorgqr_(const MAT_SIZE_T *m   , const MAT_SIZE_T *n    , const MAT_SIZE_T *k  ,
                   double     *A   , const MAT_SIZE_T *ld   , const double     *tau,
                   double     *work, const MAT_SIZE_T *lwork,       MAT_SIZE_T *info);
