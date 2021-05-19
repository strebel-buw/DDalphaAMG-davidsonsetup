
#ifndef DAVIDSON_SETUP_PRECISION_HEADER
  #define DAVIDSON_SETUP_PRECISION_HEADER
  
  #define nrm2_float scnrm2_
  #define nrm2_double dznrm2_
  #define rscl_float csrscl_
  #define rscl_double zdrscl_
  #define lacpy_float clacpy_
  #define lacpy_double zlacpy_
  #define lapmt_float clapmt_
  #define lapmt_double zlapmt_
  #define dlapmt_float slapmt_
  #define dlapmt_double dlapmt_
  #define gemm_float cgemm_
  #define gemm_double zgemm_
  #define geqrf_float cgeqrf_
  #define geqrf_double zgeqrf_
  #define ungqr_float cungqr_
  #define ungqr_double zungqr_
  #define geev_float cgeev_
  #define geev_double zgeev_
  #define ggev_float cggev_
  #define ggev_double zggev_
  #define getrf_float cgetrf_
  #define getrf_double zgetrf_
  #define getri_float cgetri_
  #define getri_double zgetri_

  
  // LAPACK routines
  extern PRECISION nrm2_PRECISION( int *n, complex_PRECISION *x, int *incx );                                                  // vector 2-norm
  extern void rscl_PRECISION( int *n, PRECISION *alpha, complex_PRECISION *x, int *incx );                                     // vector real scale
  extern void lacpy_PRECISION( char *uplo, int *m, int *n, complex_PRECISION *A, int *lda, complex_PRECISION *B, int *ldb );   // matrix copy
  extern void lapmt_PRECISION( int *fwd, int *m, int *n, complex_PRECISION *X, int *ldx, int *k );                             // reorder columns
  extern void dlapmt_PRECISION( int *fwd, int *m, int *n, PRECISION *X, int *ldx, int *k );
  extern void gemm_PRECISION( char* transA, char* transB, int *m, int *n, int *k, complex_PRECISION *alpha,                    // C = a*A*B+b*C 
                              complex_PRECISION *A, int *lda, complex_PRECISION *B, int *ldb, complex_PRECISION *beta, complex_PRECISION *C, int *ldc );
  extern void geqrf_PRECISION( int *m, int *n, complex_PRECISION *A, int *lda, complex_PRECISION *tau,  // computes HH reflectors and R for a QR decomposition
                               complex_PRECISION *work, int *lwork, int *info );
  extern void ungqr_PRECISION( int *m, int *n, int *k, complex_PRECISION *A, int *lda, complex_PRECISION *tau,                 // recover Q from HH reflectors
                               complex_PRECISION *work, int *lwork, int *info );
  extern void geev_PRECISION( char *jobvl, char *jobvr, int *n, complex_PRECISION *A, int *lda, complex_PRECISION *w, complex_PRECISION *vl, 
                              int *ldvl, complex_PRECISION *vr, int *ldvr, complex_PRECISION *work, int *lwork, PRECISION *rwork, int *info );
  extern void ggev_PRECISION( char *jobvl, char *jobvr, int *n, complex_PRECISION *A, int *lda, complex_PRECISION *B, int *ldb, 
                              complex_PRECISION *alpha, complex_PRECISION *beta, complex_PRECISION *vl, int *ldvl, 
                              complex_PRECISION *vr, int *ldvr, complex_PRECISION *work, int *lwork, PRECISION *rwork, int *info ); // solve generalized eigenvalue problem
  extern void getrf_PRECISION( int *m, int *n, complex_PRECISION *A, int *lda, int *piv, int* info );                                                                        // A= L*D*L^H
  extern void getri_PRECISION( int *n, complex_PRECISION *A, int *lda, int *piv, complex_PRECISION* work, int* lwork, int* info ); // A^-1 = (L*D*L^H)^-1
  
  void davidson_setup_PRECISION_struct_init( davidson_setup_PRECISION_struct *d );
  void davidson_setup_PRECISION_struct_alloc( int setup_iter, int maxsize, PRECISION tol, davidson_setup_PRECISION_struct *d, level_struct *l );
  void davidson_setup_PRECISION_struct_free( davidson_setup_PRECISION_struct *d, level_struct *l );

  void initial_davidson_setup_PRECISION( int setup_iter, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading );
  void iterative_davidson_setup_PRECISION( int iterative_iter, level_struct *l, Thread *threading );
#endif
