
#include "main.h"

void davidson_setup_PRECISION_struct_init( davidson_setup_PRECISION_struct *d ) {
  
  d->lambda = NULL;
  d->theta = NULL;
  d->H1 = NULL;
  d->H2 = NULL;
  d->s = NULL;
  d->Au = NULL;
  d->u = NULL;
  d->res = NULL;
  d->X = NULL;
  d->V = NULL;
  d->AV = NULL;
  d->W = NULL;
  d->iter = 0;
  d->size = 0;
  d->nconv = 0;
  d->time = 0.0;
  d->initial_phase = 1;
  d->arnoldi = 0;
}


void davidson_setup_PRECISION_struct_alloc( int setup_iter, int maxsize, PRECISION tol, davidson_setup_PRECISION_struct *d, level_struct *l ) {
  
  int vs = (!l->depth)?l->inner_vector_size:l->vector_size;
  
  d->setup_iter = setup_iter;
  d->maxsize = maxsize;
  d->tol = tol;
  
  MALLOC( d->lambda, complex_PRECISION, maxsize );
  MALLOC( d->theta, complex_PRECISION, maxsize );
  MALLOC( d->H1, complex_PRECISION*, maxsize ); d->H1[0] = NULL;
  MALLOC( d->H1[0], complex_PRECISION, SQUARE(maxsize));
  MALLOC( d->H2, complex_PRECISION*, maxsize ); d->H2[0] = NULL;
  MALLOC( d->H2[0], complex_PRECISION, SQUARE(maxsize));
  MALLOC( d->s, complex_PRECISION, SQUARE(maxsize) );
  MALLOC( d->Au, complex_PRECISION, vs );
  MALLOC( d->u, complex_PRECISION, vs );
  MALLOC( d->res, complex_PRECISION, vs );
  MALLOC( d->X, complex_PRECISION*, maxsize ); d->X[0] = NULL;
  MALLOC( d->X[0], complex_PRECISION, maxsize*vs );
  MALLOC( d->V, complex_PRECISION*, maxsize ); d->V[0] = NULL;
  MALLOC( d->V[0], complex_PRECISION, maxsize*vs );
  MALLOC( d->AV, complex_PRECISION*, maxsize ); d->AV[0] = NULL;
  MALLOC( d->AV[0], complex_PRECISION, maxsize*vs );
  MALLOC( d->W, complex_PRECISION*, maxsize ); d->W[0] = NULL;
  MALLOC( d->W[0], complex_PRECISION, maxsize*vs );
  
  for( int i = 1; i<maxsize; i++ ) {
    d->H1[i] = d->H1[0]+i*maxsize;
    d->H2[i] = d->H2[0]+i*maxsize;
  }
  
  for( int i = 1; i<maxsize; i++ ) {
    d->X[i] = d->X[0]+ i*vs;
    d->V[i] = d->V[0]+ i*vs;
    d->AV[i] = d->AV[0]+ i*vs;
    d->W[i] = d->W[0]+ i*vs;
  }
}


void davidson_setup_PRECISION_struct_free( davidson_setup_PRECISION_struct *d, level_struct *l ) {
  
  int vs = (!l->depth)?l->inner_vector_size:l->vector_size;
  
  FREE( d->lambda, complex_PRECISION, d->maxsize );
  FREE( d->theta, complex_PRECISION, d->maxsize );
  FREE( d->H1[0], complex_PRECISION, SQUARE(d->maxsize));
  FREE( d->H1, complex_PRECISION*, d->maxsize );
  FREE( d->H2[0], complex_PRECISION, SQUARE(d->maxsize));
  FREE( d->H2, complex_PRECISION*, d->maxsize );
  FREE( d->s, complex_PRECISION, SQUARE(d->maxsize) );
  FREE( d->Au, complex_PRECISION, vs );
  FREE( d->u, complex_PRECISION, vs );
  FREE( d->res, complex_PRECISION, vs );
  FREE( d->X[0], complex_PRECISION, d->maxsize*vs );
  FREE( d->X, complex_PRECISION*, d->maxsize );
  FREE( d->V[0], complex_PRECISION, d->maxsize*vs );
  FREE( d->V, complex_PRECISION*, d->maxsize );
  FREE( d->AV[0], complex_PRECISION, d->maxsize*vs );
  FREE( d->AV, complex_PRECISION*, d->maxsize );
  FREE( d->W[0], complex_PRECISION, d->maxsize*vs );
  FREE( d->W, complex_PRECISION*, d->maxsize );
  
}


void test_eigenspace_PRECISION( davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  int vs = l->inner_vector_size, nev = d->size;
  PRECISION norm = 0;
  complex_PRECISION H[nev][nev];
  
  printf0("Checking davidson setup for Level: %d\n", l->level);
  
  // check ||V||
  for( int i=0; i<nev; i++ ) {
    norm = global_norm_PRECISION( d->V[i], 0, vs, l, threading );
    printf0("||V[%2d]|| = %le\n", i, norm );
  }
  
  // test V* V = I
  norm = 0;
  for( int i=0; i<nev; i++ ) {
    for( int j=0; j<nev; j++ ) {
      H[i][j] = global_inner_product_PRECISION( d->V[i], d->V[j], 0, vs, l, threading );
      norm += cabs_PRECISION( H[i][j] );
    }
  }
  norm -= nev;
  printf0("||V* V - I|| = %le ", norm);  
}


void qr_PRECISION( complex_PRECISION *Q, complex_PRECISION *A, int n, int lda, int ldq ) {
  
  int lwork = n*n, info;
  char uplo = 'A';
  complex_PRECISION work[lwork], tau[n];
  
  lacpy_PRECISION( &uplo, &n, &n, A, &lda, Q, &ldq );               // copy A into Q to avoid overwriting of input
  geqrf_PRECISION( &n, &n, Q, &ldq, tau, work, &lwork, &info );     // compute QR factorization
  ungqr_PRECISION( &n, &n, &n, Q, &ldq, tau, work, &lwork, &info ); // recover Q from Householder reflections
}


void matrix_full_inverse_PRECISION( complex_PRECISION *Ainv, complex_PRECISION *A, int n, int lda, int ldainv ) {
  
  int lwork = n*n, piv[n], info;
  char uplo = 'A';
  complex_PRECISION work[lwork];

  lacpy_PRECISION( &uplo, &n, &n, A, &lda, Ainv, &ldainv );    // copy A into Ainv to avoid overwriting of input
  getrf_PRECISION( &n, &n, Ainv, &lda, piv, &info );           // compute A= L*D*L^H
  getri_PRECISION( &n, Ainv, &lda, piv, work, &lwork, &info ); // A^-1 = (L*D*L^H)^-1

}


void dense_ev_PRECISION_driver( void *eig_val, complex_PRECISION *eig_vec, complex_PRECISION **H1, complex_PRECISION **H2, int n ) {
  
  int lwork = 2*n*n, pi[n], info;
  PRECISION rwork[8*n];
  char N = 'N', V = 'V';
  complex_PRECISION *A=NULL, *B=NULL, *Binv_A=NULL, *alpha=NULL, *beta=NULL, *work=NULL, one = 1, zero = 0;
  
  ASSERT(lwork >= n*n);
  
  MALLOC( A, complex_PRECISION, n*n );
  MALLOC( B, complex_PRECISION, n*n );
  MALLOC( Binv_A, complex_PRECISION, n*n );
  MALLOC( alpha, complex_PRECISION, n );
  MALLOC( beta, complex_PRECISION, n );
  MALLOC( work, complex_PRECISION, lwork );

  // local copy for ggev
  for( int i = 0; i<n; i++ )
    for( int j = 0; j<n; j++ ) {
      A[i*n+j] = H1[j][i];
      if( H2 != NULL ) B[i*n+j] = H2[j][i];
    }

  if( H2 == NULL ) {
    complex_PRECISION *ev = (complex_PRECISION*) eig_val;
    
     // LAPACK routine for standard eigenvalue problem
    geev_PRECISION( &N, &V, &n, A, &n, ev, NULL, &n, eig_vec, &n, work, &lwork, rwork, &info );
    if( info ) error0("LAPACK error in geev_PRECISION! Errorcode: %d\n", info);
    
    // sort eigenvalues ascending in absolute value
    for( int i=0; i<n; i++ ) pi[i] = i+1;
    for( int i=n; i>1; i-- )
      for (int j=0; j<i-1; j++ )
        if ( cabs(ev[pi[j]-1]) > cabs(ev[pi[j+1]-1]) ) {
          info = pi[j];
          pi[j] = pi[j+1];
          pi[j+1] = info;
        }
    info = 1;
    lapmt_PRECISION( &info, &n, &n, eig_vec, &n, pi );
    lapmt_PRECISION( &info, &info, &n, ev, &info, pi );
  }
  else {
    complex_PRECISION *ev = (complex_PRECISION*) eig_val;

    // transform A x = lambda B x to B^-1 A x = lambda x and solve standard eigenvalue problem
    matrix_full_inverse_PRECISION( work, B, n, n, n); // buffer B^-1 into work
    gemm_PRECISION( &N, &N, &n, &n, &n, &one, work, &n, A, &n, &zero, Binv_A, &n ); //Binv_A = B^-1 A
    geev_PRECISION( &N, &V, &n, Binv_A, &n, ev, NULL, &n, eig_vec, &n, work, &lwork, rwork, &info );
    if( info ) error0("LAPACK error in geev_PRECISION! Errorcode: %d\n", info);
    
//     // LAPACK routine for generalized eigenvalue problem
//     ggev_PRECISION( &N, &V, &n, A, &n, B, &n, alpha, beta, NULL, &n, eig_vec, &n, work, &lwork, rwork, &info );
//     if( info ) error0("LAPACK error in ggev_PRECISION! Errorcode: %d\n", info);
//     // obtain eigenvalues and cast to PRECISION
//     for ( int i=0; i<n; i++ ) ev[i] = (complex_PRECISION) alpha[i]/beta[i]; //NOTE: division might be dangerous (need kappa(A) ~ kappa(B) to be safe)
    
//   printf0("s=\n");
//   for (int i = 0; i<5; i++) {
//     for (int j = 0; j<5; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(*(eig_vec+i*n+j)));
//     }
//     printf0("\n");
//   }
//   getchar();
    
    // sort eigenvalues ascending in absolute value
    for( int i=0; i<n; i++ ) pi[i] = i+1;
    for( int i=n; i>1; i-- )
      for (int j=0; j<i-1; j++ )
        if ( cabs(ev[pi[j]-1]) > cabs(ev[pi[j+1]-1]) ) {
          info = pi[j];
          pi[j] = pi[j+1];
          pi[j+1] = info;
        }
    info = 1;
    lapmt_PRECISION( &info, &n, &n, eig_vec, &n, pi );
    lapmt_PRECISION( &info, &info, &n, ev, &info, pi );
  }
  
FREE( A, complex_PRECISION, n*n );
FREE( B, complex_PRECISION, n*n );
FREE( Binv_A, complex_PRECISION, n*n );
FREE( alpha, complex_PRECISION, n );
FREE( beta, complex_PRECISION, n );
FREE( work, complex_PRECISION, lwork );
}


void define_next_level_PRECISION( davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {

  int vs = l->inner_vector_size, n = l->num_eig_vect; 
  
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
  for ( int k=0; k<n; k++ ) {
    vector_PRECISION_copy( l->is_PRECISION.test_vector[k], d->V[k], 0, vs, l );
//     vector_PRECISION_copy( l->is_PRECISION.test_vector[k], k < d->nconv?d->X[k]:d->V[k], 0, vs, l );
    vector_PRECISION_copy( l->is_PRECISION.interpolation[k], l->is_PRECISION.test_vector[k], 0, vs, l );
  }
#endif
  
  testvector_analysis_PRECISION( l->is_PRECISION.test_vector, l, threading );
if( d->initial_phase ) {
  #ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
    define_interpolation_PRECISION_operator( l->is_PRECISION.test_vector, l, threading );
    gram_schmidt_on_aggregates_PRECISION_vectorized( l->is_PRECISION.operator, n, l, threading );
  #else
    gram_schmidt_on_aggregates_PRECISION( l->is_PRECISION.interpolation, n, l, threading );
    define_interpolation_PRECISION_operator( l->is_PRECISION.interpolation, l, threading );
  #endif
  }
}


void get_residual_PRECISION( vector_PRECISION res, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  /* Computes residual for next preconditioner step */
  
  int k, found, m = d->size, vs = l->inner_vector_size;
  PRECISION norm, theta[m], lambda[m];
  complex_PRECISION *s = d->s, v[m*m];
  vector_PRECISION u = d->u, Au = d->Au;
  
  FILE *fop;
  fop = fopen("H1H2.m","w");
  fprintf(fop,"H1 = [\n");
  for( int i = 0; i<m; i++ ) {
    for( int j = 0; j<m; j++ ) {
      fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H1[i][j]));
    }
    fprintf(fop,";\n");
  }
  fprintf(fop,"];\n");
  fprintf(fop,"H2 = [\n");
  for( int i = 0; i<m; i++ ) {
    for( int j = 0; j<m; j++ ) {
      fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H2[i][j]));
    }
    fprintf(fop,";\n");
  }
  fprintf(fop,"];\n");
//   printf0("done\n");
//   getchar();
  
  // compute eigenvalues
  dense_ev_PRECISION_driver( d->theta, v, d->H1, d->H2, m );
  qr_PRECISION( s, v, m, m, d->maxsize );
  
  fprintf(fop,"VV=[\n");
  for (int i = 0; i<m; i++) {
    for (int j = 0; j<m; j++) {
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
      fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(*(v+i*m+j)));
    }
    fprintf(fop,";\n");
  }
  fprintf(fop,"];\n VV = VV.';\n");
  fprintf(fop,"theta = [\n");
  for (int i = 0; i<m; i++) {
    fprintf(fop,"%+.16lf%+.16lfi\n", CSPLIT(d->theta[i]));
  }
  fprintf(fop,"];\n");
  
  printf0("dense ev end\n");
  fclose(fop);
  getchar();
  
  for( int i=0; i<m; i++) {
    theta[i] = creal(d->theta[i]);
    lambda[i] = creal(d->lambda[i]);
  }
//   printf0("s=\n");
//   for (int i = 0; i<m; i++) {
//     for (int j = 0; j<m; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
//     }
//     printf0("\n");
//   }
//   getchar();
//   printf0("initial phase = %d\n",d->initial_phase);
//   if (!d->initial_phase) {
//     memset(s, 0, sizeof(complex_PRECISION)*d->maxsize*d->maxsize);
//     for (int i = 0; i<m; i++) {
//       s[i*d->maxsize+i] = 1;
//       theta[i] = i+m+1;
//       getchar();
//     }
//   }  
  // locate next non-converged eigenvalue
  for( k=0; k<m; k++ ) {
    found = 0;
    for( int i=0; i<d->nconv; i++ ) {
        if( fabs(lambda[i] - theta[k]) < d->tol ) {
          found = 1;
          break;
        }
    }
    if( !found )
      break;
  }

//   for (int i = 0; i<m; i++) {
//     printf0("theta = %+.16lf\n", theta[i]);
//   }
//   if (!d->initial_phase) getchar();
  
  // project ritz vector up
  vector_PRECISION_define( u, 0, 0, vs, l );
  vector_PRECISION_multi_saxpy( u, d->V, s+k*d->maxsize, 1, m, 0, vs, l );
  norm = global_norm_PRECISION( u, 0, vs, l, threading );
  vector_PRECISION_scale( u, u, 1/norm, 0, vs, l );
  vector_PRECISION_define( Au, 0, 0, vs, l );
  vector_PRECISION_multi_saxpy( Au, d->AV, s+k*d->maxsize, 1, m, 0, vs, l );
  gamma5_PRECISION( Au, Au, l, threading );
  
  printf0("k= %d, theta[k] = %+.16lf\n", k, theta[k]);
  
//   printf0("s[k]=\n");
//   for (int j = 0; j<5; j++) {
//   printf0("%+.16lf%+.16lfi\n", CSPLIT(*(s+k*d->maxsize+j)));
//     }
//   printf0("\n");
//     
//   printf0("d->V=\n");
//   for (int i = 0; i<3; i++) {
//     for (int j = 0; j<5; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(d->V[i][j]));
//     }
//     printf0("\n");
//   }
//   getchar();
//   
//     printf0("d->AV=\n");
//   for (int i = 0; i<2; i++) {
//     for (int j = 0; j<5; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(d->AV[i][j]));
//     }
//     printf0("\n");
//   }
//   getchar();
//   
//   printf0("u=\n");
//   for (int j = 0; j<5; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(u[j]));
//   }
//   printf0("\n");
//   
//   printf0("Au=\n");
//   for (int j = 0; j<5; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(Au[j]));
//   }
//   printf0("\n");
  
  // compute and check residual
  vector_PRECISION_saxpy( res, Au, u, -theta[k], 0, vs, l );
  norm = global_norm_PRECISION( res, 0, vs, l, threading );
//   printf0("\nk = %d, theta = %.14lf, norm(r) = %.14le\n", k, theta[k], norm);
//   printf0("V=[");
//   for (int j = 0; j<3072; j++) {
//     printf0("%+.16lf%+.16lfi; ", CSPLIT(Au[j]));
//   }
//   printf0("];\n");
  
//   getchar();
  
  if( norm < d->tol ) {
    printf0("Lvl %d, iter: %d: Vector %d reached desired accuracy, eigenvalue: %+.16lf\n", l->level, m, d->nconv, theta[k]);
    d->lambda[d->nconv] = creal(theta[k]);
    vector_PRECISION_copy( d->X[d->nconv], u, 0, vs, l );
    d->nconv++; 

//   printf0("d->lambda=\n");
//   for (int j = 0; j<d->nconv; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(d->lambda[j]));
//   }
//   printf0("\n");
//   
//   printf0("theta=\n");
//   for (int j = 0; j<m; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(theta[j]));
//   }
//   printf0("\n");
//     
    // locate next non-converged eigenvalue
    for( ; k<m; k++ ) {
      found = 0;
      for( int i=0; i<d->nconv; i++ ) {
        if( fabs(d->lambda[i]-theta[k]) < d->tol ) {
          found = 1;
          break;
        }
      }
      if( !found )
        break;
    }
  
//   printf0("s_full=\n");
//   for (int i = 0; i<m; i++) {
//     for (int j = 0; j<m; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
//     }
//     printf0("\n");
//   }
//   printf0("\n");
//   
//   printf0("s_new=\n");
//   for (int j = 0; j<5; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(*(s+k*d->maxsize+j)));
//   }
//   printf0("\n");

    // project ritz vector up
    vector_PRECISION_define( u, 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( u, d->V, s+k*d->maxsize, 1, m, 0, vs, l );
    norm = global_norm_PRECISION( u, 0, vs, l, threading );
    vector_PRECISION_scale( u, u, 1/norm, 0, vs, l );
    
    // recompute ritz estimate from ritz vector (instead of using harmonic ritz value)
//     apply_operator_PRECISION( Au, u, &(l->p_PRECISION), l, threading );
    norm = theta[k]; //global_inner_product_PRECISION( u, Au, 0, vs, l, threading );
    
    vector_PRECISION_define( Au, 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( Au, d->AV, s+k*d->maxsize, 1, m, 0, vs, l );  
    gamma5_PRECISION( Au, Au, l, threading );
    vector_PRECISION_saxpy(res, Au, u, -norm, 0, vs, l );
    norm = global_norm_PRECISION( res, 0, vs, l, threading );
//     printf0("\nnew: k = %d, theta = %.14lf, norm(r) = %.14le\n", k, theta[k], norm);
//     getchar();
  }
//   printf0("residual end\n");
}


void get_residual_PRECISION_on_D( vector_PRECISION res, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  /* Computes residual for next preconditioner step */
  
  int k, found, m = d->size, vs = l->inner_vector_size;
  PRECISION norm;
  complex_PRECISION *s = d->s, v[m*m];
  vector_PRECISION u = d->u, Au = d->Au;
  
  // compute eigenvalues
  dense_ev_PRECISION_driver( d->theta, v, d->H1, NULL, m );
  qr_PRECISION( s, v, m, m, d->maxsize );

//   printf0("s=\n");
//   for (int i = 0; i<m; i++) {
//     for (int j = 0; j<m; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
//     }
//     printf0("\n");
//   }
//   getchar();
  
  // locate next non-converged eigenvalue
  for( k=0; k<m; k++ ) {
    found = 0;
    for( int i=0; i<d->nconv; i++ ) {
        if( cabs(d->lambda[i] - d->theta[k]) < d->tol ) {
          found = 1;
          break;
        }
    }
    if( !found )
      break;
  }

//   for (int i = 0; i<m; i++) {
//     printf0("theta = %+.16lf%+.16lfi\n", CSPLIT(d->theta[i]));
//   }
//   if (!d->initial_phase) getchar();
  
  // project ritz vector up
  vector_PRECISION_define( u, 0, 0, vs, l );
  vector_PRECISION_multi_saxpy( u, d->V, s+k*d->maxsize, 1, m, 0, vs, l );
  norm = global_norm_PRECISION( u, 0, vs, l, threading );
  vector_PRECISION_scale( u, u, 1/norm, 0, vs, l );
  vector_PRECISION_define( Au, 0, 0, vs, l );
  vector_PRECISION_multi_saxpy( Au, d->AV, s+k*d->maxsize, 1, m, 0, vs, l );
//   gamma5_PRECISION( Au, Au, l, threading );
  
//   printf0("s[k]=\n");
//   for (int j = 0; j<m; j++) {
//   printf0("%+.16lf%+.16lfi ", CSPLIT(*(s+k*d->maxsize+j)));
//     }
//   printf0("\n");
//     
//   printf0("d->V=\n");
//   for (int i = 0; i<3; i++) {
//     for (int j = 0; j<5; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(d->V[i][j]));
//     }
//     printf0("\n");
//   }
//   getchar();
//   
//     printf0("d->AV=\n");
//   for (int i = 0; i<2; i++) {
//     for (int j = 0; j<5; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(d->AV[i][j]));
//     }
//     printf0("\n");
//   }
//   getchar();
//   
//   printf0("u=\n");
//   for (int j = 0; j<5; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(u[j]));
//   }
//   printf0("\n");
//   
//   printf0("Au=\n");
//   for (int j = 0; j<5; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(Au[j]));
//   }
//   printf0("\n");
  
  // compute and check residual
  vector_PRECISION_saxpy( res, Au, u, -d->theta[k], 0, vs, l );
  norm = global_norm_PRECISION( res, 0, vs, l, threading );
//   printf0("\nk = %d, theta = %+.14lf%+.14lfi, norm(r) = %.14le\n", k, CSPLIT(d->theta[k]), norm);
//   printf0("V=[");
//   for (int j = 0; j<3072; j++) {
//     printf0("%+.16lf%+.16lfi; ", CSPLIT(Au[j]));
//   }
//   printf0("];\n");
  
//   getchar();
  
  if( norm < d->tol ) {
    printf0("Lvl %d, iter: %d: Vector %d reached desired accuracy, eigenvalue: %+.16lf%+.16lfi\n", l->level, m, d->nconv, CSPLIT(d->theta[k]));
    d->lambda[d->nconv] = d->theta[k];
    vector_PRECISION_copy( d->X[d->nconv], u, 0, vs, l );
    d->nconv++; 

//   printf0("d->lambda=\n");
//   for (int j = 0; j<d->nconv; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(d->lambda[j]));
//   }
//   printf0("\n");
//   
//   printf0("theta=\n");
//   for (int j = 0; j<m; j++) {
//     printf0("%+.16lf%+.16lfi\n", CSPLIT(theta[j]));
//   }
//   printf0("\n");
//     
    // locate next non-converged eigenvalue
    for( ; k<m; k++ ) {
      found = 0;
      for( int i=0; i<d->nconv; i++ ) {
        if( cabs(d->lambda[i]-d->theta[k]) < d->tol ) {
          found = 1;
          break;
        }
      }
      if( !found )
        break;
    }
  
//   printf0("s_full=\n");
//   for (int i = 0; i<m; i++) {
//     for (int j = 0; j<m; j++) {
//       printf0("%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
//     }
//     printf0("\n");
//   }
//   printf0("\n");
//   
  printf0("s_new=\n");
  for (int j = 0; j<5; j++) {
    printf0("%+.16lf%+.16lfi\n", CSPLIT(*(s+k*d->maxsize+j)));
  }
  printf0("\n");

    // project ritz vector up
    vector_PRECISION_define( u, 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( u, d->V, s+k*d->maxsize, 1, m, 0, vs, l );
    norm = global_norm_PRECISION( u, 0, vs, l, threading );
    vector_PRECISION_scale( u, u, 1/norm, 0, vs, l );
    
    // recompute ritz estimate from ritz vector (instead of using harmonic ritz value)
//     apply_operator_PRECISION( Au, u, &(l->p_PRECISION), l, threading );
    norm = d->theta[k]; //global_inner_product_PRECISION( u, Au, 0, vs, l, threading );
    
    vector_PRECISION_define( Au, 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( Au, d->AV, s+k*d->maxsize, 1, m, 0, vs, l );  
//     gamma5_PRECISION( Au, Au, l, threading );
    vector_PRECISION_saxpy(res, Au, u, -norm, 0, vs, l );
    norm = global_norm_PRECISION( res, 0, vs, l, threading );
    printf0("\nnew: k = %d, theta = %+.14lf%+.14lfi, norm(r) = %.14le\n", k, CSPLIT(d->theta[k]), norm);
//     getchar();
  }
}


void extend_subspace_PRECISION( vector_PRECISION *V_add, int blocksize, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  /* Updates system matrices after extending the basis */
  
  int m = d->size, vs = l->inner_vector_size;
  PRECISION norm;
  complex_PRECISION y[m];
  vector_PRECISION w = d->W[0], ww = d->W[1], *V = d->V, *AV = d->AV;
//   printf0("extend begin\n");
  if( V_add != NULL ) {
    for( int b = 0; b<blocksize; b++ ) {
      // Add new basis vector to V and AV
      vector_PRECISION_projection( ww, V_add[b], m+b, V, 1, l, threading );
      vector_PRECISION_projection( w, ww, m+b, V, 1, l, threading );
      norm = global_norm_PRECISION( w, 0, vs, l, threading );
      vector_PRECISION_scale( V[m+b], w, 1.0/norm, 0, vs, l );
      apply_operator_PRECISION( AV[m+b], V[m+b], &(l->p_PRECISION), l, threading );

      // Update H1
      process_multi_inner_product_PRECISION( m, y, AV, AV[m+b], 0, vs, l, threading );
      if ( g.num_processes > 1 ) {
        PROF_PRECISION_START( _ALLR );
        MPI_Allreduce( y, d->H1[m+b], m, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
        PROF_PRECISION_STOP( _ALLR, 1 );
        for( int i = 0; i<m; i++ ) {
          d->H1[i][m+b] = d->H1[m+b][i];
          d->H1[m+b][i] = conj_PRECISION(d->H1[m+b][i]);
        }
      } else {
        for( int i = 0; i<m; i++ ) {
          d->H1[i][m+b] = y[i];
          d->H1[m+b][i] = conj_PRECISION(d->H1[i][m+b]);
        }
      }
      d->H1[m+b][m+b] = global_inner_product_PRECISION( AV[m+b], AV[m+b], 0, vs, l, threading );

      // Update H2
      gamma5_PRECISION( w, AV[m+b], l, threading );
      process_multi_inner_product_PRECISION( m, y, V, w, 0, vs, l, threading );
      if ( g.num_processes > 1 ) {
        PROF_PRECISION_START( _ALLR );
        MPI_Allreduce( y, d->H2[m+b], m, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
        PROF_PRECISION_STOP( _ALLR, 1 );
        for( int i = 0; i<m; i++ ) {
          d->H2[i][m+b] = d->H2[m+b][i];
          d->H2[m+b][i] = conj_PRECISION(d->H2[m+b][i]);
        }
      } else {
        for( int i = 0; i<m; i++ ) {
          d->H2[i][m+b] = y[i];
          d->H2[m+b][i] = conj_PRECISION(d->H2[i][m+b]);
        }
      }
      d->H2[m][m] = global_inner_product_PRECISION( V[m+b], w, 0, vs, l, threading );
    }
    d->size += blocksize;
  }
//   printf0("extend end\n");
  
}


void extend_subspace_on_D_PRECISION( vector_PRECISION *V_add, int blocksize, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  /* Updates system matrices after extending the basis */
  
  int m = d->size, vs = l->inner_vector_size;
  PRECISION norm;
  complex_PRECISION y[m];
  vector_PRECISION w = d->W[0], *V = d->V, *AV = d->AV;
  
  if( V_add != NULL ) {
    for( int b = 0; b<blocksize; b++ ) {
      // Add new basis vector to V and AV
      vector_PRECISION_projection( w, V_add[b], m+b, V, 1, l, threading );
      norm = global_norm_PRECISION( w, 0, vs, l, threading );
      vector_PRECISION_scale( V[m+b], w, 1.0/norm, 0, vs, l );
      apply_operator_PRECISION( AV[m+b], V[m+b], &(l->p_PRECISION), l, threading );
      
//       printf0("V=[");
//       for (int j = 0; j<5; j++) {
//         trans_back_PRECISION( AV[(m+1)%20], AV[m], l->s_PRECISION.op.translation_table, l, threading );
//         printf0("%+.16lf%+.16lfi; ", CSPLIT(AV[(m+1)%20][j]));
//       }
//       printf0("];\n");
//       getchar();
        
      // Update H
      process_multi_inner_product_PRECISION( m, y, V, AV[m+b], 0, vs, l, threading );
      if ( g.num_processes > 1 ) {
        PROF_PRECISION_START( _ALLR );
        MPI_Allreduce( y, d->H1[m+b], m, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
        PROF_PRECISION_STOP( _ALLR, 1 );
        for( int i = 0; i<m; i++ ) {
          d->H1[i][m+b] = d->H1[m+b][i];
          d->H1[m+b][i] = conj_PRECISION(d->H1[m+b][i]);
        }
      } else {
        for( int i = 0; i<m; i++ ) {
          d->H1[i][m+b] = y[i];
          d->H1[m+b][i] = conj_PRECISION(d->H1[i][m+b]);
        }
      }
      d->H1[m+b][m+b] = global_inner_product_PRECISION( V[m+b], AV[m+b], 0, vs, l, threading );
    }
    d->size += blocksize;
  }
}


void davidson_step_PRECISION( int steps, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  int m;

  for( int k = 0; k<steps; k++) {
    m = d->size;

    if( d->arnoldi==0) {
      if( g.davidson_setup==1 ) {
        get_residual_PRECISION( d->res, d, l, threading );
        gamma5_PRECISION( d->res, d->res, l, threading );
      } else
        get_residual_PRECISION_on_D( d->res, d, l, threading );
    }
    else
      vector_PRECISION_copy( d->res, d->V[m-1], 0, l->inner_vector_size, l );
    
    if( d->initial_phase ) {
      smoother_PRECISION( d->V[m], NULL, d->res, 4, _NO_RES, _NO_SHIFT, l, threading );
//       smoother_PRECISION( d->res, NULL, d->V[m], 2, _NO_RES, _NO_SHIFT, l, threading );
//       smoother_PRECISION( d->V[m], NULL, d->res, 3, _NO_RES, _NO_SHIFT, l, threading );
    }
    else {
      vcycle_PRECISION( d->V[m], NULL, d->res, _NO_RES, l, threading );
    }
    PRECISION norm = global_norm_PRECISION( d->res, 0, l->inner_vector_size, l, threading );
    printf0("m = %d norm r = %+.16lf    ", m, norm );
    norm = global_norm_PRECISION( d->V[m], 0, l->inner_vector_size, l, threading );
    printf0("norm V[%d] = %+.16lf\n", m, norm );
    for( int i = 0; i<5; i++ ) 
      printf0("r[%d] = %+.16lf%+.16lfi\n", i, CSPLIT(d->res[i]));
    for( int i = 0; i<5; i++ ) 
      printf0("V[%d,%d] = %+.16lf%+.16lfi\n", m, i, CSPLIT(d->V[m][i]));
    if( g.davidson_setup==2 ) extend_subspace_on_D_PRECISION( &(d->V[m]), 1, d, l, threading );
    else extend_subspace_PRECISION( &(d->V[m]), 1, d, l, threading );
  }
  
//   for( int j = 0; j<l->num_eig_vect; j++ ) {
//     if( creal(d->V[j][0]) < 0 )
//       vector_PRECISION_scale( d->V[j], d->V[j], -1, 0, l->inner_vector_size, l );
//   }
//   FILE *fop;
//     
//   fop = fopen("V1x.m","w");
//   fprintf(fop,"V1 = [");
//   for ( int i = 0; i<l->num_eig_vect; i++ ) {
//     trans_back_PRECISION( d->W[i], d->V[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<l->inner_vector_size; j++ )
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(d->W[i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"V1=V1.';\n");
//   fclose(fop);
  
}


void restart_subspace_PRECISION( davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  int n = d->size, m = MIN(l->num_eig_vect,d->size), vs = l->inner_vector_size;
  char N = 'N', C = 'C', T = 'T'; // GEMM: "nothing", "complex conjugate", "transpose"
  complex_PRECISION one = 1, zero = 0, *s = d->s, v[n*n], H1[m*m], H2[m*m];
//   vector_PRECISION *V = d->V, *AV = d->AV, *W = d->W;
  
//   FILE *fop;
//   fop = fopen("H1H2.m","w");
//   fprintf(fop,"H1 = [\n");
//   for( int i = 0; i<n; i++ ) {
//     for( int j = 0; j<n; j++ ) {
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H1[i][j]));
//     }
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"H2 = [\n");
//   for( int i = 0; i<n; i++ ) {
//     for( int j = 0; j<n; j++ ) {
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H2[i][j]));
//     }
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fclose(fop);
//   printf0("\nrestart begin\n");
//   getchar();
  
  // get eigenvalues and orthogonalize eigenvectors
  dense_ev_PRECISION_driver( d->theta, v, d->H1, g.davidson_setup==2?NULL:d->H2, n ); // NOTE: Replace by zgges? (QZ algorithm)
  qr_PRECISION( s, v, n, n, d->maxsize );
  
//   fprintf(fop,"VV=[\n");
//   for (int i = 0; i<m; i++) {
//     for (int j = 0; j<m; j++) {
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(*(s+i*d->maxsize+j)));
//     }
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n VV = VV.';\n");
//   fprintf(fop,"theta = [\n");
//   for (int i = 0; i<m; i++) {
//     fprintf(fop,"%+.16lf%+.16lfi\n", CSPLIT(d->theta[i]));
//   }
//   fprintf(fop,"];\n");
//   
//   printf0("restart EV end\n");
//   fclose(fop);
//   getchar();

  
  // rebuild V from smallest ritz values 
  for ( int i=0; i<m; i++ ) {
    vector_PRECISION_define( d->W[i], 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( d->W[i], d->V, s+i*d->maxsize, 1, m, 0, vs, l );
  }

  // rebuild AV from smallest ritz values 
  for ( int i=0; i<m; i++ ) {
    vector_PRECISION_copy( d->V[i], d->W[i], 0, vs, l );
    vector_PRECISION_define( d->W[i], 0, 0, vs, l );
    vector_PRECISION_multi_saxpy( d->W[i], d->AV, s+i*d->maxsize, 1, m, 0, vs, l );
  }
  for ( int i=0; i<m; i++ ) {
    vector_PRECISION_copy( d->AV[i], d->W[i], 0, vs, l );
  }
//   printf0("V\n");
//   for( int j = 0; j<5; j++ ) {
//     printf0("%+.8lf%+.8lfi ", CSPLIT(V[0][j]));
//   }
//   printf0("\n");
//   getchar();
//   gram_schmidt_PRECISION( d->V, v, 0, m, l, threading );
//   test_eigenspace_PRECISION( d, l, threading );
  
  // update H1 and H2
  gemm_PRECISION( &C, &T, &m, &n, &n, &one, s, &(d->maxsize), d->H1[0], &(d->maxsize), &zero, v, &n ); // v <- s* H1
  gemm_PRECISION( &N, &N, &m, &m, &n, &one, v, &n, s, &(d->maxsize), &zero, H1, &m );                  // H1 <- v s
  gemm_PRECISION( &C, &T, &m, &n, &n, &one, s, &(d->maxsize), d->H2[0], &(d->maxsize), &zero, v, &n ); // v <- s* H2
  gemm_PRECISION( &N, &N, &m, &m, &n, &one, v, &n, s, &(d->maxsize), &zero, H2, &m );                  // H2 <- v s
  for( int i=0; i<m; i++) {
    for( int j=0; j<m; j++) {
      d->H1[j][i] = H1[j+m*i];
      d->H2[j][i] = H2[j+m*i];
    }
  }
  
//     FILE *fop;
//   printf0("vs=%d\n",vs);
//   fop = fopen("Vrestart.m","w");
//   fprintf(fop,"V_restart = [");
//   for ( int i = 0; i<20; i++ ) {
//     trans_back_PRECISION( d->V[20+i], d->V[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<vs; j++ )
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(d->V[20+i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"V_restart=V_restart.';\n");
//   
//   fprintf(fop,"AV_restart = [");
//   for ( int i = 0; i<20; i++ ) {
//     trans_back_PRECISION( d->AV[20+i], d->AV[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<vs; j++ )
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(d->AV[20+i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"AV_restart=AV_restart.';\n");
//   fclose(fop);
  
//   FILE *fop;
//   fop = fopen("H1H2.m","w");
//   fprintf(fop,"H1 = [\n");
//   for( int i = 0; i<m; i++ ) {
//     for( int j = 0; j<m; j++ ) {
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H1[i][j]));
//     }
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"H2 = [\n");
//   for( int i = 0; i<m; i++ ) {
//     for( int j = 0; j<m; j++ ) {
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->H2[i][j]));
//     }
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fclose(fop);
//   printf0("done\n");
//   getchar();
//   complex_PRECISION buffer[m];
//   gram_schmidt_PRECISION( d->V, buffer, 0, m, l, threading );

  d->size = m;
}


void initial_davidson_setup_PRECISION( int setup_iter, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  int vs = l->inner_vector_size, nev = l->num_eig_vect;
  
//   vector_PRECISION_define( d->V[0], 1, 0, vs, l );
  vector_PRECISION_define( d->V[0], 0, 0, vs, l );
  if (g.my_rank == 0) vector_PRECISION_define( d->V[0], 1, 0, 1, l );
  setup_iter = 2;
  
  if( g.davidson_setup==2 )
    extend_subspace_on_D_PRECISION( &(d->V[0]), 1, d, l, threading );
  else
    extend_subspace_PRECISION( &(d->V[0]), 1, d, l, threading );
  
  for( int k = 0; k<setup_iter; k++ ) {
    davidson_step_PRECISION( k==0?nev-1:nev, d, l, threading );
//     for( int i = 0; i<nev*setup_iter; i++ ) {
//       vector_PRECISION_define_random(l->d_PRECISION.V[i], 0, l->inner_vector_size, l );
//     }
    restart_subspace_PRECISION( d, l, threading );
    printf0("Level: %d -- Initial Setup Step %d/%d finished\n", l->level, k+1, setup_iter );
  }
//     restart_subspace_PRECISION( d, l, threading );
  
//   complex_PRECISION buffer[2*nev];
//   gram_schmidt_PRECISION( d->V, buffer, 0, nev, l, threading );
    
  for( int j = 0; j<nev; j++ ) {
    if( creal(d->V[j][0]) < 0 )
      vector_PRECISION_scale( d->V[j], d->V[j], -1, 0, l->inner_vector_size, l );
  }
  
  printf0("V_init\n");
  for( int i = 0; i<1; i++ ) {
    for( int j = 0; j<5; j++ ) {
      printf0("%+.16lf%+.16lfi ", CSPLIT(d->V[i][j]));
    }
    printf0("\n");
  }
//   getchar();

  define_next_level_PRECISION( d, l, threading );
//   test_eigenspace_PRECISION( d, l, threading );
//   {
//       int ivs = l->next_level->inner_vector_size, vs = l->next_level->vector_size, bvs = l->inner_vector_size;
//     double tmp;
//     vector_PRECISION ee = NULL, v = NULL;
//     
//     MALLOC( ee, complex_PRECISION, bvs );
//     MALLOC( v, complex_PRECISION, vs );
//     vector_PRECISION_define( ee, 1, 0, bvs, l );
//     vector_PRECISION_define( v, 0, 0, vs, l );
    
//     restrict_PRECISION( v, ee, l, threading );
//     tmp = global_norm_PRECISION( v, 0, ivs, l->next_level, threading );
//     printf0("\nnorm restriction on ones after init: %+.14lf\n", tmp);
//     for (int i=0; i<5; i++)
//       printf0("%+.16lf%+.16lfi ",CSPLIT(v[i]));
//     getchar();
//     FREE( ee, complex_PRECISION, bvs );
//     FREE( v, complex_PRECISION, vs );
//   }
  // print results in MATLAB format  
//     FILE *fop;
//     
//   fop = fopen("V_init.m","w");
//   fprintf(fop,"Vinit = [");
//   for ( int i = 0; i<nev*setup_iter; i++ ) {
//     trans_back_PRECISION( d->W[i], d->V[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<vs; j++ )
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(d->W[i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"Vinit=Vinit.';\n");
//   fclose(fop);
//     getchar();
  
  
  // print results in MATLAB format
//   FILE *fop;
  
  
  
//   fop = fopen("P_init.m","w");
//   fprintf(fop,"Pinit = [");
//   for ( int i = 0; i<nev; i++ ) {
//     trans_back_PRECISION( d->W[i], l->is_PRECISION.interpolation[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<vs; j++ )
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(d->W[i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"Pinit=Pinit.';\n");
//   fclose(fop);
//     getchar();

}


void hierarchy_update_PRECISION( vector_PRECISION *V, davidson_setup_PRECISION_struct *d, level_struct *l, Thread *threading ) {
  
  if( g.davidson_setup==2 ) 
    extend_subspace_on_D_PRECISION( V, d->size, d, l, threading );
  else
    extend_subspace_PRECISION( V, d->size, d, l, threading );
  
  restart_subspace_PRECISION( d, l, threading );
  define_next_level_PRECISION( d, l, threading );
  re_setup_PRECISION( l, threading );
  if( l->next_level->level > 0 ) {
    for( int i = 0; i<d->size; i++ )
      restrict_PRECISION( l->next_level->d_PRECISION.V[i], d->V[i], l, threading );
    l->next_level->d_PRECISION.size = d->size;
    hierarchy_update_PRECISION( l->next_level->d_PRECISION.V, &(l->next_level->d_PRECISION), l->next_level, threading );
  }
  test_eigenspace_PRECISION( d, l, threading );
}


void iterative_davidson_setup_PRECISION( int iterative_iter, level_struct *l, Thread *threading ) {
  
  int nev = l->num_eig_vect;
  
  l->d_PRECISION.initial_phase = 0;
  
  for( int i = 0; i < iterative_iter; i++ ) {
    davidson_step_PRECISION( nev, &(l->d_PRECISION), l, threading );
//     for( int i = 0; i<nev; i++ ) {
//       vector_PRECISION_define_random(l->d_PRECISION.V[i], 0, l->inner_vector_size, l );
//     }
  hierarchy_update_PRECISION( NULL, &(l->d_PRECISION), l, threading );
    
//   FILE *fop;
//     
//   fop = fopen("V_iter.m","w");
//   fprintf(fop,"Viter = [");
//   for ( int i = 0; i<nev; i++ ) {
//     trans_back_PRECISION( l->d_PRECISION.W[i], l->d_PRECISION.V[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<l->inner_vector_size; j++ )
//       fprintf(fop,"%+.16lf%+.16lfi ", CSPLIT(l->d_PRECISION.W[i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"Viter=Viter.';\n");
//   fclose(fop);
//     
//   fop = fopen("P_iter.m","w");
//   fprintf(fop,"Piter = [");

//   for ( int i = 0; i<nev; i++ ) {
//     trans_back_PRECISION( l->d_PRECISION.W[i], l->is_PRECISION.interpolation[i], l->s_PRECISION.op.translation_table, l, threading );
//     for ( int j = 0; j<l->inner_vector_size; j++ )
//       fprintf(fop,"%+.14lf%+.14lfi ", CSPLIT(l->d_PRECISION.W[i][j]));
//     fprintf(fop,";\n");
//   }
//   fprintf(fop,"];\n");
//   fprintf(fop,"Piter=Piter.';\n");
//   fclose(fop);
//     printf0("V_iter\n");
//     for( int i = 0; i<1; i++ ) {
//       for( int j = 0; j<5; j++ ) {
//         printf0("%+.16lf%+.16lfi ", CSPLIT(l->d_PRECISION.V[i][j]));
//       }
//       printf0("\n");
//     }
//     getchar();
//     if( l->next_level->level > 0 )
//       iterative_davidson_setup_PRECISION( iterative_iter, l->next_level, threading );
//     printf0("Level: %d -- Iterative Setup Step %d/%d finished\n", l->level, i+1, iterative_iter);
  }
}
