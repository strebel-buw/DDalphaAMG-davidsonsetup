#ifndef _SU3_H
#define _SU3_H
#include "complex.h"

/* SU(3) */
typedef struct
{
    complex_fast_t e[3][3];
} su3_matrix;
typedef struct
{
    float a[8];
} su3_matrix_comp;
typedef struct
{
    complex_fast_t e[3][3][3][3];
} su3_hypermatrix;
typedef struct
{
    complex_fast_t c[3];
} su3_vector;
typedef struct
{
    complex_fast_t m01, m02, m12;
    double_fast_t m00im, m11im, m22im;
} anti_hermitmat;

typedef struct
{
    su3_vector d[4];
} wilson_vector;

typedef struct
{
    su3_vector h[2];
} half_wilson_vector;
typedef struct
{
    wilson_vector d[4];
} spin_wilson_vector;
typedef struct
{
    spin_wilson_vector c[3];
} wilson_propagator;
typedef struct
{
    wilson_vector c[3];
} color_wilson_vector;
typedef struct
{
    color_wilson_vector d[4];
} wilson_matrix;

/* Clover vectors */
typedef struct
{
    complex_fast_t tr[2][15];
} triangular;
typedef struct
{
    double_fast_t di[2][6];
} diagonal;
typedef struct
{
    complex_fast_t b[2][6];
} wilson_block_vector;

/* SU2 for gauge fixing */
typedef struct
{
    complex_fast_t esu2[2][2];
} su2_matrix;
typedef struct
{
    double_fast_t a[4];
} su2_matr_comp;

/* HMC evolution */
typedef struct
{
    complex_hmc_t e[3][3];
} hmc_su3_matrix;
typedef struct
{
    complex_hmc_t m01, m02, m12;
    double_hmc_t m00im, m11im, m22im;
} hmc_anti_hermitmat;

#define GAMMAFIVE -1		/* some integer which is not a direction */
#define PLUS 1			/* flags for selecting M or M_adjoint */
#define MINUS -1
/* Macros to multiply complex numbers by +-1 and +-i */
#define TIMESPLUSONE(a,b) { (b).real =  (a).real; (b).imag = (a).imag; }
#define TIMESMINUSONE(a,b) { (b).real =  -(a).real; (b).imag = -(a).imag; }
#define TIMESPLUSI(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
#define TIMESMINUSI(a,b) { (b).real =  (a).imag; (b).imag = -(a).real; }

#define FORMAT(a,b) for(a=0;a<3;a++) for (b=0;b<3;b++)

/* for lattice i/o */
void su3_to_comp( su3_matrix * U, su3_matrix_comp * alpha );
void comp_to_su3( su3_matrix_comp * alpha, su3_matrix * result );


#endif
