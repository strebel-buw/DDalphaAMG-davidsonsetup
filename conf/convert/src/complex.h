#ifndef _COMPLEX_H
#define _COMPLEX_H

#define COLORSPINOR(i,j)	d[j].c[i]
#define ROWCOL(i,j)	e[j][i]
#define double_hmc_t double
#define double_fast_t double

typedef struct
{
    double real;
    double imag;
} complex;

typedef struct
{
    double_hmc_t real;
    double_hmc_t imag;
} complex_hmc_t;

typedef struct
{
    double_fast_t real;
    double_fast_t imag;
} complex_fast_t;

/* Function Prototypes for Complex Numbers */
complex cmplx( double x, double y );
complex_fast_t cmplx_fast( double x, double y );
complex cadd( complex * a, complex * b );
complex cmul( complex * a, complex * b );
complex csub( complex * a, complex * b );
complex cdiv( complex * a, complex * b );
complex conjg( complex * a );
complex cexp_milc( complex * a );
complex clog_milc( complex * a );
complex csqrt_milc( complex * z );
complex ce_itheta( double theta );

/* Macros for Complex Numbers */
#define set_complex_equal(a,b) { (*b).real=(*a).real; (*b).imag=(*a).imag; }
/*    |*a|    */
#define cabs(a) (sqrt( (*a).real*(*a).real + (*a).imag*(*a).imag ) )
/*  *a * *a*  */
#define dcabs cabs
#define cabs_sq(a) ( (*a).real*(*a).real + (*a).imag*(*a).imag )
/* phase(*a)  */
#define carg(a) (atan2((double)(*a).imag, (double)(*a).real ) )
/*   b = a*   */
#define dcarg carg
#define CONJG(a,b) { (b).real = (a).real; (b).imag = -(a).imag; }
/*  c = a + b */
#define CADD(a,b,c) { (c).real = (a).real + (b).real;  \
	(c).imag = (a).imag + (b).imag; }
	/*  a += b    */
#define CSUM(a,b) { (a).real += (b).real; (a).imag += (b).imag; }
	/*  c = a - b */
#define CSUB(a,b,c) { (c).real = (a).real - (b).real;  \
	(c).imag = (a).imag - (b).imag; }
	/*  c = a * b */
#define CMUL(a,b,c) { (c).real = (a).real*(b).real - (a).imag*(b).imag; \
	(c).imag = (a).real*(b).imag + (a).imag*(b).real; }
	/* c = a / b  */
#define CDIV(a,b,c) { double t = (b).real*(b).real + (b).imag*(b).imag; \
	(c).real = ((a).real*(b).real + (a).imag*(b).imag)/t; \
		(c).imag = ((a).imag*(b).real - (a).real*(b).imag)/t; }
		/* c = a * b* */
#define CMUL_J(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
	(c).imag = (a).imag*(b).real - (a).real*(b).imag; }
	/* c = a* * b */
#define CMULJ_(a,b,c) { (c).real = (a).real*(b).real + (a).imag*(b).imag; \
	(c).imag = (a).real*(b).imag - (a).imag*(b).real; }
	/* c = (a*b)* */
#define CMULJJ(a,b,c) { (c).real =  (a).real*(b).real - (a).imag*(b).imag; \
	(c).imag = -(a).real*(b).imag - (a).imag*(b).real; }
	/* b = - a    */
#define CNEGATE(a,b) { (b).real = -(a).real; (b).imag = -(a).imag; }
	/* b =  ia    */
#define CMUL_I(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
	/* b = -ia    */
#define CMUL_MINUS_I(a,b) { (b).real = (a).imag; (b).imag = -(a).real; }
	/* c = ba     */
#define CMULREAL(a,b,c) { (c).real = (b) * (a).real; (c).imag = (b)*(a).imag; }
	/* c = a/b    */
#define CDIVREAL(a,b,c) { (c).real = (a).real/(b); (c).imag = (a).imag/(b); }

	/* a += i*b */
#define CSUM_TPI(a,b) { (a).real -= (b).imag; (a).imag +=  (b).real; }

	/* a += -i*b */
#define CSUM_TMI(a,b) { (a).real += (b).imag; (a).imag -=  (b).real; }

#define CABSSQ(c) (((c).real*(c).real+(c).imag*(c).imag))
#define CABS(c) (sqrt((c).real*(c).real+(c).imag*(c).imag))
/* c += a*b */
#define CMULSUM(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag + (a).imag*(b).real; }
/*  c += conj(a)*b */
#define CMULJ_SUM(a,b,c) { (c).real += (a).real*(b).real + (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag - (a).imag*(b).real; }
#endif
