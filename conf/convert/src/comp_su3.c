#include <includes.h>

#define SQRT3R .577350269189625764509148780501957

#ifndef FORTRAN_UNDERSCORE
//#define cgeev_ cgeev
//#define zgeev_ zgeev
#endif


const complex el_zero = { 0.0, 0.0 };



/* Prototype of lapack function */
void zgeev_( char *, char *, int *, complex *, int *, complex *,
	     complex *, int *, complex *, int *, complex *, int *, double *, int * );
/* Variables used by lapack */
char el_JOBVL = 'V';
char el_JOBVR = 'V';
int el_N = 3;
int el_LDA = 3;
complex el_W[3];
complex el_VL[3 * 3];
int el_LDVL = 3;
complex el_VR[3 * 3];
int el_LDVR = 3;
complex el_WORK[16 * 3];
int el_LWORK = 16 * 3;
double el_RWORK[2 * 3];
int el_INFO;
 /**/
    /*  Returns the complex exponential of the complex number `x'  */
    complex zexp( complex x )
{
    complex c;
    double d = exp( x.real );
    c.real = d * cos( x.imag );
    c.imag = d * sin( x.imag );
    return c;
}


/*  Returns the phase of the complex number `x'.  -pi < carg(x) <= pi 
   Returns 0 if x==0.  */
double zarg( complex x )
{
    if( x.imag == 0.0 )
    {
	if( x.real >= 0.0 )
	    return 0.0;
	return PI;
    }
    if( x.real == 0.0 )
    {
	if( x.imag > 0.0 )
	    return PI / 2.0;
	return -PI / 2.0;
    }
    if( x.real > 0 )
	return atan( x.imag / x.real );
    if( x.imag > 0 )
	return atan( -x.real / x.imag ) + PI / 2.0;
    return atan( x.imag / x.real ) - PI;
}


/*  Returns the complex logarithm of `x'.  -pi < clog(x).imag <= pi  
   Returns 0 if x==0.  */
complex zlog( complex x )
{
    complex c = { 0.0, 0.0 };
    double d = CABS( x );
    if( d == 0.0 )
	return c;
    c.real = log( d ) / 2.0;
    c.imag = zarg( x );
    return c;
}





/* Calculates the matrix logarithm of `A' and puts the result back into `A' */
void mlog( complex * A )
{
    complex tmp1, tmp2, tmp3;
    int i, j, k;

    /* Transposing `A' because of lapack's fortran convention */
    tmp1 = A[1];
    A[1] = A[3];
    A[3] = tmp1;
    tmp1 = A[2];
    A[2] = A[6];
    A[6] = tmp1;
    tmp1 = A[5];
    A[5] = A[7];
    A[7] = tmp1;

    /* Calling lapack to calculate the eigenvalues and eigenvectors */
    el_INFO = 1;
    zgeev_( &el_JOBVL, &el_JOBVR, &el_N, A, &el_LDA, el_W, el_VL, &el_LDVL, el_VR, &el_LDVR, el_WORK, &el_LWORK, el_RWORK, &el_INFO );

    /* Normalizing the eigenvectors */
    for ( i = 0; i < 3; i++ )
    {
	tmp1 = el_zero;
	for ( j = 0; j < 3; j++ )
	{
	    CONJG( el_VL[3 * i + j], tmp2 );
	    CMUL( tmp2, el_VR[3 * i + j], tmp3 );
	    CADD( tmp3, tmp1, tmp1 );
	}
	for ( j = 0; j < 3; j++ )
	{
	    CDIV( el_VR[3 * i + j], tmp1, tmp2 );
	    el_VR[3 * i + j] = tmp2;
	}
    }

    /* Calculating the result */
    /* Zeroing `A' */
    for ( i = 0; i < 9; i++ )
	A[i] = el_zero;

    /* Calculating the logarithm of the eigenvalues */
    for ( i = 0; i < 3; i++ )
	el_W[i] = zlog( el_W[i] );

    /*  Making '-pi < Im(Tr(log(A))) <= pi' true */
    while( el_W[0].imag + el_W[1].imag + el_W[2].imag > PI )
	el_W[2].imag -= 2.0 * PI;
    while( el_W[0].imag + el_W[1].imag + el_W[2].imag <= -PI )
	el_W[2].imag += 2.0 * PI;

    /* Loop over the eigenvalues */
    for ( i = 0; i < 3; i++ )
    {
	/* Adding to `A' the projector corresponding to the eigenvalue */
	for ( j = 0; j < 3; j++ )
	    for ( k = 0; k < 3; k++ )
	    {
		CONJG( el_VL[3 * i + k], tmp1 );
		CMUL( el_VR[3 * i + j], tmp1, tmp2 );
		CMUL( tmp2, el_W[i], tmp1 );
		CADD( tmp1, A[3 * j + k], A[3 * j + k] );
	    }
    }
    /* log(A) is ready */

}




/* Calculates the matrix exponential of `A' and puts the result back into `A' */
void mexp( complex * A )
{
    complex tmp1, tmp2, tmp3;
    int i, j, k;

    /* Transposing `A' because of lapack's fortran convention */
    tmp1 = A[1];
    A[1] = A[3];
    A[3] = tmp1;
    tmp1 = A[2];
    A[2] = A[6];
    A[6] = tmp1;
    tmp1 = A[5];
    A[5] = A[7];
    A[7] = tmp1;

    /* Calling lapack to calculate the eigenvalues and eigenvectors */
    el_INFO = 1;
    zgeev_( &el_JOBVL, &el_JOBVR, &el_N, A, &el_LDA, el_W, el_VL, &el_LDVL, el_VR, &el_LDVR, el_WORK, &el_LWORK, el_RWORK, &el_INFO );

    /* Normalizing the eigenvectors */
    for ( i = 0; i < 3; i++ )
    {
	tmp1 = el_zero;
	for ( j = 0; j < 3; j++ )
	{
	    CONJG( el_VL[3 * i + j], tmp2 );
	    CMUL( tmp2, el_VR[3 * i + j], tmp3 );
	    CADD( tmp3, tmp1, tmp1 );
	}
	for ( j = 0; j < 3; j++ )
	{
	    CDIV( el_VR[3 * i + j], tmp1, tmp2 );
	    el_VR[3 * i + j] = tmp2;
	}
    }

    /* Calculating the result */
    /* Zeroing `A' */
    for ( i = 0; i < 9; i++ )
	A[i] = el_zero;

    /* Loop over the eigenvalues */
    for ( i = 0; i < 3; i++ )
    {
	/* Calculating the exponential of the eigenvalue */
	el_W[i] = zexp( el_W[i] );
	/* Adding to `A' the projector corresponding to the eigenvalue */
	for ( j = 0; j < 3; j++ )
	    for ( k = 0; k < 3; k++ )
	    {
		CONJG( el_VL[3 * i + k], tmp1 );
		CMUL( el_VR[3 * i + j], tmp1, tmp2 );
		CMUL( tmp2, el_W[i], tmp1 );
		CADD( tmp1, A[3 * j + k], A[3 * j + k] );
	    }
    }
    /* exp(A) is ready */

}




/*  Returns Tr( A * lambda_n )  */
complex dot_gellmann( complex * A, int n )
{
    complex c;

    switch ( n )
    {
	case 1:
	    c.real = A[1].real + A[3].real;
	    c.imag = A[1].imag + A[3].imag;
	    break;
	case 2:
	    c.real = -A[1].imag + A[3].imag;
	    c.imag = A[1].real - A[3].real;
	    break;
	case 3:
	    c.real = A[0].real - A[4].real;
	    c.imag = A[0].imag - A[4].imag;
	    break;
	case 4:
	    c.real = A[2].real + A[6].real;
	    c.imag = A[2].imag + A[6].imag;
	    break;
	case 5:
	    c.real = -A[2].imag + A[6].imag;
	    c.imag = A[2].real - A[6].real;
	    break;
	case 6:
	    c.real = A[5].real + A[7].real;
	    c.imag = A[5].imag + A[7].imag;
	    break;
	case 7:
	    c.real = -A[5].imag + A[7].imag;
	    c.imag = A[5].real - A[7].real;
	    break;
	case 8:
	    c.real = ( A[0].real + A[4].real - A[8].real - A[8].real ) * SQRT3R;
	    c.imag = ( A[0].imag + A[4].imag - A[8].imag - A[8].imag ) * SQRT3R;
	    break;
	default:
	    c.real = c.imag = 0.0;
    }
    return c;
}




void su3_to_comp( su3_matrix * U, su3_matrix_comp * alpha )
{
    complex A[9];
    int i, j;

    for ( i = 0; i < 3; i++ )
	for ( j = 0; j < 3; j++ )
	{
	    A[3 * i + j].real = ( double ) ( ( *U ).ROWCOL( i, j ).real );
	    A[3 * i + j].imag = ( double ) ( ( *U ).ROWCOL( i, j ).imag );
	}

    mlog( A );

    for ( i = 0; i < 8; i++ )
    {
	alpha->a[i] = ( float ) ( ( dot_gellmann( A, i + 1 ) ).imag );
    }
}



void comp_to_su3( su3_matrix_comp * alpha, su3_matrix * result )
{
    complex A[9];
    int i, j;

    A[0].real = 0.0;
    A[0].imag = ( ( double ) alpha->a[7] * SQRT3R + ( double ) alpha->a[2] ) / 2.0;
    A[1].real = ( double ) alpha->a[1] / 2.0;
    A[1].imag = ( double ) alpha->a[0] / 2.0;
    A[2].real = ( double ) alpha->a[4] / 2.0;
    A[2].imag = ( double ) alpha->a[3] / 2.0;
    A[3].real = -( double ) alpha->a[1] / 2.0;
    A[3].imag = ( double ) alpha->a[0] / 2.0;
    A[4].real = 0.0;
    A[4].imag = ( -( double ) alpha->a[2] + ( double ) alpha->a[7] * SQRT3R ) / 2.0;
    A[5].real = ( double ) alpha->a[6] / 2.0;
    A[5].imag = ( double ) alpha->a[5] / 2.0;
    A[6].real = -( double ) alpha->a[4] / 2.0;
    A[6].imag = ( double ) alpha->a[3] / 2.0;
    A[7].real = -( double ) alpha->a[6] / 2.0;
    A[7].imag = ( double ) alpha->a[5] / 2.0;
    A[8].real = 0.0;
    A[8].imag = -( double ) alpha->a[7] * SQRT3R;

    mexp( A );

    for ( i = 0; i < 3; i++ )
	for ( j = 0; j < 3; j++ )
	{
	    ( *result ).ROWCOL( i, j ).real = ( double_fast_t ) A[3 * i + j].real;
	    ( *result ).ROWCOL( i, j ).imag = ( double_fast_t ) A[3 * i + j].imag;
	}
}
