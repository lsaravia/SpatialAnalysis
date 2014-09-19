// This program calculates the two-dimensional fast 
// fourier transform by combining the two subprograms 
// fastf and ft2d. channels used are as follows: 
//   ft01 - input matrix size and format statements 
//   ft02 - input data matrix 
//   ft03 - complete spectrum  i=0,...,n1-1 ; j=0,...,n2-1 
//   ft04 - x fourier coefficients to enable inversion 
//   ft05 - y fourier coefficients to enable inversion 
//   ft06 - output matrix size and variance (i.e. ss/n1*n2) 
#include "Spectral.h"
//#include "fortify.h"
#include <math.h>

using namespace std;

int fastf_(double *xreal, double *ximag, int* isize);
int ft2d_(double *a, double *b, double *work1, double *work2, int *n1, int *n2);

double SpectralAnalysis::Spec2D(simplmat<double>& r,simplmat<double>& xout)
{
    int mtot;
    double total,mean,var; // percent;

	int m1 = r.getRows();
	int m2 = r.getCols();
	int n2 = m2;
	int n1 = m1;
//    int mm1 = m1 / 2 + 1;
//    int mm2 = m2 / 2 + 1;
	
	xout.resize(m1,m2);

    double * x = new double[m1*m2];
    double * y = new double[m1*m2];
    double * work1 = new double[m2];
    double * work2 = new double[m2];

	if( x==NULL )
	{
		cerr << "Error memory allocation, variable x" << endl;
		return -1;
	}

	int i,j,l;
    for (j=0; j< m2; ++j)
    {
		for (i=0; i< m1; ++i)
		{
		    l = j * m1 + i;
		    x[l] = r(i,j);
		}
    }
    
//
//          Subtract mean and calculate variance
//
    total = 0.0;
    mtot = m1 * m2;

    for (i=0; i<mtot; ++i)
		total += x[i];

    mean = total / double(mtot);
    
    for (i=0; i<mtot; ++i)
    {
		x[i] -= mean;
		y[i] = 0.0;
    }

    var = 0.0;
    
    for (i=0; i < mtot; ++i)
		var += x[i] * x[i];

    var /= (double) mtot;

    ft2d_(x, y, work1, work2, &n1, &n2);
    
//
//          Rearrange in standard form and output
//
	
    for (j = 0; j <m2; ++j)
    {
		for (i= 0; i< m1; ++i)
		{
		    l = j * m1 + i;
		    xout(i,j) = mtot * (x[l] * x[l] + y[l] * y[l]);
		}
    }
    
    delete[] x;
    delete[] y;
    delete[] work1;
	delete[] work2;

	return var;
}


int fastf_(double *xreal, double *ximag, int* isize)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12, i__13, i__14, i__15, i__16, i__17, i__18, i__19, 
	    i__20, i__21, i__22, i__23, i__24, i__25, i__26, i__27, i__28, 
	    i__29, i__30, i__31, i__32, i__33, i__34, i__35, i__36, i__37, 
	    i__38, i__39;
    double d__1;
    static int equiv_19[20];

    /* Local variables */
    static double bcos, bsin;
    static int ifaca, ifcab, ifacc, k;
#define l (equiv_19)
    static int n;
    static double z__;
    static int litla, itime;
    static double tempr;
    static int i0, i1, i2, i3;
#define l1 (equiv_19)
#define l2 (equiv_19 + 1)
#define l3 (equiv_19 + 2)
#define l4 (equiv_19 + 3)
#define l5 (equiv_19 + 4)
#define l6 (equiv_19 + 5)
#define l7 (equiv_19 + 6)
#define l8 (equiv_19 + 7)
#define l9 (equiv_19 + 8)
    static int j1, j2, j3;
    static double x1, x2, x3, y1, y2, y3;
    static int j4, j5, j6, j7, j8, j9, j10, j11;
#define l10 (equiv_19 + 9)
#define l11 (equiv_19 + 10)
#define l12 (equiv_19 + 11)
#define l13 (equiv_19 + 12)
#define l14 (equiv_19 + 13)
#define l15 (equiv_19 + 14)
#define l16 (equiv_19 + 15)
#define l17 (equiv_19 + 16)
#define l18 (equiv_19 + 17)
#define l19 (equiv_19 + 18)
#define l20 (equiv_19 + 19)
    static int ii, j12, j13, j14, j15, j16, j17, j18, j19, j20;
    static double cw1;
    static float cw2, cw3;
    static double sw1, xs0, xs1, xs2, xs3, ys0, ys1, ys2, ys3;
    static float sw2, sw3;
    static double pie;


/* RADIX 4 COMPLEX DISCRETE FAST FOURIER TRANSFORM */

/* XREAL = ARRAY WHICH ON INPUT CONTAINS REAL PART OF DATA */
/* FOR TRANSFORMATION AND ON OUTPUT GIVES REAL PART OF THE */
/* RESULT,TYPE REAL, DIMENSION IABS(ISIZE) OR GREATER. */
/* XIMAG = ARRAY WHICH ON INPUT CONTAINS IMAGINARY PART */
/* FOR TRANSFORMATION AND ON OUTPUT GIVES IMAGINARY PART OF */
/* THE RESULT, TYPE REAL, DIMENSION IABS(ISIZE) OR GREATER. */
/* ISIZE = INTEGER VARIABLE OR CONSTANT SPECIFYING LENGTH AND */
/* TYPE OF TRANSFORM. THE LENGTH IS IABS(ISIZE) WHICH MUST BE */
/* A POWER OF 2. MINIMUM LENGTH 4. MAXIMUM 2**20. IF ISIZE IS */
/* POSITIVE THE FORWARD TRANSFORM IS CALCULATED. IF NEGATIVE */
/* THE INVERSE TRANSFORM IS FOUND. */

    /* Parameter adjustments */
    --ximag;
    --xreal;

    /* Function Body */
    pie = (float)3.141592653589793;
    n = abs(*isize);
    if (n - 4 >= 0) {
	goto L1;
    } else {
	goto L24;
    }
/* SET UP INITIAL VALUES OF TRANSFORM SPLIT */
L1:
    ifacc = 1;
    ifaca = n / 4;
    if (*isize >= 0) {
	goto L4;
    } else {
	goto L2;
    }
/* IF THIS IS TO BE AN INVERSE TRANSFORM, CONJUGATE THE DATA */
L2:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L3: */
	ximag[k] = -ximag[k];
    }
L4:
    itime = 0;
L5:
    ifcab = ifaca << 2;
    itime += 2;
/* DO THE TRANSFORMS REQUIRED BY THIS STAGE */
    z__ = pie / (double) ifcab;
/* Computing 2nd power */
    d__1 = sin(z__);
    bcos = d__1 * d__1 * (float)-2.;
    bsin = sin(z__ * (float)2.);
    cw1 = (float)1.;
    sw1 = (float)0.;
    i__1 = ifaca;
    for (litla = 1; litla <= i__1; ++litla) {
	i__2 = n;
	i__3 = ifcab;
	for (i0 = litla; i__3 < 0 ? i0 >= i__2 : i0 <= i__2; i0 += i__3) {
/* THIS IS THE MAIN CALCULATION OF RADIX 4 TRANSFORMS */
	    i1 = i0 + ifaca;
	    i2 = i1 + ifaca;
	    i3 = i2 + ifaca;
	    xs0 = xreal[i0] + xreal[i2];
	    xs1 = xreal[i0] - xreal[i2];
	    ys0 = ximag[i0] + ximag[i2];
	    ys1 = ximag[i0] - ximag[i2];
	    xs2 = xreal[i1] + xreal[i3];
	    xs3 = xreal[i1] - xreal[i3];
	    ys2 = ximag[i1] + ximag[i3];
	    ys3 = ximag[i1] - ximag[i3];
	    xreal[i0] = xs0 + xs2;
	    ximag[i0] = ys0 + ys2;
	    x1 = xs1 + ys3;
	    y1 = ys1 - xs3;
	    x2 = xs0 - xs2;
	    y2 = ys0 - ys2;
	    x3 = xs1 - ys3;
	    y3 = ys1 + xs3;
	    if ((i__4 = litla - 1) < 0) {
		goto L24;
	    } else if (i__4 == 0) {
		goto L6;
	    } else {
		goto L7;
	    }
L6:
	    xreal[i2] = x1;
	    ximag[i2] = y1;
	    xreal[i1] = x2;
	    ximag[i1] = y2;
	    xreal[i3] = x3;
	    ximag[i3] = y3;
	    goto L8;
/* MULTIPLY BY TWIDDLE FACTORS IF REQUIRED */
L7:
	    xreal[i2] = x1 * cw1 + y1 * sw1;
	    ximag[i2] = y1 * cw1 - x1 * sw1;
	    xreal[i1] = x2 * cw2 + y2 * sw2;
	    ximag[i1] = y2 * cw2 - x2 * sw2;
	    xreal[i3] = x3 * cw3 + y3 * sw3;
	    ximag[i3] = y3 * cw3 - x3 * sw3;
L8:
	    ;
	}
	if ((i__3 = litla - ifaca) < 0) {
	    goto L9;
	} else if (i__3 == 0) {
	    goto L10;
	} else {
	    goto L24;
	}
/* CALCULATE A NEW SET OF TWIDDLE FACTORS */
L9:
	z__ = cw1 * bcos - sw1 * bsin + cw1;
	sw1 = bcos * sw1 + bsin * cw1 + sw1;
	tempr = (float)1.5 - (z__ * z__ + sw1 * sw1) * (float).5;
	cw1 = z__ * tempr;
	sw1 *= tempr;
	cw2 = cw1 * cw1 - sw1 * sw1;
	sw2 = cw1 * (float)2. * sw1;
	cw3 = cw1 * cw2 - sw1 * sw2;
	sw3 = cw1 * sw2 + cw2 * sw1;
L10:
	;
    }
    if (ifaca - 1 <= 0) {
	goto L14;
    } else {
	goto L11;
    }
/* SET UP THE TRANSFORM SPLIT FOR THE NEXT STAGE */
L11:
    ifacc <<= 2;
    ifaca /= 4;
    if (ifaca < 0) {
	goto L24;
    } else if (ifaca == 0) {
	goto L12;
    } else {
	goto L5;
    }
/* THIS IS THE CALCULATION OF A RADIX TWO STAGE */
L12:
    i__1 = n;
    for (k = 1; k <= i__1; k += 2) {
	tempr = xreal[k] + xreal[k + 1];
	xreal[k + 1] = xreal[k] - xreal[k + 1];
	xreal[k] = tempr;
	tempr = ximag[k] + ximag[k + 1];
	ximag[k + 1] = ximag[k] - ximag[k + 1];
/* L13: */
	ximag[k] = tempr;
    }
    ++itime;
L14:
    if (*isize < 0) {
	goto L15;
    } else if (*isize == 0) {
	goto L19;
    } else {
	goto L17;
    }
/* IF THIS WAS AN INVERSE TRANSFORM, CONJUGATE THE RESULT */
L15:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* L16: */
	ximag[k] = -ximag[k];
    }
    goto L19;
/* IF THIS WAS A FORWARD TRANSFORM, SCALE THE RESULT */
L17:
    z__ = (float)1. / (double) n;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	xreal[k] *= z__;
/* L18: */
	ximag[k] *= z__;
    }
/* UNSCRAMBLE THE RESULT */
L19:
    i1 = 20 - itime;
    i__1 = i1;
    for (k = 1; k <= i__1; ++k) {
/* L20: */
	l[k - 1] = 1;
    }
    ii = 1;
    ++i1;
    for (k = i1; k <= 20; ++k) {
	ii <<= 1;
/* L21: */
	l[k - 1] = ii;
    }
    ii = 1;
    i__1 = *l1;
    for (j1 = 1; j1 <= i__1; ++j1) {
	i__3 = *l2;
	i__2 = *l1;
	for (j2 = j1; i__2 < 0 ? j2 >= i__3 : j2 <= i__3; j2 += i__2) {
	    i__4 = *l3;
	    i__5 = *l2;
	    for (j3 = j2; i__5 < 0 ? j3 >= i__4 : j3 <= i__4; j3 += i__5) {
		i__6 = *l4;
		i__7 = *l3;
		for (j4 = j3; i__7 < 0 ? j4 >= i__6 : j4 <= i__6; j4 += i__7) 
			{
		    i__8 = *l5;
		    i__9 = *l4;
		    for (j5 = j4; i__9 < 0 ? j5 >= i__8 : j5 <= i__8; j5 += 
			    i__9) {
			i__10 = *l6;
			i__11 = *l5;
			for (j6 = j5; i__11 < 0 ? j6 >= i__10 : j6 <= i__10; 
				j6 += i__11) {
			    i__12 = *l7;
			    i__13 = *l6;
			    for (j7 = j6; i__13 < 0 ? j7 >= i__12 : j7 <= 
				    i__12; j7 += i__13) {
				i__14 = *l8;
				i__15 = *l7;
				for (j8 = j7; i__15 < 0 ? j8 >= i__14 : j8 <= 
					i__14; j8 += i__15) {
				    i__16 = *l9;
				    i__17 = *l8;
				    for (j9 = j8; i__17 < 0 ? j9 >= i__16 : 
					    j9 <= i__16; j9 += i__17) {
					i__18 = *l10;
					i__19 = *l9;
					for (j10 = j9; i__19 < 0 ? j10 >= 
						i__18 : j10 <= i__18; j10 += 
						i__19) {
					    i__20 = *l11;
					    i__21 = *l10;
					    for (j11 = j10; i__21 < 0 ? j11 >=
						     i__20 : j11 <= i__20; 
						    j11 += i__21) {
			  i__22 = *l12;
			  i__23 = *l11;
			  for (j12 = j11; i__23 < 0 ? j12 >= i__22 : j12 <= 
				  i__22; j12 += i__23) {
			      i__24 = *l13;
			      i__25 = *l12;
			      for (j13 = j12; i__25 < 0 ? j13 >= i__24 : j13 
				      <= i__24; j13 += i__25) {
				  i__26 = *l14;
				  i__27 = *l13;
				  for (j14 = j13; i__27 < 0 ? j14 >= i__26 : 
					  j14 <= i__26; j14 += i__27) {
				      i__28 = *l15;
				      i__29 = *l14;
				      for (j15 = j14; i__29 < 0 ? j15 >= 
					      i__28 : j15 <= i__28; j15 += 
					      i__29) {
					  i__30 = *l16;
					  i__31 = *l15;
					  for (j16 = j15; i__31 < 0 ? j16 >= 
						  i__30 : j16 <= i__30; j16 +=
						   i__31) {
			i__32 = *l17;
			i__33 = *l16;
			for (j17 = j16; i__33 < 0 ? j17 >= i__32 : j17 <= 
				i__32; j17 += i__33) {
			    i__34 = *l18;
			    i__35 = *l17;
			    for (j18 = j17; i__35 < 0 ? j18 >= i__34 : j18 <= 
				    i__34; j18 += i__35) {
				i__36 = *l19;
				i__37 = *l18;
				for (j19 = j18; i__37 < 0 ? j19 >= i__36 : 
					j19 <= i__36; j19 += i__37) {
				    i__38 = *l20;
				    i__39 = *l19;
				    for (j20 = j19; i__39 < 0 ? j20 >= i__38 :
					     j20 <= i__38; j20 += i__39) {
					if (ii - j20 >= 0) {
		      goto L23;
					} else {
		      goto L22;
					}
L22:
					tempr = xreal[ii];
					xreal[ii] = xreal[j20];
					xreal[j20] = tempr;
					tempr = ximag[ii];
					ximag[ii] = ximag[j20];
					ximag[j20] = tempr;
L23:
					++ii;
				    }
				}
			    }
			}
					  }
				      }
				  }
			      }
			  }
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
L24:
    return 0;
} /* fastf_ */

#undef l20
#undef l19
#undef l18
#undef l17
#undef l16
#undef l15
#undef l14
#undef l13
#undef l12
#undef l11
#undef l10
#undef l9
#undef l8
#undef l7
#undef l6
#undef l5
#undef l4
#undef l3
#undef l2
#undef l1
#undef l


int ft2d_(double *a, double *b, double *work1, double *work2, int *n1, int *n2)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int k, l;
    static int m1, m2, kk;


/* COMPLEX DISCRETE FAST FOURIER TRANSFORM IN TWO DIMENSIONS */

/* A = MATRIX REPRESENTING REAL PART OF TRANSFORM DATA, TYPE */
/* REAL, DIMENSIONS EXACTLY IABS(N1) X IABS(N2) */
/* B = MATRIX OF IMAGINARY PARTS OF TRANSFORMS DATA,TYPE */
/* REAL, DIMENSIONS EXACTLY IABS(N1) X IABS(N2) */
/* WORK1,WORK2 = WORKING ARRAYS IN ONE DIMENSION, TYPE REAL. */
/* BOTH DIMENSIONED IABS(N2) OR GREATER. */
/* N1,N2 = INTEGER VARIABLES OR CONSTANTS WHICH SPECIFY THE */
/* SIZE AND TYPE OF TRANSFORM. SIZE IS IABS(N1) X IABS(N2). */
/* BOTH MUST BE POWERS OF 2. MINIMUM 4. MAXIMUM 2**20. IF */
/* BOTH POSITIVE THE FORWARD TRANSFORM IS FOUND. IF BOTH */
/* NEGATIVE THE INVERSE. IF SIGNS ARE MIXED. SO IS TRANSFORM. */

/* THE REAL AND IMAGINARY RESULTS ARE RETURNED IN A AND B */

    /* Parameter adjustments */
    --work2;
    --work1;
    --b;
    --a;

    /* Function Body */
    m1 = abs(*n1);
    m2 = abs(*n2);
/* DO THE ROWS IN WORKING ARRAYS */
    i__1 = m1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = m2;
	for (l = 1; l <= i__2; ++l) {
	    kk = k + m1 * (l - 1);
	    work1[l] = a[kk];
/* L1: */
	    work2[l] = b[kk];
	}
	fastf_(&work1[1], &work2[1], n2);
	i__2 = m2;
	for (l = 1; l <= i__2; ++l) {
	    kk = k + m1 * (l - 1);
	    a[kk] = work1[l];
/* L2: */
	    b[kk] = work2[l];
	}
    }
/* DO THE COLUMNS IN PLACE */
    i__2 = m2;
    for (k = 1; k <= i__2; ++k) {
	kk = m1 * (k - 1) + 1;
/* L3: */
	fastf_(&a[kk], &b[kk], n1);
    }
    return 0;
} /* ft2d_ */

