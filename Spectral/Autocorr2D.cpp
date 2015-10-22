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
#include "fortify.h"

int fastf_(double *xreal, double *ximag, int* isize);
int ft2d_(double *a, double *b, double *work1, double *work2, int *n1, int *n2);

double SpectralAnalysis::AutoFFT(simplmat<double>& r,simplmat<double>& xout)
{
    int mtot;
    double total,mean,var, percent;

	int m1 = r.getRows();
	int m2 = r.getCols();
    int mm1 = m1*2;
    int mm2 = m2*2;
	int n2 = mm2;
	int n1 = mm1;
	
	xout.resize(m1,m2);

    double * x = new double[mm1*mm2];
    double * y = new double[mm1*mm2];
    double * work1 = new double[mm2];
    double * work2 = new double[mm2];

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
		    l = j * mm1 + i;
		    x[l] = r(i,j);
		}
    }


    mtot = m1 * m2;
    mtot4 = mtot*4;
    for (i=mtot; i<mtot4; ++i)
		x[i]=0;


//
//          Subtract mean and calculate variance
//
    total = 0.0;

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
//  Derive inverse fft of periodogram
//
    for (i=0; i<mtot4; ++i)
    {
		x[i]= x[i]*x[i] + y[i]*y[i];
		y[i]=0;
	}
    n1 = -mm1
    n2 = -mm2

    ft2d_(x, y, work1, work2, &n1, &n2);
    
    
//
//  Rearrange in standard form and output
//
	int na = m2+1;
	int nb = m2-1;
	int nc = mm2-1;
	
    for (j = 0; j <mm2; ++j)
    {
		for (i= 0; i< m1; ++i)
		{
		    l = j * mm1 + i;
		    r(i,j) = x[l]*4.0/var;
		}
    }

    for (i = 0; i <m1; ++i)
    {
		for (k= 0; k< nb; ++k)
		{
			ka = k+na;
		    p(k) = r(i,ka);
		}
		for (k= 0; k< m2; ++k)
		{
			ka = k+nb;
		    p(ka) = r(i,k);
		}
		for (k= 0; k< nc; ++k)
		{
			xout(i,k) = p(k);
		}
    }



    delete[] x;
    delete[] y;
    delete[] work1;
	delete[] work2;

	return var;
}

