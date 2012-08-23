// THIS PROGRAM EVALUATES THE PERIODOGRAM DIRECTLY FROM THE FOURIER
// COEFFICIENTS AND OUTPUTS VALUES OF I(J,K) FOR J=0,...,M/2 AND K=0,...,N-1.
// CHANNELS USED ARE AS FOLLOWS:

//   FT01 - INPUT MATRIX SPECIFICATION
//   FT02 - INPUT DATA MATRIX
//   FT03 - OUTPUT PERIODOGRAM
//   FT04 - OUTPUT MATRIX SIZE AND VARIANCE (I.E. SS/M*N)

#include "Spectral.h"
#include "fortify.h"
#include <math.h>

using namespace std;

double SpectralAnalysis::Period2D(simplmat<double>& x,simplmat<double>& xout)
{
	simplmat <double> a,b,p;
	double total,con1,con2,pie,ang,inv,ss,isize;

//	readSeed(filein,x);
	
	int n = x.getCols();  // ReadSeed lee x,y 
	int m = x.getRows();
	if( !((m%2)==0 && (n%2)==0) )
	{
		cerr << "Dimension must be multiplo of 2" << endl;
		exit(1);
	}
		
//	Cols = n;
//	Rows = m;
	
	int n1 = n/2+1;
	int m1 = m/2+1;
	a.resize(n);
	b.resize(n);
	p.resize(n);
	xout.resize(m1,n);
//
//          SET PARAMETERS
//
	total = 0.0;
	ss = 0.0;
	isize = m*n;
	inv = 1.0/double(isize);
//	pie = 3.141592653589793;
	con1 = (2.0*M_PI)/double(m);
	con2 = (2.0*M_PI)/double(n);
//
//          CORRECT FOR MEAN AND EVALUATE SS
//
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			total += x(i,j);
			
	total *= inv;
	
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			x(i,j) -= total;
			ss += x(i,j)*x(i,j);
		}
	ss *= inv;
			
//
//          EVALUATE FOURIER COEFFICIENTS AND PERIODOGRAM
//
	int k,l;
	for(i=0;i<m1;i++)
	{
		for(j=0;j<n;j++)
		{
			a(j)=0.0;
			b(j)=0.0;
			for(k=0;k<m;k++)
				for(l=0;l<n;l++)
				{
					ang = con1*(i*(k+1)) + con2*(j*(l+1));
					a(j) = a(j) + x(k,l)*cos(ang);
					b(j) = b(j) + x(k,l)*sin(ang);
				}
			p(j) = ( a(j)*a(j) + b(j)*b(j) )*inv;
		}

		for(j=0;j<n;j++)
			xout(i,j)=p(j);
	}
	return ss;
}


void SpectralAnalysis::Transpose( simplmat<double> & data )
{
	simplmat<double> temp(data);
	data.resize(temp.getCols(),temp.getRows());
	int cols= data.getCols();
	int rows= data.getRows();
	int i,j;
	for( i=0; i<rows; i++)
		for( j=0; j<cols; j++)
			data(i,j) = temp(j,i);
	
}

