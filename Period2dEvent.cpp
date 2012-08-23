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

void SpectralAnalysis::Period2DEvent(simplmat<double>& x,double longX,double longY,int range,simplmat<double>& xout)
{
	simplmat <double> a,b;
	double con1,ang;

//	readSeed(filein,x);
	
	int n = x.getCols();  // ReadSeed lee x,y 
	int m = x.getRows();
	if( n!=2 )
	{
		cerr << "Dimension must be multiplo of 2" << endl;
		exit(1);
	}
		
//	Cols = n;
//	Rows = m;
	
	int q = range; //q=n1
	int p = range/2+1;   //p=m1
	a.resize(q); 
	b.resize(q);
	xout.resize(p,q);
//
//          SET PARAMETERS
//
	con1 = 2.0*M_PI;
//
//          CORRECT FOR MEAN AND EVALUATE SS
//
	int i,j;
	
	for(i=0;i<m;i++)
	{
		x(i,0) /= longX;
		x(i,1) /= longY;
	}
			
//
//          EVALUATE FOURIER COEFFICIENTS AND PERIODOGRAM
//
	int k,l;
	for(i=0;i<p;i++)
	{
		for(j=0;j<q;j++)
		{
			a(j)=0.0;
			b(j)=0.0;
			for(k=0;k<m;k++)
				{
					ang = con1*(i*x(k,0) + j*x(k,1));
					a(j) += cos(ang);
					b(j) += sin(ang);
				}
		}

		for(j=0;j<q;j++)
		{
//			xout(i,j)=( a(j)*a(j) + b(j)*b(j) );
			if(j==0 && i==0 )
				xout(i,j)=( a(j)*a(j) + b(j)*b(j)-m  )/m;
			else
				xout(i,j)=( a(j)*a(j) + b(j)*b(j) )/m;
		}
	}
}



