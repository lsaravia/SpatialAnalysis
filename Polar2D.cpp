// Transforms the rearranged spectrum (m/2+1)*n into polars
// and expresses it in polar segments:
//                    radius = 1,2,...,9(maxsqrt(p**2+q**2) = sqroot(64+16)=8.94).
//                    angle = 0,10,...,170 degrees.
//  Output
//		unweighted polar spectrum
//  	weighted polar spectrum
//  	number of elements and weights in each segment
//
//  NOTE: It modifies the f matrix

#include "Spectral.h"
#include <algorithm>
#include <iomanip>
#include <math.h>
#include "fortify.h"
using namespace std;

void SpectralAnalysis::Polar2D(
					int const &m,				// Rows (dimX) of the data matrix
					int const &n,				// Cols (dimY) of the data matrix 
					simplmat<double>& f,		// Rearranged spectrum
					simplmat<double>& pol		// Polar spectra (4 x Diag) R,#elems R, Theta, # Theta
					)
{
	if( n != f.getCols() )
	{
		cerr << "Error Cols wrong" << endl;
		exit(1);
	}
		
	int m1 = f.getRows();
	if( m1 != m/2+1 )
	{
		cerr << "Error Rows wrong" << endl;
		exit(1);
	}
	
	int na = n/2;
	int ma = m/2;
	int n1 = na+1;

    int ang,r,outt;

    double diag = m*m + n*n;
    
    outt = int(0.5*sqrt(diag)+1);
    outt = max(outt,20);

	simplmat <double> rad(outt);
	simplmat <double> theta(outt);
	simplmat <double> zr(outt);
	simplmat <double> zt(outt);
	simplmat <double> cr(outt);
	simplmat <double> ct(outt);
    
//
//        EXPRESS SPECTRUM AS A PERCENTAGE OF VARIANCE
//
    f(0,0) 			*=0.5;
    f(0,n1-1) 		*=0.5;
    f(m1-1,0)		*=0.5;
    f(m1-1,n1-1) 	*=0.5;

	double sum = 0.0;

	int i,j;
	for(i=1; i<ma; i++)
	{
		for(j=0; j<n; j++)
			sum += f(i,j);
	}

	for(j=0; j<n1; j++)
   		sum += f(0,j) + f(ma,j);

	double div = 100.0/sum;

	for(i=0; i<m1; i++)
	{
		for(j=0; j<n; j++)
			f(i,j) *= div;
	}
//
//       SET VECTORS TO ZERO
//

    rad.fill(0.0);
    theta.fill(0.0);
    cr.fill(0.0);
	ct.fill(0.0);
//
//       GROUP IN POLAR COORDINATES
//

	double con = 180.0/M_PI;

	int jz;
	for(j=0; j<na; j++)
	{
		jz = na - j -1;
		rad(jz) += f(0,j);
		theta(19) += f(0,j);
		cr(jz) += 1.0;
		ct(19) += 1.0;
	}

	double zi = 0.0;
	int n2 = na+1;
	int zj;
	for(i=1; i<ma; i++)
	{
		for(j=n2; j<n; j++)
		{
			zi = i;
			zj = j-na;
			ang = 1 + int((con*atan(zi/zj) + 5.0)*0.1 - 0.0001);
			r = 1 + int(sqrt(zi*zi + zj*zj) - 0.0001);
			theta(ang-1) += f(i,j);
			rad(r-1) += f(i,j);
			cr(r-1) += 1.0;
			ct(ang-1) +=1.0;
		}
	}
	
	int iz;
	for(i=1;i<m1;i++)
	{
		iz = i-1;
		rad(iz) += f(i,na);
		theta(9) += f(i,na);
		cr(iz) += 1.0;
		ct(9) += 1.0;
	}

	for(i=1;i<m1;i++)
	{
		for(j=0;j<na;j++)
		{
			zi = i;
			zj = na-j;
			ang = 1 + int((180.0 - con*atan(zi/zj) + 5.0)*0.1 - 0.0001);
			r = 1 + int(sqrt(zi*zi + zj*zj) - 0.0001);
			theta(ang-1) += f(i,j);
			rad(r-1) += f(i,j);
			cr(r-1)  += 1.0;
			ct(ang-1) += 1.0;
		}
	}


	int l;
/* 	cout << "\n\n\tRadius\tAngle" << endl;
	for(l=0;l<outt;l++)
	{
		cout << setw(2) << l << "\t";
 		cout.form("%8.5f", rad(l)) << "\t";
		cout.form("%8.5f", theta(l)) << endl;
	}
*/
//
//      EXPRESS WEIGHTS AS A % OF TOTAL WEIGHT
//
	double rtot = 0.0;
	double thetatot = 0.0;

	for(l=0;l<outt;l++)
	{
		rtot += cr(l);
		thetatot += ct(l);
	}
	double facr = 100.0/(rtot + 1.0e-30);
	double fact = 100.0/(thetatot + 1.0e-30);

	for(l=0;l<outt;l++)
	{
		zr(l) = cr(l)*facr;
		zt(l) = ct(l)*fact;
	}
	theta(0) += theta(18)+theta(19);
	
	zt(0) += zt(18)+zt(19);
	ct(0) += ct(18)+ct(19);
	

	for(l=0;l<18;l++)
		theta(l) = theta(l)/(zt(l) + 1.0e-30);

	for(l=0;l<outt;l++)
		rad(l) = rad(l)/(zr(l) + 1.0e-30);
//
//       OUTPUT WEIGHTS AND SCALED POLAR SPECTRA
//
	pol.resize(outt,4,0.0);
	for(l=0;l<18;l++)
	{
		pol(l,0)=rad(l);
		pol(l,1)=cr(l);
		pol(l,2)=theta(l);
		pol(l,3)=ct(l);
	}
	for(l=18;l<outt;l++)
	{
		pol(l,0)=rad(l);
		pol(l,1)=cr(l);
	}
/*	
 	cout << "\n\tRadius\tAngle" << endl;
	for(l=0;l<18;l++)
	{
		cout << setw(2) << l << "\t";
 		cout.form("%8.5f", rad(l)) << "\t";
		cout.form("%8.5f", theta(l)) << endl;
	}
	for(l=18;l<outt;l++)
	{
		cout << setw(2) << l << "\t";
 		cout.form("%8.5f", rad(l)) << endl;
	}

	cout << "\nElements" << "\t\t\tWeights" << endl;
 	cout << "\tRadius\tAngle\tRadius\tAngle" << endl;
	for(l=0;l<18;l++)
	{
		cout << setw(2) << l << "\t";
 		cout.form("%8.5f", cr(l)) << "\t";
 		cout.form("%8.5f", ct(l)) << "\t";
 		cout.form("%8.5f", zr(l)) << "\t";
 		cout.form("%8.5f", zt(l)) << endl;
	}
	for(l=18;l<outt;l++)
	{
		cout << setw(2) << l << "\t";
 		cout.form("%8.5f", cr(l)) << "\t        \t";
 		cout.form("%8.5f", zr(l)) << endl;
	}
*/	
}
