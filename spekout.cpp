//	This program inputs the top half of the complete spectrum
//	(i=1,...,m/2+1 ; j=1,...,n) ; interchanges the top left and
//	top right quadrants ; expresses the spectrum as a proportion
//	of total variance ; and , includes only those elements which
//	contribute to at least a% of total variance.
//  returns the explained variace fvariance
//
//  Rearranged periodogram I(j,k) with j=0,...,m/2 and k=-n/2,...,n/2-1;
#include "Spectral.h"
//#include "fortify.h"

using namespace std;

double SpectralAnalysis::Spekout(
					int const &m,				// Rows of the original data
					int const &n,				// Cols of the original data
					simplmat<double>& f,		// Complete spectrum
					simplmat<double>& g,  		// Rearranged Spectrum
					double var,				// Variance of the original data
					double a					// Percentage sobre el cual se incluye
												// 			el valor del spectrum
					)
{
	double tvar=0.0;
	
	if( n != f.getCols() )
	{
		cerr << "Error Cols wrong (n)" << endl;
		exit(1);
	}
		
	int m1 = f.getRows();
	if( m1 != m/2+1 )
	{
		if( m1 == m )
			m1 = m/2+1;
		else
		{
			cerr << "Error Rows wrong (m)" << endl;
			exit(1);
		}
	}
	
	int na = n/2;
	int ma = m/2;
	int n1 = na+1;

	double inv = 1/static_cast<double>(m*n);
	g.resize(m1,n);
//	simplmat <int> l(m1,n);

//
//         REARRANGE THE SPECTRUM
//
	int i,j,ja,jb;
	
	for(i=0; i<m1; i++)
	{
		for(j=0; j<na; j++)
		{
	      	ja = j+na;
			g(i,j) = f(i,ja);
		}
		
		for(j=na;j<n; j++)
		{
			jb = j-na;
	 		g(i,j) = f(i,jb);
		}
	}
// The rearranged spectra is in g(0..m1,0..n)

 
//
//         EXPRESS AS A PERCENTAGE OF VARIANCE
//
	if( var>0 )
	{
		double c = 200.0*inv/var;
		for(i=0;i<m1;i++)
		{
			for(j=0; j<n; j++)
	   			g(i,j) = g(i,j)*c;
		}
	    g(0,0) 			*=0.5;
	    g(0,n1-1) 		*=0.5;
	    g(m1-1,0) 		*=0.5;
	    g(m1-1,n1-1) 	*=0.5;
	}
	
// The percentagge spectra is in h(0..m1,0..n)

//
//         ELEMENTS CONTRIBUTING TO AT LEAST A% OF TOTAL VAR
//
	if( var>0 && a>0 )
	{
		for(j=0; j<n1; j++)
		{
	      	if(g(0,j)>=a)
			{
	      		tvar += g(1,j);
	      		g(0,j) = static_cast<int>(g(0,j)*10.0 );
			}
			else
				g(0,j) = 0.0;
		}

		for(i=1;i<ma;i++)
		{
			for(j=0;j<n;j++)
			{
		      	if( g(i,j)>=a)
	      		{
	      			tvar += g(i,j);
					g(i,j) = static_cast<int>( g(i,j)*10.0 );
	      		}
				else
					g(i,j) = 0.0;

			}
		}

		for(j=0;j<n1;j++)
		{
	      	if(g(ma,j)>=a)
	      	{
	      		tvar += g(ma,j);
				g(ma,j) = static_cast<int>( g(ma,j)*10.0 );
	      	}
			else
				g(ma,j) = 0.0;

		}
		
		for(j=n1;j<n;j++)
		{
			g(0,j) = 0.0;
			g(ma,j) = 0.0;
		}
		
	}
	return tvar;
}
