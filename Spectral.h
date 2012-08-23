#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "smattpl.h"
//#include "r250.h"
#include <time.h>
#include "randlib.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

class SpectralAnalysis
{
//	int Rows;  // (dimX) Esta alrevez
//	int Cols;  // (dimY)

	public:

    SpectralAnalysis(int rSeed=0) {
    					if(rSeed==0)
                        	rSeed=time(0);
                       	setall(rSeed,rSeed+1);
                        };

// OJO AL CAMBIO DE GENERADOR DE NROS AL AZAR!!!!!!!!!!
    unsigned Rand(unsigned r) { return ignuin(0,r-1); };

	void Transpose( simplmat<double> & data );
	
	// This program evaluates the periodogram directly from the fourier
	// coefficients and outputs values of i(j,k) for j=0,...,m/2 and k=0,...,n-1.
	// Retorna la varianza de los datos
	//
	double Period2D(simplmat<double>& x,simplmat<double>& xout);

	void Period2DEvent(simplmat<double>& x,double longX,double longY,int range,simplmat<double>& xout);


	double Spekout(int const &m,
				int const &n,
				simplmat<double>& f,
				simplmat<double>& g,
				double var=0,
				double a=0 );


	void Polar2D(
					int const &m,				// Rows of the data matrix
					int const &n,				// Cols of the data matrix
					simplmat<double>& f,			// Rearranged spectrum
					simplmat<double>& pol			// Rearranged spectrum
					);

	void Polar2DEvent(
					unsigned const &events, // Number of events
					simplmat<double>& f,		// Rearranged spectrum
					simplmat<double>& pol		// Polar spectra (4 x Diag) R,#elems R, Theta, # Theta
					);
		


	void Polar2DDir(
					int const &m,				// Rows of the data matrix
					int const &n,				// Cols of the data matrix
					simplmat<double>& f,			// Rearranged spectrum
					simplmat<double>& pol,			// Polar spectra
					simplmat<double>& sPs,			// Polar spectra in bins
					simplmat<double>& csPs			// Counts of Polar spectra in bins
					);

	double Spec2D(simplmat<double>& r,simplmat<double>& xout);


	void spectMin(double const &pspect, simplmat<double> & pMin, int & l, int nExtr );
	void spectMax(double const &pspect, simplmat<double> & pMax, int & l, int nExtr );
	void SetPow2(simplmat<double> & data,simplmat<double> & origData, int initialPos=0);

};


#endif  // SPECTRAL_H


