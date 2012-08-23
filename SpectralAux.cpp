#include <math.h>
#include "Spectral.h"
#include "smattpl.h"

// Copy from origData to data and resize data to a power of 2 side lengths
//
void SpectralAnalysis::SetPow2(simplmat<double> & data,simplmat<double> & origData, int initialPos)
{
	int oRows=origData.getRows();
	int oCols=origData.getCols();

	double e2rows = int(log2(oRows));
	double p2rows = pow2(e2rows);
	if( p2rows > oRows )
	{
		e2rows -=1;
		p2rows = pow2(e2rows);
	}

	double e2cols = int(log2(oCols));
	double p2cols = pow2(e2cols);
	if( p2cols > oCols )
	{
		e2cols -=1;
		p2cols = pow2(e2cols);
	}
	
	data.resize(p2rows,p2cols,0.0);
    int rows = static_cast<int>(p2rows);
    int cols = static_cast<int>(p2cols);
    int rIni = 0;
    int rEnd = 0;
    int cIni = 0;
    int cEnd = 0;

	switch(initialPos)
    {
    case 0:
    	rIni = 0;
        rEnd = rows;
        cIni = 0;
        cEnd = cols;
        break;
    case 1:
    	rIni = 0;
        rEnd = rows;
        cIni = oCols-cols;
        cEnd = oCols;
        break;

    case 2:
    	rIni = oRows-rows;
        rEnd = oRows;
        cIni = 0;
        cEnd = cols;
        break;
    
    case 3:
    	rIni = oRows-rows;
        rEnd = oRows;
        cIni = oCols-cols;
        cEnd = oCols;
        break;
    }
    
	int r,c;		
	for(r=rIni;r<rEnd;r++)
		for(c=cIni;c<cEnd;c++)
			data(r-rIni,c-cIni) = origData(r,c);

}

// Mantiene los "nExtr" valores minimos de "pspect"
//
void SpectralAnalysis::spectMin(double const &pspect, simplmat<double> & pMin, int & l, int nExtr )
{
	if( pspect < pMin(l,0) )
	{
		for(int k=nExtr-1; k>=0; k--)
		{
			if( pspect < pMin(l,k) )
			{
				for(int kk=1; kk<=k; kk++)
					pMin(l,kk-1) = pMin(l,kk);
				pMin(l,k) = pspect;
				break;
			}
		}
	}
}

void SpectralAnalysis::spectMax(double const &pspect, simplmat<double> & pMax, int & l, int nExtr )
{
	if( pspect > pMax(l,0) )
	{
		for(int k=nExtr-1; k>=0; k--)
		{
			if( pspect > pMax(l,k) )
			{
				for(int kk=1; kk<=k; kk++)
					pMax(l,kk-1) = pMax(l,kk);
				pMax(l,k) = pspect;
				break;
			}

		}
	}
}


