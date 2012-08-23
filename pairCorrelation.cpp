#include <math.h>
#include <algorithm>
#include "RWFile.h"
#include "mf.h"

struct data
{
	double x;
	double y;
	double m;
};

//	pixval: Measure
//  fileout: Archivo de Output
//
int PairCorrelation(simplmat <double> &pixval, char * outFile,
	double minBoxSize, double maxBoxSize, double deltaBoxSize, int maxPoints)
{
	simplmat <double> box;
	simplmat <double> g;
	simplmat <double> kmm;
	simplmat <double> Kmm;
	simplmat <double> K;
	simplmat <data> d;
	
	
	double numBoxes;

	int  boxSize, i, j, ix, iy, jx, jy, xDim, yDim;
	int  actBoxSize=0;

	xDim = pixval.getRows();
	yDim = pixval.getCols();
    
	if( maxBoxSize > yDim/2 || maxBoxSize > xDim/2  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/2;
	int numBoxSizes = (maxBoxSize/deltaBoxSize) + 1;
	box.resize(numBoxSizes);
	box(0)= minBoxSize;

	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=box(i-1)+deltaBoxSize;

		if(box(i) > maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
		
	}
	
	g.resize(numBoxSizes);
	kmm.resize(numBoxSizes);
	K.resize(numBoxSizes);
	Kmm.resize(numBoxSizes);

	g.fill(0.0);
	kmm.fill(0.0);
	K.fill(0.0);
	Kmm.fill(0.0);

	double totalValue = 0.0;
	long totalPoints = 0;
	
	for(iy=0; iy<yDim; iy++ )
		for(ix=0; ix<xDim; ix++ )
		{
			if( pixval(ix,iy) > 0 )
			{
				totalPoints ++;
				totalValue+=pixval(ix,iy);			
			}
		}

	double markMean = totalValue / (xDim*yDim);
	double pointDensity = static_cast<double>(totalPoints) / (xDim*yDim);
	double kernelWidth = 0.2 / sqrt(pointDensity);
	
	i=0,
	d.resize(totalPoints);
	for(iy=0; iy<yDim; iy++ )
		for(ix=0; ix<xDim; ix++ )
		{
			if( pixval(ix,iy) > 0 )
			{
				d(i).x = ix;
				d(i).y = iy;
				d(i).m = pixval(ix,iy);			
				i++;
			}
		}

    double origMarkMean = markMean;
    double origPointDensity = pointDensity;
    long   origTotalPoints  = totalPoints;

	if(maxPoints>0 && maxPoints<totalPoints)
	{
		data * ptr= d.pointer();
	    random_shuffle(ptr, ptr + totalPoints);
		totalPoints = maxPoints;
		pointDensity = static_cast<double>(totalPoints) / (xDim*yDim);
        totalValue=0;
		for(i=0; i<totalPoints; i++)
				totalValue+=d(i).m;
        markMean = totalValue / (xDim*yDim);
        kernelWidth = 0.2 / sqrt(pointDensity);
	}
		
	double norm=0,x=0,w=0,s=0,r=0;


	for(i=0; i<totalPoints; i++)
		for(j=0; j<totalPoints; j++)
            if( j!=i )
            {
				for(boxSize=0; boxSize<numBoxSizes; boxSize++)
				{
					r = box(boxSize);
		    	    if( r >0 )
					{
                        ix = d(i).x;
                        jx = d(j).x;
                        iy = d(i).y;
                        jy = d(j).y;
						norm = sqrt( (ix-jx)*(ix-jx) + (iy-jy)*(iy-jy) );
						x = norm - r;
						if( fabs(x) < kernelWidth )
						{
							w = 3/(4*kernelWidth)*(1-(x*x) /(kernelWidth*kernelWidth));
							s = xDim*yDim - r*(2*xDim+2*yDim - r) / M_PI;
							g(boxSize) += w /(pointDensity*pointDensity*2*M_PI*r*s);

						}
						if( norm <= r )
						{
							s = xDim*yDim - norm*(2*xDim+2*yDim - norm) / M_PI;
							K(boxSize) += 1/( pointDensity*pointDensity*s );
//							Kmm(boxSize) += d(i).m*d(j).m/( pointDensity*pointDensity*s );
							Kmm(boxSize) += d(i).m*d(j).m/s;
						}
					}
				}
            }
	for(i=0; i<totalPoints; i++)
		for(j=0; j<totalPoints; j++)
            if( j!=i )
            {
				for(boxSize=0; boxSize<numBoxSizes; boxSize++)
				{
					r = box(boxSize);
		    	    if( r >0 )
					{
						norm = sqrt( (d(i).x-d(j).x)*(d(i).x-d(j).x) + (d(i).y-d(j).y)*(d(i).y-d(j).y) );
						x = norm - r;
						if( fabs(x) < kernelWidth )
						{
							w = 3/(4*kernelWidth)*(1-(x*x) /(kernelWidth*kernelWidth));
							s = xDim*yDim - r*(2*xDim+2*yDim - r) / M_PI;
//							kmm(boxSize) +=w*d(i).m*d(j).m /(pointDensity*pointDensity*markMean*markMean*2*M_PI*r*s*g(boxSize));
							kmm(boxSize) +=w*d(i).m*d(j).m /(markMean*markMean*2*M_PI*r*s*g(boxSize));							
						}
					}
				}
            }

    double L,Lmm;
	ofstream oFile(outFile);
	oFile << "dim:\t" << xDim <<  "\t" << yDim << endl;
	oFile << "Point Density:\t" << origPointDensity << endl;
	oFile << "Mean Mark    :\t" << origMarkMean << endl;
    oFile << "Points sampled\t" << totalPoints << "\t";
    oFile << "% sampled\t" << static_cast<double>(totalPoints)/origTotalPoints*100 << endl;
	oFile << "r\tg(r)\tkmm(r)\tK(r)\tKmm(r)\tL(r)\tLmm(r)\n";
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		if(box(boxSize) != 0 )
		{
			oFile << box(boxSize) << "\t" << g(boxSize) << "\t" << kmm(boxSize) << "\t";
			oFile << K(boxSize) << "\t" << Kmm(boxSize) << "\t";
            L=sqrt(K(boxSize)/M_PI);
            Lmm=sqrt(Kmm(boxSize)/(M_PI*markMean*markMean));
            oFile << L << "\t" << Lmm << endl;
            
            
		}
	}
	
}


