#include <math.h>
#include "RWFile.h"
#include "mf.h"

//	pixval: Measure
//  fileout: Archivo de Output
//
int CoherenceLength(simplmat <double> &pixval, char * outFile,
	int minBoxSize, int maxBoxSize, int deltaBoxSize)
{
	simplmat <double> box;
	simplmat <double> var;
	
	double numBoxes, globalMean,cnt, tempVar;

	int  boxSize, i, j, iRow, iCol, ix, iy, xDim, yDim, yMax,xMax;
	int  actBoxSize=0,yResto=0,xResto=0;

	xDim = pixval.getRows();
	yDim = pixval.getCols();
    
	if( maxBoxSize > yDim/2 || maxBoxSize > xDim/2  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/2;
	int numBoxSizes = (maxBoxSize/deltaBoxSize) + 1;
	box.resize(numBoxSizes);
	box(0)= minBoxSize;
	int kk;
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=box(i-1)+deltaBoxSize;
        kk = box(i);

		if(box(i) > maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
		
	}
	
	var.resize(numBoxSizes);
	var.fill(0.0);

	double totalValue = 0.0;
	double totalValue1 = 0.0;
	
	for(iy=0; iy<yDim; iy++ )
		for(ix=0; ix<xDim; ix++ )
			totalValue+=pixval(ix,iy);			

	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
//	!	WINDOW MOVEMENT
		actBoxSize = box(boxSize);
        if( actBoxSize == 0 )
        	goto endWinMovement;
            
		yResto = yDim % actBoxSize;
		xResto = xDim % actBoxSize;
		numBoxes = ((yDim-yResto)/actBoxSize) * ((xDim-xResto)/actBoxSize);

		if( !(yResto==0 && xResto==0) )
		{
			totalValue1 = 0.0;
	        yMax = yDim - yResto;
    	    xMax = xDim - xResto;
			for(iy=0; iy<yMax; iy++)
				for(ix=0; ix<xMax; ix++ )
					totalValue1+=pixval(ix,iy);

			globalMean = totalValue1 / numBoxes;
		}
		else
			globalMean = totalValue / numBoxes;

		tempVar =0;
		kk=0;

        yMax = yDim - yResto - actBoxSize;
        xMax = xDim - xResto - actBoxSize;

		for(iRow=0; iRow <= yMax; iRow+=actBoxSize )
			for(iCol=0; iCol <= xMax; iCol+=actBoxSize )
			{
				cnt = 0.0;

				// WINDOW EXAMINATION
				for(iy=iRow; iy<iRow+actBoxSize; iy++)
					for(ix=iCol; ix<iCol+actBoxSize; ix++)
						cnt+=pixval(ix,iy);

				tempVar += (globalMean-cnt)*(globalMean-cnt);
				kk++;
				
			}
		endWinMovement:
		var(boxSize) = tempVar/numBoxes;
	}

	ofstream oFile(outFile);
	oFile << "bSize\tvar\tcLen\tlog10(cLen)\tcLen1\tlog10(cLen1)\n";
	double cLen=0;
	double cLen1=0;
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		if(box(boxSize) != 0 )
		{
			cLen = var(boxSize)/(box(boxSize)*box(boxSize));
            cLen1= sqrt(var(boxSize)*box(boxSize)*box(boxSize));
			oFile << box(boxSize) << "\t" << var(boxSize) << "\t" << cLen << "\t";
			oFile << log10(cLen) << "\t" << cLen1 << "\t" << log10(cLen1) << endl;
		}
	}
	
}


