// 	PROGRAM MFRAC
//	BASED ON JORDI MACH'S
//	C++ PROGRAM AND MULTIFRACTAL APPROACH
#include <math.h>
#include "RWFile.h"
#include "mf.h"

// 	Q: Valores de Q posibles 
//	pixval: Measure
//  fileout: Archivo de Output
//
//  Referencia: Numerical Estimates of Generalized Dimensions D(q) For negative q
//              Romualdo Pastor-Satorras
//              Rudolf H. Riedi
//              
//  EBA: Enlarged box algorithm
//
int MultifractalEBA(simplmat <double> &pixval,simplmat <double> &q, char * outFile,
	int minBoxSize, int maxBoxSize, int numBoxSizes, char normalize)
{
//	simplmat <double> piQ;
//	simplmat <double> piHat;
	simplmat <double> box;
	simplmat <double> tauQ;
	simplmat <double> alphaQ;
	simplmat <double> fQ;
	simplmat <int> boxIni;
	
	double sumAlphaQ, sumFQ,cnt,cnt1,qT,piQT,piHatT,tauQT;

	int  boxSize, iq, i, j, iRow, iCol, ix, iy, qNum, xDim, yDim;
	int  actBoxSize=0,yResto=0,xResto=0,iRowFinal=0,iColFinal=0;

	qNum = q.getRows();
	xDim = pixval.getRows();
	yDim = pixval.getCols();

    cnt=0;
    if(toupper(normalize)=='S')
    {
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				cnt+=pixval(ix,iy);
   
		for(iy=0; iy<yDim; iy++)
			for(ix=0; ix<xDim; ix++)
				pixval(ix,iy)/=cnt;
	}
	
	if( maxBoxSize > yDim/3 || maxBoxSize > xDim/3  )
		maxBoxSize = (xDim<yDim ? xDim : yDim)/3;
/*	int numBoxSizes = maxBoxSize/deltaBoxSize;
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
*/

	double deltaBoxSize = (log10(static_cast<double>(maxBoxSize)) - log10(static_cast<double>(minBoxSize)))/ numBoxSizes;
	box.resize(numBoxSizes);

	box(0)= minBoxSize;
	for(i=1; i<numBoxSizes; i++)
	{
		box(i)=pow10(log10(box(i-1))+deltaBoxSize);
		if(box(i) > maxBoxSize )
		{
			numBoxSizes = i;
			break;
		}
		
	}
	for(i=0; i<numBoxSizes; i++)
		box(i) = ceil(box(i));

	for(i=numBoxSizes-1; i>1; i--)
	{
		if( box(i) == box(i-1) )
		{
			for(int ii=i; ii<numBoxSizes-1; ii++)
				box(ii) = box(ii+1);
			numBoxSizes--;
		}
	}
		

	double numRep=0;
		
	alphaQ.resize(numBoxSizes,qNum,0.0);
	tauQ.resize(numBoxSizes,qNum,0.0);
	fQ.resize(numBoxSizes,qNum,0.0);
	boxIni.resize(4,2,0);
	
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		for(iq=0;iq<qNum; iq++)
		{
//			
//			TauQ
//
			qT = q(iq);

//		!	WINDOW MOVEMENT
			actBoxSize = box(boxSize);
			yResto = yDim % actBoxSize;
			xResto = xDim % actBoxSize;
			
			boxIni.fill(actBoxSize);
 			numRep=1;
			boxIni(1,0)= actBoxSize;
			boxIni(3,0)= actBoxSize;
			boxIni(2,1)= actBoxSize;
			boxIni(3,1)= actBoxSize;

            iRowFinal=1+yDim-actBoxSize*2;
            iColFinal=1+xDim-actBoxSize*2;

/*			if( xResto>0 && yResto>0  )
			{
				numRep=4;
				boxIni(1,0)= actBoxSize+xResto;
				boxIni(3,0)= actBoxSize+xResto;
				boxIni(2,1)= actBoxSize+yResto;
				boxIni(3,1)= actBoxSize+yResto;

			}
 			if( xResto>0  && yResto==0)
			{
				numRep=2;
				boxIni(1,0)= actBoxSize+xResto;
			}

			if( yResto>0  && xResto==0)
			{
				numRep=2;
				boxIni(1,1)= actBoxSize+yResto;
			}
*/
			for(int rep=0;rep<numRep;rep++)
			{

				piQT=0;
				tauQT=0;

				for(iRow=boxIni(rep,1); iRow < iRowFinal; iRow+=actBoxSize )
					for(iCol=boxIni(rep,0); iCol < iColFinal; iCol+=actBoxSize )
					{
						cnt = 0.0;
						cnt1= 0.0;

	//				!	WINDOW EXAMINATION
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
							{
								cnt1+=pixval(ix-actBoxSize,iy-actBoxSize);
								cnt1+=pixval(ix-actBoxSize,iy);
								cnt1+=pixval(ix-actBoxSize,iy+actBoxSize);
								cnt1+=pixval(ix           ,iy-actBoxSize);
								cnt1+=pixval(ix           ,iy+actBoxSize);
								cnt1+=pixval(ix+actBoxSize,iy-actBoxSize);
								cnt1+=pixval(ix+actBoxSize,iy);
								cnt1+=pixval(ix+actBoxSize,iy+actBoxSize);
								cnt+=pixval(ix,iy);
							}
						cnt1 += cnt;

//						if( !(cnt==0.0 && qT <= 0 ) )
						if( !(cnt==0.0) ) 
							piQT+=pow(cnt1,qT);
						
					}
				
				if( piQT > 0.0 )
					tauQT=log10(piQT);
			
				// TO DO AlphaQ AND FQ

				sumAlphaQ=0;
				sumFQ=0;
				piHatT=0;

				//	WINDOW MOVEMENT
				for(iRow=boxIni(rep,1); iRow < iRowFinal; iRow +=actBoxSize )
					for(iCol=boxIni(rep,0); iCol < iRowFinal; iCol +=actBoxSize )
					{

						cnt=0.0;
						cnt1= 0.0;

	//				!	WINDOW EXAMINATION
						for(iy=iRow; iy<iRow+actBoxSize; iy++)
							for(ix=iCol; ix<iCol+actBoxSize; ix++)
							{
								cnt1+=pixval(ix-actBoxSize,iy-actBoxSize);
								cnt1+=pixval(ix-actBoxSize,iy);
								cnt1+=pixval(ix-actBoxSize,iy+actBoxSize);
								cnt1+=pixval(ix           ,iy-actBoxSize);
								cnt1+=pixval(ix           ,iy+actBoxSize);
								cnt1+=pixval(ix+actBoxSize,iy-actBoxSize);
								cnt1+=pixval(ix+actBoxSize,iy);
								cnt1+=pixval(ix+actBoxSize,iy+actBoxSize);
								cnt+=pixval(ix,iy);
							}
						cnt1 += cnt;
						
//						if( (piQT != 0.0) && ( (cnt!=0.0) || (qT > 0.0) ) )
						if( (piQT != 0.0) && (cnt!=0.0) )								
							piHatT=pow(cnt1,qT)/piQT;						
						else
							piHatT=0;

						if( cnt > 0.0 )
							sumAlphaQ+=piHatT*log10(cnt1);
						

						if(piHatT > 0.0)							
							sumFQ+=piHatT*log10(piHatT);
					}

				alphaQ(boxSize,iq)+=sumAlphaQ;
				fQ(boxSize,iq)+=sumFQ;
				tauQ(boxSize,iq)+=tauQT;
//				cout << "q: " << qT << "Box: " << actBoxSize << "TauQ: " << tauQT << "\t" << tauQ(boxSize,iq) << endl;
			}
			tauQ(boxSize,iq)/=static_cast<double>(numRep);
			alphaQ(boxSize,iq)/=static_cast<double>(numRep);
			fQ(boxSize,iq)/=static_cast<double>(numRep);
//			cout << "q: " << qT << "Box: " << actBoxSize << "TauQ: " << tauQ(boxSize,iq) << endl;
			
		}
	}

    string::size_type pos=0;
	string outF = outFile;
	if( (pos=outF.find("\\"))!= string::npos )
		outF = outF.substr(pos);
		
	string tFileName("t.");
	tFileName+=outF;
	ofstream tFile(tFileName.c_str());
   
	string aFileName("a.");
	aFileName+= outF;
	ofstream aFile(aFileName.c_str());
   
	string fFileName("f.");
	fFileName+= outF;
	ofstream fFile(fFileName.c_str());

	tFile << "Box Size" << "\t" << "Log Box";
	aFile << "Box Size" << "\t" << "Log Box";
	fFile << "Box Size" << "\t" << "Log Box";
	for(i=0;i<qNum;i++)
	{
		tFile << "\t" << q(i);
		aFile << "\t" << q(i);
		fFile << "\t" << q(i);
	}
	tFile << endl;
	aFile << endl;
	fFile << endl;
	for(boxSize=0; boxSize<numBoxSizes; boxSize++)
	{
		tFile << box(boxSize) << "\t" << log10(box(boxSize));
		for(i=0;i<qNum;i++)
			tFile << "\t" << tauQ(boxSize,i);
		tFile << endl;
		
		aFile << box(boxSize) << "\t" << log(box(boxSize));
		for(i=0;i<qNum;i++)
			aFile << "\t" << alphaQ(boxSize,i);
		aFile << endl;
		
		fFile << box(boxSize) << "\t" << log(box(boxSize));
		for(i=0;i<qNum;i++)
			fFile << "\t" << fQ(boxSize,i);
		fFile << endl;
	}
	
	// Regression of the log - log variables
	//
	//
	string sFileName("s.");
	sFileName+=outF;
	ofstream sFile(sFileName.c_str());
	sFile << "q\tTau\talfa\tf(alfa)\tR-Tau\tR-alfa\tR-f\tSD-Tau\tSD-alfa\tSD-f" << endl;
	for(i=0; i<qNum; i++)
	{
		double sumx=0,
			sumyt = 0,
			sumya = 0,
			sumyf = 0,
			sumxyt = 0,
			sumxya = 0,
			sumxyf = 0,
			sumxsq = 0,
			sumysqt = 0,
			sumysqa = 0,
			sumysqf = 0,
			sdbt=0,
			sdba=0,
			sdbf=0;
		
		for(boxSize=0; boxSize<numBoxSizes; boxSize++)
		{
			double x = log10(box(boxSize));
			double yt = tauQ(boxSize,i);
			double ya = alphaQ(boxSize,i);
			double yf = fQ(boxSize,i);
			
			sumx += x;
			sumyt += yt;
			sumya += ya;
			sumyf += yf;
			sumxyt += x * yt;
			sumxya += x * ya;
			sumxyf += x * yf;
			sumxsq += x * x;
			sumysqt += yt * yt;
			sumysqa += ya * ya;
			sumysqf += yf * yf;
		}
		double xbar  = sumx/numBoxSizes ;
		double ybart = sumyt/numBoxSizes ;
		double ybara = sumya/numBoxSizes ;
		double ybarf = sumyf/numBoxSizes ;
	
		double bt = (sumxyt - numBoxSizes *xbar*ybart)/(sumxsq - numBoxSizes *xbar*xbar);
		double ba = (sumxya - numBoxSizes *xbar*ybara)/(sumxsq - numBoxSizes *xbar*xbar);
		double bf = (sumxyf - numBoxSizes *xbar*ybarf)/(sumxsq - numBoxSizes *xbar*xbar);
	//	a = ybar - (b*xbar);
		double rt = (sumxyt - numBoxSizes *xbar*ybart)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqt - numBoxSizes *ybart*ybart));
		double ra = (sumxya - numBoxSizes *xbar*ybara)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqa - numBoxSizes *ybara*ybara));
		double rf = (sumxyf - numBoxSizes *xbar*ybarf)/sqrt((sumxsq - numBoxSizes  * xbar * xbar) * (sumysqf - numBoxSizes *ybarf*ybarf));
		if( numBoxSizes  > 4 )
		{
			sdbt = sqrt((1-rt*rt)*(sumysqt-numBoxSizes*ybart*ybart))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdba = sqrt((1-ra*ra)*(sumysqa-numBoxSizes*ybara*ybara))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
			sdbf = sqrt((1-rf*rf)*(sumysqf-numBoxSizes*ybarf*ybarf))/((numBoxSizes -4)*(sumxsq-numBoxSizes *xbar*xbar));
		}
		sFile << q(i) << "\t" << bt << "\t" << ba << "\t" << bf << "\t"
				<< rt*rt << "\t" << ra*ra << "\t" << rf*rf << "\t"
				<< sdbt  << "\t" << sdba << "\t" << sdbf << endl;
	}
}


