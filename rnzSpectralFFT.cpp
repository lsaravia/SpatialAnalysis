#include "RWFile.h"
#include "Randomizations.h"
#include <iomanip>
#include <math.h>
//#include "fortify.h"
#include "cdflib.h"

using namespace std;

int main(int argc, char * argv[])
{
	simplmat <double> data;
	simplmat <double> dout;
	SpectralAnalysis sa;
	
	// Do 199 Simulations and takes the 5 lowest and 5 highest values
	// to make the confidence envelopes. If the calculated value is
	// between these highest or lowest values it is significative at 5% level
	//
	int numSimul=199,numExtreme=5,windowPos=0,numExtremeT=5;
	char randz;
	double probConf=0.05;
	double probConfT=0.05;
	double probOrig=0.05;
	int bonfCorr=0,signif=0;
	// OJO Falta agregar en parametros
	int reaPer=0;
	string fType="BI";
	
	if( argc <3 )
    {
		cerr << "Usage: rnzSpectral inputFile.sed outFile fileType winPos{0,1,2,3} [R/A] [prob] [bonfCorr{0,1}] [numSimulations]" << endl;
        exit(1);
	}
	if( argc >= 4)
	{
		fType = argv[3];
    	windowPos  = atoi(argv[4]);
		randz = toupper(argv[5][0]); // R=Randomizations, A=Asintotic chisqr.
        if( windowPos > 3 || windowPos <0 )
        {
			cerr << "Out of range: 3 > windowPos > 0" << endl;
        	exit(1);
        }
	}
	
	if( randz=='R' )
	{
		if( argc == 9)
		{
			probOrig = atof(argv[6]);
			bonfCorr = atoi(argv[7]);
			numSimul   = atoi(argv[8]);
		}
		else
		{
			cerr << "Invalid number of arguments" << endl;
			exit(1);
		}
	}
	else
	{
		if( argc == 8 )
		{
			probOrig = atof(argv[6]);
			bonfCorr = atoi(argv[7]);
		}
		else
		{
			cerr << "Invalid number of arguments" << endl;
			exit(1);
		}
	}
	
	
	RWFile file;
	string fName = argv[1];
	string outFName = argv[2];
	if( fName.find(".sed")!=string::npos )
	{
		if(!file.ReadSeed(fName.c_str(), data, fType.c_str() ))
			exit(1);
	}
	else if( fName.find(".rst")!=string::npos )
	{
		if(!file.ReadIdrisi(fName.c_str(), data))
			exit(1);
	}
	else if( fName.find(".tif")!=string::npos )
	{
		if(!file.ReadTiff(fName.c_str(), data))
			exit(1);
	}
	else
	{
		cerr << "File type not supported: " << fName << endl;
		exit(1);
	}
		
   	sa.Transpose(data); // Los datos leidos estan en formato (x,y) y las funciones
	                    // tienen (y,x) o sea (row,col)
	int rows=data.getRows();
	int cols=data.getCols();
	double var;
	int s,l;

	simplmat <double> origData(data);	
	sa.SetPow2(data,origData,windowPos);
	int origRows = rows;
	int origCols = cols;
	rows=data.getRows();
	cols=data.getCols();

	var = sa.Spec2D(data,dout); // Calculates the periodogram usando FFT

	simplmat<double> rper;		// Rearranged periodogram
	simplmat<double> polper;	// Polar
	sa.Spekout(rows,cols,dout,rper); // Rearranges the periodogram


   	ofstream fOut;
	if(reaPer)
	{
		string rperOut="r."+outFName;
    	fOut.open(rperOut.c_str());
    	fOut <<	rper.getCols() << "\t" << rper.getRows() << endl;
    	fOut << "Rearranged periodogram" << endl;
    	for(l=0;l<rper.getRows();l++)
    	{
    		for(s=0;s<rper.getCols();s++)
//    			fOut.form("%10.6f",rper(l,s)) <<  "\t";
				fOut <<  setw(10) << fixed << showpoint << setprecision(6) << rper(l,s) << "\t";
    		fOut << endl;
    	}
    	fOut.close();
	}
	
	sa.Polar2D(rows,cols,rper,polper); // Calculates the Polar Spectrum
									   // OJO modifica rper!!!!!!!!!!!!
	int d =	int(0.5*sqrt(rows*rows + cols*cols)+1);

/*
	// R-Spectra in intervals of 10 degrees
	//
	simplmat<double> polSp;		// Polar in Bins
	simplmat<double> cpolSp;	// Counts of Polar in Bins
	sa.Polar2DDir(rows,cols,rper,polper,polSp,cpolSp);

	rperOut="b."+outFName;
	fOut.open(rperOut.c_str());
	fOut <<	18 << "\t" << d << endl;
	fOut << "Polar spectrum in Bins" << endl;
	for(l=0;l<d;l++)
	{
		for(s=0;s<18;s++)
			fOut.form("%10.6f",polSp(l,s)) <<  "\t";
		fOut << endl;
	}
	fOut << endl;
	
	fOut << "Counts Polar spectrum in Bins" << endl;
	for(l=0;l<d;l++)
	{
		for(s=0;s<18;s++)
			fOut.form("%10.6f",cpolSp(l,s)) <<  "\t";
		fOut << endl;
	}
	fOut.close();

*/

	// Outputs the rearranged periodogram
	//
	//string rperOut("rper.");
	//rperOut+=outFName;
	//if(!file.WriteSeed(rperOut.c_str(), rper))
	//	exit(1);

	fOut.open(outFName.c_str());

	if( bonfCorr )
	{
		probConf=probOrig/d;
		probConfT=probOrig/18;
	}
	else
	{
		probConf=probOrig;
		probConfT=probOrig;
	}

	if(	randz=='R' )
	{
		numExtreme = 1+probConf*numSimul/2;
        if( numExtreme <=1 )
        {
			cerr << "numExtreme <= 1" << endl;
        	exit(1);
        }
        numExtremeT = 1+probConfT*numSimul/2;
        if( numExtremeT <=1 )
        {
			cerr << "numExtremeT <= 1" << endl;
        	exit(1);
        }

		Randomizations rz;
		simplmat<double> thetamin;
		simplmat<double> thetamax;
		simplmat<double> rmin;
		simplmat<double> rmax;
		simplmat<double> rpol;		// Polar Spectrum for randomizations

		thetamin.resize(18,numExtremeT,1000.0);
		thetamax.resize(18,numExtremeT,0.0);
		rmin.resize(d,numExtreme,1000.0);
		rmax.resize(d,numExtreme,0.0);
	
		for(s=0; s<numSimul; s++)
		{
			rz.Randomize(origData);
			sa.SetPow2(data,origData);
			sa.Spec2D(data,dout);
			sa.Spekout(rows,cols,dout,rper);
			sa.Polar2D(rows,cols,rper,rpol);
			for(l=0;l<d;l++)
			{
				sa.spectMax(rpol(l,0),rmax,l,numExtreme);
				sa.spectMin(rpol(l,0),rmin,l,numExtreme);
	
	//			cout << setw(2) << l+1 << "\t";
	// 			cout.form("%8.5f", rpol(l,0)) << endl;
			}
			for(l=0;l<18;l++)
			{
				sa.spectMax(rpol(l,2),thetamax,l,numExtremeT);
				sa.spectMin(rpol(l,2),thetamin,l,numExtremeT);
	//			cout << "\t\t" << setw(4) << l*10 << "\t";
	//	 		cout.form("%8.5f", rpol(l,2)) << endl;
			}
		}

		fOut << "R Spectrum "  << "\t" << outFName << endl;
		fOut << "Randomizations: "<< "\t" << numSimul << "\t" << numExtreme << "\t" << probConf << endl;
		fOut << "Data dim:" << "\t" << origRows << "\t" << origCols << endl;
		fOut << "Window dim:" << "\t" << rows << "\t" << cols << "\tPosition:\t" << windowPos << endl;

		for(l=0;l<d;l++)
		{
			signif = 0;
			fOut << setw(2) << l+1 << "\t";
	 		fOut.form("%8.5f", polper(l,0)) << "\t";
	 		fOut.form("%8.5f", rmin(l,0)) << "\t";
	 		fOut.form("%8.5f", rmax(l,0)) << "\t";
	 		fOut.form("%8.5f", polper(l,1)) << "\t";
	 		if( polper(l,0)<rmin(l,0) )
	 			signif = -1;
	 		if( polper(l,0)>rmax(l,0) )
	 			signif = 1;
	 		fOut << signif << endl;
		}
		fOut << "\nTheta Spectrum" << endl;
		fOut << "Randomizations: "<< "\t" << numSimul << "\t" << numExtremeT << "\t" << probConfT << endl;
		for(l=0;l<18;l++)
		{
			signif = 0;
			fOut << setw(4) << l*10 << "\t";
	 		fOut.form("%8.5f", polper(l,2)) << "\t";
	 		fOut.form("%8.5f", thetamin(l,0)) << "\t";
	 		fOut.form("%8.5f", thetamax(l,0)) << "\t";
	 		fOut.form("%8.5f", polper(l,3)) << "\t";
	 		if( polper(l,2)<thetamin(l,0) )
	 			signif = -1;
	 		if( polper(l,2)>thetamax(l,0) )
	 			signif = 1;
	 		fOut << signif << endl;
		}
	}
	else
	{
		int which=2;
		double p=0.95;
		double q=0.05;
		double x=0.00;
		double df = 0,df2=0;
		int status=0;
		double bound=0;

		fOut << "R Spectrum "  << "\t" << outFName << endl;
		if( bonfCorr )
		{
			probConf=probOrig/d;
			fOut << "Two Sides ChiSqr Confidence interval Experiment-wise"
					<< "\t" << probOrig <<  "\t Bonf. corrected\t" << probConf << endl;
		}
		else
			fOut << "Two Sides ChiSqr Confidence interval" << "\t" << probConf << endl;

		fOut << "Data dim:" << "\t" << origRows << "\t" << origCols << endl;
		fOut << "Window dim:" << "\t" << rows << "\t" << cols << "\tPosition:\t" << windowPos << endl << endl;
		for(l=0;l<d;l++)
		{
 			signif = 0;
			fOut << setw(2) << l+1 << "\t";
	 		fOut.form("%8.5f", polper(l,0)) << "\t";

	 		p = probConf/2.0;	// Low confidence limit
	 		q = 1 - p;
	 		df = polper(l,1);
			df2 = df*2;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,0)<x )
	 			signif = -1;
	 		
	 		p = q;					// High confidence limit
			q = 1 - p;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,0)>x )
	 			signif = 1;

	 		fOut.form("%8.5f", df) << "\t";
	 		fOut << signif << endl;
	 		
		}
		fOut << "\nTheta Spectrum" << endl;
		fOut << "Two tail ChiSqr Confidence interval" << "\t" << probConfT << endl;
		for(l=0;l<18;l++)
		{
 			signif = 0;
			fOut << setw(4) << l*10 << "\t";
	 		fOut.form("%8.5f", polper(l,2)) << "\t";

	 		p = probConfT/2.0;		// Low confidence limit
	 		q = 1 - p;
	 		df = polper(l,3);
			df2 = df*2;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,2)<x )
	 			signif = -1;
	 		
	 		p = q;					// High confidence limit
	 		q = 1 - p;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,2)>x )
	 			signif = 1;
	 		
	 		fOut.form("%8.5f", df) << "\t";
	 		fOut << signif << endl;
		}
	}
	
	return 0;
}


void spectMin(double const &pspect, simplmat<double> & pMin, int & l )
{
	if( pspect < pMin(l,0) )
	{
		for(int k=4; k>=0; k--)
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

void spectMax(double const &pspect, simplmat<double> & pMax, int & l )
{
	if( pspect > pMax(l,0) )
	{
		for(int k=4; k>=0; k--)
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

