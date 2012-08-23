#include "Spectral.h"
#include "RWFile.h"
#include "Randomizations.h"
#include <iomanip>
#include "fortify.h"
#include "ctype.h"
#include "cdflib.h"

int main(int argc, char * argv[])
{
	simplmat <double> data;
	simplmat <double> dout;
	SpectralAnalysis sa;
	
	// By default does 199 Simulations and takes the 5 lowest and 5 highest values
	// to make the confidence envelopes. If the calculated value is
	// between these highest or lowest values it is significative at 5% level
	//
	int numSimul=199,numExtreme=5;
	char randz;
	double probConf=0.05;
	int bonfCorr=0;
	string fType="BI";
	
	
	if( argc <7 )
    {
		cerr << "Usage: rnzSpectral inputFile.sed outFile [R/A] [numSimul/prob] [numExtreme/bonfCorr] [Type SP/AG/BI/AJ]" << endl;
        exit(1);
	}

	randz = toupper(argv[3][0]); // R=Randomizations, A=Asintotic chisqr.

	if( randz=='R' )
	{
		numSimul = atoi(argv[4]);
		numExtreme = atoi(argv[5]);
	}
	else
	{
		probConf = atof(argv[4]);
		bonfCorr = atoi(argv[5]);
	}
	fType = argv[6];
			

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
	
//	file.ReadSeed(argv[1], data,argv[6]);
    sa.Transpose(data); // Los datos leidos estan en formato (x,y) y las funciones
						// tienen (y,x) o sea (row,col)
	int rows=data.getRows();
	int cols=data.getCols();
	if( (rows % 2)!=0 || (cols % 2)!=0)
	{
		cerr << "Error, dimensions must be multiplo of 2" << endl;
		exit(1);
	}
	
	double var;
	var = sa.Period2D(data,dout); // Calculates the periodogram

	simplmat<double> rper;		// Rearranged periodogram
	simplmat<double> polper;	// Polar
	sa.Spekout(rows,cols,dout,rper); // Rearranges the periodogram

	string rperOut=outFName + ".rpe";
	ofstream fOut(rperOut.c_str());
	fOut <<	rper.getCols() << "\t" << rper.getRows() << endl;
	fOut << "Rearranged periodogram" << endl;
	for(int l=0;l<rper.getRows();l++)
	{
		for(int s=0;s<rper.getCols();s++)
			fOut.form("%10.6f",rper(l,s)) <<  "\t";
		fOut << endl;
	}
	fOut.close();
	
	sa.Polar2D(rows,cols,rper,polper); // Calculates the Polar Spectrum

	int d =	int(0.5*sqrt(rows*rows + cols*cols)+1);
	int s,l;

	fOut.open(outFName.c_str());

	
	if(	randz=='R' )
	{
		Randomizations rz;
		simplmat<double> thetamin;
		simplmat<double> thetamax;
		simplmat<double> rmin;
		simplmat<double> rmax;
		simplmat<double> rpol;		// Polar Spectrum for randomizations

		thetamin.resize(18,numExtreme,1000.0);
		thetamax.resize(18,numExtreme,0.0);
		rmin.resize(d,numExtreme,1000.0);
		rmax.resize(d,numExtreme,0.0);
	
		for(s=0; s<numSimul; s++)
		{
			rz.Randomize(data);
			sa.Period2D(data,dout);
			sa.Spekout(rows,cols,dout,rper);
			sa.Polar2D(rows,cols,rper,rpol);
			for(l=0;l<d;l++)
			{
				sa.spectMax(rpol(l,0),rmax,l, numExtreme);
				sa.spectMin(rpol(l,0),rmin,l, numExtreme);
	
	//			fOut << setw(2) << l+1 << "\t";
	// 			fOut.form("%8.5f", rpol(l,0)) << endl;
			}
			for(l=0;l<18;l++)
			{
				sa.spectMax(rpol(l,2),thetamax,l, numExtreme);
				sa.spectMin(rpol(l,2),thetamin,l, numExtreme);
	//			fOut << "\t\t" << setw(4) << l*10 << "\t";
	//	 		fOut.form("%8.5f", rpol(l,2)) << endl;
			}
		}
		fOut << "\nRandomizations: " << numSimul << "\t" << numExtreme << endl;
		fOut << "Data dim:" << "\t" << rows << "\t" << cols << endl << endl;
		fOut << "R Spectrum" << endl;

		for(l=0;l<d;l++)
		{
			fOut << setw(2) << l+1 << "\t";
	 		fOut.form("%8.5f", polper(l,0)) << "\t";
	 		fOut.form("%8.5f", rmin(l,0)) << "\t";
	 		fOut.form("%8.5f", rmax(l,0)) << "\t";
	 		fOut.form("%8.5f", polper(l,1)) << endl;
		}
		fOut << "\nTheta Spectrum" << endl;
		for(l=0;l<18;l++)
		{
			fOut << setw(4) << l*10 << "\t";
	 		fOut.form("%8.5f", polper(l,2)) << "\t";
	 		fOut.form("%8.5f", thetamin(l,0)) << "\t";
	 		fOut.form("%8.5f", thetamax(l,0)) << "\t";
	 		fOut.form("%8.5f", polper(l,3)) << endl;
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
		int signif=0;
		double probOrig=probConf;
		if( bonfCorr )
		{
			probConf=probOrig/d;
			fOut << "Two Sides ChiSqr Confidence interval Experiment-wise"
					<< "\t" << probOrig <<  "\t Bonf. corrected\t" << probConf << endl;
		}
		else
			fOut << "Two Sides ChiSqr Confidence interval" << "\t" << probConf << endl;
		fOut << "Data dim:" << "\t" << rows << "\t" << cols << endl << endl;
		fOut << "R Spectrum " << endl;
		for(l=0;l<d;l++)
		{
			signif=0;
			fOut << setw(2) << l+1 << "\t";
	 		fOut.form("%8.5f", polper(l,0)) << "\t";

	 		p = probConf/2.0;	// Low confidence limit 0.05/2
	 		q = 1 - p;
	 		df = polper(l,1);
			df2 = df*2;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,0)<x )
	 			signif = -1;
	 		
	 		p = q;				// High confidence limit
			q = 1 - p;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";

	 		fOut.form("%8.5f", df) << "\t";
	 		if( polper(l,0)>x )
	 			signif = 1;
	 		fOut << signif << endl;
		}
			
		fOut << "\nTheta Spectrum\t";
		if( bonfCorr )
		{
			probConf=probOrig/18;
			fOut << "Experiment-wise" << "\t" << probOrig <<  "\tBonf. corrected\t" << probConf << endl;
		}
		else
			fOut << endl;
		for(l=0;l<18;l++)
		{
			signif=0;
			fOut << setw(4) << l*10 << "\t";
	 		fOut.form("%8.5f", polper(l,2)) << "\t";

	 		p = probConf/2;		// Low confidence limit 0.05/2
	 		q = 1 - p;
	 		df = polper(l,3);
			df2 = df*2;
			cdfchi(&which,&p,&q,&x,&df2,&status,&bound);
			x = 1/df2 * x;
	 		fOut.form("%8.5f", x) << "\t";
	 		if( polper(l,2)<x )
	 			signif = -1;
	 		
	 		p = q;									// High confidence limit
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



