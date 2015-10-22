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
	simplmat <double> xydata;
	simplmat <double> dout;
	SpectralAnalysis sa;
	
	char randz;
	double probConf=0.05;
	int bonfCorr=0,rows=0,cols=0;
	float lx=0,ly=0;
	
	if( argc <3 )
    {
		cerr << "Usage: SpectralEvent inputFile.sed outFile A prob bonfCorr" << endl;
		cerr << "       SpectralEvent inputFile.dat outFile A prob bonfCorr longx longy cantValues" << endl;
        exit(1);
	}
	if( argc >= 4)
	{
		randz = toupper(argv[3][0]); // R=Randomizations, A=Asintotic chisqr.
	}
	if( randz=='R' )
	{
		exit(1);
	}
	else
		if( argc >= 6 )
		{
			probConf = atof(argv[4]);
			bonfCorr = atoi(argv[5]);
		}
			

	RWFile file;
	string fName = argv[1];
	string outFName = argv[2];
	
	if( fName.find(".sed")!=string::npos )
	{
		if(!file.ReadSeed(fName.c_str(), data))
			exit(1);
		file.Conv2XY(data,xydata);
		rows=data.getRows();
		cols=data.getCols();
		lx=cols;
		ly=rows;
	}
	else if( fName.find(".dat")!=string::npos )
	{
		if(!file.ReadXYVec(fName.c_str(), xydata))
			exit(1);
		lx = atof(argv[6]);
		ly = atof(argv[7]);
		rows = atoi(argv[8]);
		cols = rows;
	}
	

	if( (rows % 2)!=0 || (cols % 2)!=0)
	{
		cerr << "Error, dimensions must be multiplo of 2" << endl;
		exit(1);
	}
	
	sa.Period2DEvent(xydata,lx,ly,rows,dout); // Calculates the periodogram

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
	
	//sa.Polar2DEvent(xydata.getRows(),rper,polper); // Calculates the Polar Spectrum
	sa.Polar2D(rows,cols,rper,polper); // Calculates the Polar Spectrum
									   // Ojo modifica rper
									   
	int d =	int(0.5*sqrt(rows*rows + cols*cols)+1);
	int s,l;

	fOut.open(outFName.c_str());

	
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



