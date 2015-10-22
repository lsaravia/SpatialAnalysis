#include "Spectral.h"
#include "RWFile.h"
#include <iomanip>
#include "fortify.h"

int main(int argc, char * argv[])
{
	simplmat <double> data;
	simplmat <double> dout;
	SpectralAnalysis sa;
	
	if( argc <2 )
    {
		cerr << "Usage: Spectral inputFile.sed " << endl;
        exit(1);
	}
	RWFile file;
	file.ReadSeed(argv[1], data);
   
   sa.Transpose(data); // Los datos leidos estan en formato (x,y) y las funciones
						// tienen (y,x) o sea (row,col)
	int rows=data.getRows();
	int cols=data.getCols();
	double var;
	
	var = sa.Period2D(data,dout); // Calculates the periodogram

	int xd=dout.getRows();
	int yd=dout.getCols();

	cout << "Data variance: " << var << endl;
	cout << "Periodogram I(j,k) for j=0,...,m/2 and k=0,...,n-1" << endl << endl;

	int i,j;
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%10.6f",dout(i,j)) << "\t";
			
		cout << endl;
	}

	var = sa.Spec2D(data,dout); // Calculates the periodogram

	int x1=dout.getRows();
	int y1=dout.getCols();
	
	cout << "Calculate the periodogram using Fast Fourirer transform" << endl;
	cout << "Data variance: " << var << endl;
	cout << "Periodogram I(j,k) for j=0,...,m/2 and k=0,...,n-1" << endl << endl;
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%10.6f",dout(i,j)) << "\t";
			
		cout << endl;
	}



	sa.Spekout(rows,cols,dout,data); // Rearranges the periodogram
								
	cout << endl << endl << "Rearranged periodogram I(j,k) with j=0,...,m/2 and k=-n/2,...,n/2-1." <<	endl;
	xd=data.getRows();
	yd=data.getCols();
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%10.6f",data(i,j)) << "\t";
			
		cout << endl;
	}

	simplmat<double> rper(data);
	
	sa.Spekout(rows,cols,dout,data,var); // Rearranges the periodogram and express as a percentage of var
	cout << endl << endl << "Rearranged periodogram as a percentage of data variance" <<	endl;
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%5.2f",data(i,j)) << "\t";
			
		cout << endl;
	}

	double sens = 400.0/static_cast<double>(rows*cols);
	var = sa.Spekout(rows,cols,dout,data,var,sens); // Rearranges the periodogram and express as a percentage of var
	cout << endl << endl << "Rearranged periodogram as percentage of data variance" <<	endl;
	cout << "Sensoring value: " << sens << endl;
	cout << "Variance explained: " << var << endl;
	
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%4.0f",data(i,j)) << "\t";
			
		cout << endl;
	}

	sa.Polar2D(rows,cols,rper);

    return 0;
}



