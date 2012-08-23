#include <math.h>
#include <algorithm>
#include "RWFile.h"
#include "mf.h"
#include "cdflib.h"
using namespace std;

int MoranIRook(simplmat <double> &data, const char * outFile, const char * ident)
{
    double bioSum=0,bioMean=0,bioVar=0,mor1=0,moran=0,z=0,z4=0,mor3=0;
	double s1=0,s2=0,moranVar=0,b2=0,n=0,W=0,W2=0,n2=0;
	double mv1,mv2,mv3,mv4;
	int x,y;
	int fmprivez=0;


	fstream fm( outFile, ios::in );
	if(!fm )
		fmprivez=1;
	fm.close();
	fm.open( outFile, ios::app );
	if(!fm)
	{
		cerr << "Cannot open Moran's I file, " << outFile << endl;
		return 0;
	}

	if( fmprivez )
		fm << "File\t" << "Moran's I\t" << "Var.\t" << "Expected\t" << "Confidence Interval" << endl;

	int xDim = data.getRows();
	int yDim = data.getCols();

	for( x=0; x<xDim; x++ )
		for( y=0; y<yDim; y++ )
		{
			bioSum += data(x,y);
		}

	n = xDim*yDim;
	bioMean = bioSum/n;

	for(x=0; x<xDim; x++)
		for(y=0; y<yDim; y++)
		{
			z = data(x,y) - bioMean;
			bioVar += z*z;
			z4 = z*z*z*z;
			mor3 += z4;

			if(y<yDim-1)
			{
				mor1 += z*(data(x,y+1)-bioMean);
			}
			if(y>0)
			{
				mor1 += z*(data(x,y-1)-bioMean);
			}
			if(x<xDim-1)
			{
				mor1 += z*(data(x+1,y)-bioMean);
			}
			if(x>0)
			{
				mor1 += z*(data(x-1,y)-bioMean);
			}
		}
	W  = 8+6*(yDim-2)+6*(xDim-2)+(xDim-2)*(yDim-2)*4;
	s2    = 4*16 + 36*(xDim-2+yDim-2) + (xDim-2)*(yDim-2)*64;
	s1    = 2*W;
	moran = n*mor1/(W*bioVar);
	b2    = n*mor3/(bioVar*bioVar);
	n2	   = n*n;
	W2    = W*W;

	mv1 = n*((n2-3*n+3)*s1-n*s2+3*W2);
	mv2 = b2*((n2-n)*s1-2*n*s2+6*W2);
	mv3 = (n-1)*(n-2)*(n-3)*W2;
	mv4 = 1.0/((n-1)*(n-1));
	moranVar = (mv1-mv2)/mv3 - mv4;

	fm << ident << "\t" << moran << "\t" << moranVar << "\t" << -1.0/(n-1);


    // Falta poner como parametro el valor de alfa = q
    //
    int which=2,status=0;
    double p=.025;
    double q=1-p;
    double mean=0, sd=1, bound=0;
    z=0;

    cdfnor(&which,&p,&q,&z,&mean,&sd,&status,&bound);

    mv3 = moran - z * sqrt(moranVar);
    mv4 = moran + z * sqrt(moranVar);
    fm << "\t" << mv4 << "\t" << mv3;

	mv4 =  z*sqrt(moranVar) - 1.0/(n-1);
	mv3 = -z*sqrt(moranVar) - 1.0/(n-1);
	if(moran>mv4 || moran<mv3)
		fm << "\t" << "*" << endl;
	else
		fm << endl;

	return 1; // NO ERROR
}
