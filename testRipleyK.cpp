// TestRWFile
// 
#include "RWFile.h"
#include "RipleyK.h"
#include <fstream.h>

int main(int argc, char* argv[])
{
	ifstream fin(argv[1]);
	if( !fin ) 	{
		cerr << "Cannot open file: " << argv[1] << endl;
		return 1;
	}
	
	int dimX,dimY,numAnn=0;
    float x,y,widAnn;
	int count=0;

	dimX = atoi(argv[3]);
	dimY = atoi(argv[4]);
	numAnn = atoi(argv[5]);
	widAnn = atoi(argv[6]);
	
	while( true  )
    {
		fin >> x;
		fin >> y;
		count++;
        if( fin.eof() )
        	break;
    }
    fin.close();
	simplmat <float> ftr(count,2);
    ftr.fill(0.0);

    fin.open(argv[1]);
	if( !fin ) 	{
		cerr << "Cannot open file: " << argv[1] << endl;
		return 1;
	}
	count = 0;
	while( true  )
    {
		fin >> x;
		fin >> y;
		ftr(count,0)=x;
		ftr(count,1)=y;
		count++;		
        if( fin.eof() )
        	break;
    }
    fin.close();
    
	RipleyK("test.out", ftr, numAnn, widAnn, count, dimX, dimY );
    
}
