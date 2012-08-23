#include "mf.h"
#include "RWFile.h"


int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	
	simplmat <double> q;	

   	string fname = argv[1];

	if( argc < 6)
	{
		cerr << "Pair Correlation\n";
		cerr << "Usage: pCorr inputFile outputFile minDist maxDist deltaDist fileType{BI,SP}" << endl;
		exit(1);
	}
	
	if( fname.find(".rst") != string::npos )
    {
		if(!file.ReadIdrisi(argv[1], data))
			exit(1);
    }
	else if( fname.find(".tif")!=string::npos )
	{
		if(!file.ReadTiff(argv[1], data))
			exit(1);
	}
	else
    {
    	if( argc>=7)
        {
			if(!file.ReadSeed(argv[1], data, argv[6]))
				exit(1);
        }
        else
			if(!file.ReadSeed(argv[1], data))
				exit(1);
    }
   	
   
	float minDist = atof(argv[3]);
	float maxDist = atof(argv[4]);
	float deltaDist = atof(argv[5]);
    
    int maxPoints=0;
    if( argc==8 )
        maxPoints = atoi(argv[7]);
        
    
	PairCorrelation(data, argv[2] ,minDist, maxDist, deltaDist, maxPoints);

	return 0;
}

