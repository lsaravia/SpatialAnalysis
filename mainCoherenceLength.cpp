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
		cerr << "Coherence length\n";
		cerr << "Usage: cl inputFile outputFile minBox maxBox deltaBox" << endl;
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
		if(!file.ReadSeed(argv[1], data))
			exit(1);
    }
   	
   
	int minBox = atoi(argv[3]);
	int maxBox = atoi(argv[4]);
	int deltaBox = atoi(argv[5]);
    
	CoherenceLength(data, argv[2] ,minBox, maxBox, deltaBox);

	return 0;
}

