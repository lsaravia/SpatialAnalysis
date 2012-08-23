#include "mf.h"
#include "RWFile.h"


int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;

   	string fname = argv[1];

	if( argc < 3)
	{
		cerr << "Patch Stats\n";
		cerr << "Usage: PStats inputFile outputFile numSpecies fileType{BI,SP}" << endl;
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
    	if( argc>4)
        {
			if(!file.ReadSeed(argv[1], data, argv[4]))
				exit(1);
        }
        else
			if(!file.ReadSeed(argv[1], data))
				exit(1);
    }

    int numSpecies = atoi( argv[3] );
	PatchStats(data, numSpecies, argv[2] ,argv[1]);

	return 0;
}

