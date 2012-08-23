#include "mf.h"
#include "RWFile.h"


int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	

   	string fname = argv[1];

	if( argc < 3)
	{
		cerr << "Moran's I Rook\n";
		cerr << "Usage: MIRook inputFile outputFile fileType{BI,SP}" << endl;
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
    	if( argc>=4)
        {
			if(!file.ReadSeed(argv[1], data, argv[6]))
				exit(1);
        }
        else
			if(!file.ReadSeed(argv[1], data))
				exit(1);
    }
   	
   
	MoranIRook(data, argv[2] ,argv[1]);

	return 0;
}

