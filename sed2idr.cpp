#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<4 )
	{
		cerr << "Usage: sed2img input Type output" << endl;
		exit(1);
	}
	simplmat<float> data;
	RWFile file;
	if(!file.ReadSeed(argv[1],data,argv[2]))
		{
			cerr << "Error " << argv[1] << endl;
			exit(1);
		}
	
	file.WriteIdrisi(argv[3],data);
}
