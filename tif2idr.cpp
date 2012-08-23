#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<3 )
	{
		cerr << "Usage: tif2idr input output" << endl;
		exit(1);
	}
	simplmat<float> data;
	RWFile file;
	file.ReadTiff(argv[1], data);
	file.WriteIdrisi(argv[2],data);
   
}
