#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<9)
	{
		cerr << "Usage: sedWin input output type"  << endl;
		cerr << "       filter xTopLeft yTopLeft xLength yLength" << endl;
		cerr << " type  : SP/BI/AG/AJ" << endl;
		exit(1);
	}
	simplmat<float> data;
	RWFile file;
	string iname = argv[1];
    string::size_type pos;

	int xTL, yTL, xLen, yLen;
    
    float filter=atof(argv[4]);
    xTL=atoi(argv[5]);
    yTL=atoi(argv[6]);
    xLen=atoi(argv[7]);
    yLen=atoi(argv[8]);

	if( (pos=iname.find(".img")) != string::npos )
    {
		file.ReadIdrisi(argv[1],data);
    }
	else
    {
		if( !file.ReadSeed(argv[1],data,argv[3]) )
            return 1;
    }
	
	file.WriteSeed(argv[2], data, argv[3], xTL,yTL,xLen,yLen);
}
