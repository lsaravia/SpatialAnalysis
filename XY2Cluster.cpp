#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<7 )
	{
		cerr << "Usage: XY2Cluster input output xsize ysize localScale readOption outOption" << endl;
		cerr << "    readOption 0: Reads only positions, format X Y" << endl;
		cerr << "    			1: Reads positions and a species, Format X Y Species" << endl;
		cerr << "    			2: Reads positions,species, and a value that may be basal area (not used)" << endl << endl;
		cerr << "    outOption 0: count the abundance in each site an puts the most abundant " << endl;
		cerr << "    			anyNumber: outputs the species anyNumber " << endl << endl;

		exit(1);
	}
    float xsize=atof(argv[3]);
    float ysize=atof(argv[4]);
    float localScale=atof(argv[5]);
	int option = atoi(argv[6]);
	int outOption = atoi(argv[7]);
	simplmat<float> data;
	RWFile file;
	file.ClusterizeXY(argv[1], xsize,ysize,localScale, data,option,outOption);
	file.WriteSeed(argv[2],data,"SP");  // data type SP species number 1
	return(0);
}
