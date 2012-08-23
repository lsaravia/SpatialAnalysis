#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<8 )
	{
		cerr << "Usage: XY2Sed input output xsize xstep ysize ystep option" << endl;
		cerr << "    Option 0: Reads only positions, format X Y" << endl;
		cerr << "    Option 1: Reads positions and value, Format X Y val" << endl;
		cerr << "    Option 2: Reads positions and value, and value2 to decide ties" << endl;
		cerr << "              If there are several ties the record with greater value2 is included" << endl;
		exit(1);
	}
    float xsize=atof(argv[3]);
    float xstep=atof(argv[4]);
    float ysize=atof(argv[5]);
    float ystep=atof(argv[6]);
	int option = atoi(argv[7]);
	simplmat<float> data;
	RWFile file;
	file.ReadMapXY(argv[1], xsize,xstep,ysize,ystep, data,option);
	file.WriteSeed(argv[2],data);
}
