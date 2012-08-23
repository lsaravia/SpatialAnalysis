#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<11)
	{
		cerr << "Usage: sedWin2XY input output type option factor" << endl;
		cerr << "       filter xTopLeft yTopLeft xLength yLength" << endl;
		cerr << " type  : SP/BI/AG/AJ" << endl;
		cerr << " option: 0 Generates an XYZ file " << endl;
		cerr << "         1 Generates an XY file with only positions" << endl;
		cerr << "           nonzero values" << endl;
		cerr << " factor: conversion factor that multiplies coords"	<< endl;
		exit(1);
	}
	simplmat<float> data;
	RWFile file;
	string iname = argv[1];
    string::size_type pos;

	int xTL, yTL, xBR, yBR;
    
    float filter=atof(argv[6]);
    xTL=atoi(argv[7]);
    yTL=atoi(argv[8]);
    xBR=atoi(argv[9]);
    yBR=atoi(argv[10]);

//	if( (pos=iname.find(".img")) != string::npos )
//    {
//		file.ReadIdrisi(argv[1],data);
//    }
//	else
//    {
		if( !file.ReadSeed(argv[1],data,argv[3],xTL,yTL,xBR,yBR) )
            return 1;
//    }
	int option = atoi(argv[4]);
	if( option <0 || option >1)
	{
		cerr << "Error: parameter option outside range" << endl;
		exit(1);
	}
	double factor = atof(argv[5]);
	if( factor<=0 )
	{
		cerr << "Error: parameter factor <= 0 " << endl;
		exit(1);
	}
	
	file.WriteMapXY(argv[2], data, option,factor,filter);
}
