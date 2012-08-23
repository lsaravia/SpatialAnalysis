#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<6 )
	{
		cerr << "Usage: sed2XY input output type option factor" << endl;
		cerr << " type  : SP/BI/AG/AJ/JU" << endl;
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
	string opt = argv[3];
	if( opt=="JU" )
		opt="AJ";
		
	if( (pos=iname.find(".img")) != string::npos )
    {
		file.ReadIdrisi(argv[1],data);
    }
	else
    {
		if(!file.ReadSeed(argv[1],data,opt.c_str()))
		{
			cerr << "Error " << argv[1] << endl;
			exit(1);
		}
    }
    opt = argv[3];
    if( opt=="JU" )
    {
    	int dx=data.getRows();
    	int dy=data.getCols();
    	for(int x=0; x<dx; x++)
    		for(int y=0; y<dy; y++)
    		{
    			int dd = data(x,y);
    			if( (dd%2)!=0)
    				data(x,y) = 0;
    		}
    }
    
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
	
	file.WriteMapXY(argv[2], data, option,factor);
}
