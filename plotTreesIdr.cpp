// Class Multi Scale Analisis
// Multifractals
// Lacunarity
// 


// Lacunarity
#include "smattpl.h"
#include "RWFile.h"
#include <iostream>
#include <fstream>
#include <string>

int plotTreesIdr(char * filein, char * fileout, int dimX,int dimY);

int main(int argc, char * argv[]) {

	if( argc <5 )
    {
		cerr << "Usage: input output xdim ydim " << endl;
        exit(1);
	}
        
	int xd = atoi(argv[3]);
	int yd = atoi(argv[4]);
	plotTreesIdr(argv[1],argv[2],xd,yd);

    return 0;
	}


int plotTreesIdr(char * filein, char * fileout, int dimX,int dimY) {

	ifstream fin(filein);
	if( !fin ) 	{
		cerr << "Cannot open file: " << filein << endl;
		return 1;
	}

    dimX++;
    dimY++;
	simplmat <int> pixval(dimX,dimY);
    pixval.fill(0);
	int dx,dy,ival=0;
    float x,y;
	string adjuv;
	float val;

	while( true  )
    {
		fin >> x;
		fin >> y;

		y = dimY - y;
		// De esta manera cuando lo muestro con el IDRISI
		// me queda el rio en la parte superior (ESTE) y el NORTE en la
		// parte izquierda de la pantalla
		
		fin >> val;
        fin >> adjuv;
        if( fin.eof() )
        	break;
            
		cerr << x << "\t" << y << "\t" << val << "\t" << adjuv << endl;
		if( x<dimX && y<dimY )
		{
	       	ival = val;
	       	
	       	if (adjuv[0]=='A')
	       		ival = val / 2;
       		else
	       		ival = val * 2;
	       	
	       	
			for( dx=x-ival; dx<=x+ival; dx++)
				for( dy=y-ival; dy<=y+ival; dy++)
	            {
	            	if( dx>=0 && dy>=0  && dx<dimX && dy<dimY )
                    {
						if( adjuv == "JLIG" )
							pixval(dx,dy) = 1;
						if( adjuv == "ALIG" )
							pixval(dx,dy) = 1;
						if( adjuv == "JTAL" )
							pixval(dx,dy) = 2;
						if( adjuv == "ATAL" )
							pixval(dx,dy) = 2;
						if( adjuv == "JCOR" )
							pixval(dx,dy) = 3;
						if( adjuv == "ACOR" )
							pixval(dx,dy) = 3;
						if( adjuv == "JJOD" )
							pixval(dx,dy) = 4;
						if( adjuv == "AJOD" )
							pixval(dx,dy) = 4;
						if( adjuv == "JSAU" )
							pixval(dx,dy) = 5;
						if( adjuv == "ASAU" )
							pixval(dx,dy) = 5;
						if( adjuv == "JLIN" )
							pixval(dx,dy) = 6;
						if( adjuv == "ALIN" )
							pixval(dx,dy) = 6;
                            
                    }
                        
	            }
                            
		}
		else
			cerr << "coord out of range " << x << "  " << y << "  " << val << endl; 
		
    }

	fin.close();

	RWFile file;

	file.WriteIdrisi( fileout, pixval);


	return 0;
}
