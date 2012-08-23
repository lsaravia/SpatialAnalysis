// Class Multi Scale Analisis
// Multifractals
// Lacunarity
// 


// Lacunarity
#include "smattpl.h"
#include <iostream>
#include <fstream>
#include "RWFile.h"

int lacunarity(char * filein, char * fileout, double threshold=0.0);

int main(int argc, char * argv[])
{
	if( argc<3 )
	{
		cerr << "Usage: lac2d inputFile outFile" << endl;
		exit(1);
	}
	
	lacunarity(argv[1],argv[2]);
	return 0;
}


int lacunarity(char * filein, char * fileout, double threshold)
{

	int  iwin, ix, iy, icol, irow;
	int  i=0,j;

    double z1,z2,lacun,icnt;

	simplmat <float> pixval;
	
	RWFile file;
   	string fname = filein;

	if( fname.find(".rst") != string::npos )
    {
		if(!file.ReadIdrisi(filein, pixval))
			exit(1);
    }
	else if( fname.find(".tif") != string::npos )
    {
		if(!file.ReadTiff(filein, pixval))
			exit(1);
    }
    else
    {
		if(!file.ReadSeed(filein, pixval))
			exit(1);
    }
   	
   	int xdim = pixval.getRows();
   	int ydim = pixval.getCols();

	ofstream out(fileout);
	if( !out ) {
		cerr << "Cannot open file: " << fileout << endl;
		return 1;
	}
    out << "box\tLacunarity\tlog10(box)\tlog10(Lacunarity)" << endl;
	
	int maxwin = (ydim < xdim ? ydim : xdim) / 2;
	
	double q,sq,s2q;
	
	maxwin++;
	
	// Window size
	for(iwin=1; iwin<maxwin; iwin++)
//	for(iwin=1; iwin<maxwin; iwin*=2)
    {
    	iwin=pow2(i);
    	i++;
		
		z1=0.0;
		z2=0.0;
		
		// Window movement
		for(irow=0; irow<ydim-(iwin-1); irow++)
        {
			for(icol=0; icol<xdim-(iwin-1); icol++)
            {

				icnt=0;
	
				for(iy=irow; iy<irow+iwin; iy++)
                {
					for(ix=icol; ix<icol+iwin; ix++)
                    {
//						if( pixval(ix,iy) > threshold )
						icnt+= pixval(ix,iy);
					}
                }
				q = icnt/((xdim-iwin+1.0)*(ydim-iwin+1.0));
                sq = icnt*q;
                s2q = icnt*icnt*q;
                z1+=sq;
			    z2+=s2q;
								
			}
		}

		lacun=z2/(z1*z1);

		out << iwin << "\t" << lacun << "\t" << log10(static_cast<double>(iwin)) << "\t" << log10(lacun) << endl;
	}

	out.close();

    return 0;

}
