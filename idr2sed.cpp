#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<4 )
	{
		cerr << "Usage: img2sed input output grad(0=No) [cutOff(0=No)] [cutDown(0=No)] " << endl;
		cerr << "[Values greater than cutOff are set to 0]" << endl;
		cerr << "[Values lesser or equal than cutDown are set to 0]" << endl;

		exit(1);
	}
	simplmat<float> data;
   
	RWFile file;
	file.ReadIdrisi(argv[1], data);
    int dimX = data.getRows();
	int dimY = data.getCols();
    int i,j;

	// Pone a 0 valores mayores a cutoff
	// 
	if( argc==5 )
	{
		int cutoff= atoi(argv[4]);
		if(cutoff>0)
		{
			for(i=0; i<dimX; i++)
				for( j=0; j<dimY; j++)
				{
					if(data(i,j)>cutoff)
						data(i,j)=0;
				}
		}
	}

	// Pone a 0 valores menores o iguales a cutoff
	// 
	if( argc==6 )
	{
		int cutoff= atoi(argv[5]);
		if(cutoff>0)
		{
			for(i=0; i<dimX; i++)
				for( j=0; j<dimY; j++)
				{
					if(data(i,j)<=cutoff)
						data(i,j)=0;
				}
		}
	}

	// Convierte a Gradiente
	//
	//
	int grad = atoi(argv[3]);
	if( grad)
	{
	   simplmat<float> ndat(dimX,dimY);
	   float xx,yy;
	   for(i=0; i<dimX; i++)
			for( j=0; j<dimY; j++)
			{
				if(i==dimX-1)
					xx = data(i, j) - data(i-1,j);
				else if(i==0)
				  xx = data(i+1, j) - data(i,j);
			else
					xx = data(i+1, j) - data(i-1,j);
					
				if(j==dimY-1)
					yy = data(i, j) - data(i,j-1);
				else if(j==0)
					yy = data(i, j+1) - data(i,j);
			else
				yy = data(i, j+1) - data(i,j-1);
				
				ndat(i,j)=sqrt( xx*xx + yy*yy );
		}
		for(i=0; i<dimX; i++)
			for( j=0; j<dimY; j++)
				data(i,j)=ndat(i,j);

   }

	file.WriteSeed(argv[2],data);
}
