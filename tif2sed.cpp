#include "smattpl.h"
#include "RWFile.h"
#include <math.h>

int main(int argc, char * argv[])
{
	if( argc<4 )
	{
		cerr << "Usage: tif2sed input output grad(0=No) win" << endl;
		exit(1);
	}
	simplmat<float> data;
	RWFile file;
	file.ReadTiff(argv[1], data);
	int dimX = data.getRows();
	int dimY = data.getCols();
	int i,j;	
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
	int win = atoi(argv[4]);   
	if( win>0 )
	{
		int noWinX= dimX/win;
		int noWinY= dimY/win;
		int ii,jj;
		
		simplmat<float> wdat(noWinX,noWinY);;
		float sum=0;
		for(i=0; i<noWinX; i++)
			for( j=0; j<noWinY; j++)
				{
					sum=0;
               int iipri= i*win;
               int jjpri= j*win;
               int iifin= (i+1)*win;
               int jjfin= (j+1)*win;
					for(ii=iipri; ii<iifin; ii++)
						for( jj=jjpri; jj<jjfin; jj++)
							sum+=data(ii,jj);
					wdat(i,j)=sum;
				}
				
		file.WriteSeed(argv[2],wdat);
	}
	else
		file.WriteSeed(argv[2],data);
}
