#include <fstream.h>
#include <iomanip>
#include <strstream>
#include "mf.h"
using namespace std;

struct patchXY{
		int x;
		int y;
		char check;
		patchXY(){ x=0;y=0;check=0; };
		patchXY & operator=(const patchXY &src){ x = src.x;  y = src.y; check=src.check; return *this;};
		patchXY(const patchXY &src){ x = src.x;  y = src.y; check=src.check; };
		};

static int PatchAux(simplmat <double> &data, int &xDim,int &yDim, int &x,int &y )
	{
	int ret=0;
	if (x>=xDim || x<0 )
		{
		ret=1;
		x= x<0 ? 0 : xDim;
		}
	if (y>=yDim || y<0 )
		{
		ret=1;
		y= y<0 ? 0 : yDim;
		}
	if( ret==1)
		return 0;
	else
		return data(x,y);
	}
	
int PatchStats(simplmat <double> &data,int numSpecies, const char * outFile,const char * ident)
{
	int x,xx,dx,y,yy,dy,e,i;
	long patSize;

	int xDim = data.getRows();
	int yDim = data.getCols();

	simplmat <unsigned long> pSize(numSpecies);
	simplmat <unsigned long> pNumb(numSpecies);
	simplmat <unsigned long> pMaxs(numSpecies);
	simplmat <patchXY> pxy(xDim*yDim);
	pSize.fill(0);
	pNumb.fill(0);
	pMaxs.fill(0);
	pxy.fill( patchXY() );

	fstream fb;
	int privez=0;

	fb.open( outFile, ios::in );
	if(!fb )
		privez=1;
	fb.close();
	fb.open( outFile, ios::app );
	if( !fb )
		{
		cerr << "Cannot open patch stats file: " << outFile << endl;
		return 0;
		}

	if( privez )
		fb << "File\tClass\tT.Size\tAvg.\tNum.\tMax.Size" << endl;

	for(x=0; x<xDim; x++)
		{
		for(y=0; y<yDim; y++)
			{
			e=data(x,y);
			if( e > 0 )
				{
                if( e>numSpecies )
                {
            		cerr << "Invalid number of categories in data\n";
		            return 0;
                }
				pNumb(e-1)++;
				patSize=0;
				xx=x;
				yy=y;
				pxy(patSize).x=xx;
				pxy(patSize).y=yy;
				pxy(patSize).check=0;
				data(x,y)=0;
				while(1)
					{
					for(i=0; i<=patSize; i++)
						{
						if(pxy(i).check==0)
							{
							xx = pxy(i).x;
							yy = pxy(i).y;
							pxy(i).check=1;
							dx = xx+1;
							dy = yy;
							if( e==PatchAux(data,xDim,yDim,dx,dy) )
								{
								data(dx,dy)=0;
								patSize++;
								pxy(patSize).x=dx;
								pxy(patSize).y=dy;
								pxy(patSize).check=0;
								}
							dx = xx-1;
							dy = yy;
							if( e==PatchAux(data,xDim,yDim,dx,dy) )
								{
								data(dx,dy)=0;
								patSize++;
								pxy(patSize).x=dx;
								pxy(patSize).y=dy;
								pxy(patSize).check=0;
								}
							dx = xx;
							dy = yy-1;
							if( e==PatchAux(data,xDim,yDim,dx,dy) )
								{
								data(dx,dy)=0;
								patSize++;
								pxy(patSize).x=dx;
								pxy(patSize).y=dy;
								pxy(patSize).check=0;
								}
							dx = xx;
							dy = yy+1;
							if( e==PatchAux(data,xDim,yDim,dx,dy) )
								{
								data(dx,dy)=0;
								patSize++;
								pxy(patSize).x=dx;
								pxy(patSize).y=dy;
								pxy(patSize).check=0;
								}

							}
						}
					pSize(e-1) += (++patSize);
					if( patSize > pMaxs(e-1) )
						pMaxs(e-1) = patSize;

					// Para distribucion de tamanios de parche
					// fb << ident << "\t" << patSize << endl
					//

					break;
					}
				}
			}
		}
	for( x=0; x<numSpecies; x++)
	{
		double tcells=xDim*yDim;
		if (pNumb(x)==0)
			fb << ident << "\t" << (x+1) << "\t"
			<< pSize(x)/tcells*100 << "\t"
            << pNumb(x)	<< "\t"
            << pNumb(x) << "\t"
            << pMaxs(x)/tcells*100 << endl;
		else
			fb << ident << "\t" << (x+1) << "\t"
			<< pSize(x)/tcells*100 << "\t"
			<< pSize(x)/double(pNumb(x))/tcells*100 << "\t"
			<< pNumb(x) << "\t"
			<< pMaxs(x)/tcells*100 << endl;
	}
		
	return 1; // No ERROR
}



