/*  Copyright 2011 Leonardo A. Saravia
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/
#ifndef RWFILE_H
#define RWFILE_H
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "smattpl.h"
#ifdef __linux
#include "tiffio.h"
#endif
using namespace std;

class RWFile
{
	unsigned DimX;
	unsigned DimY;

	public:
	RWFile() {DimX=0;DimY=0;};
	~RWFile() {};

	template <class Type> bool WriteMapXY(const char * fname, simplmat<Type>& data,int option=0, double factor=1, Type filter=0);
	template <class Type> bool Conv2XY(simplmat<Type>& data,simplmat<Type>& out);
	template <class Type> bool ReadSeed(const char * fname,simplmat<Type>& data, const char * dataType="");
	template <class Type> bool ReadSeed(const char * fname,simplmat<Type>& data, const char * dataType,int xTL, int yTL, int xLen, int yLen);
	template <class Type> bool WriteSeed(const char * fname,simplmat<Type>& data, const char * dataType="");
	template <class Type> bool WriteSeed(const char * fname,simplmat<Type>& data, const char * dataType, int xTL, int yTL, int xLen, int yLen);
	template <class Type> bool ReadIdrisi(const char * fname, simplmat<Type>& data);
	template <class Type> bool ReadMapXY(const char * finp, float xsize,float xstep,float ysize,float ystep, simplmat<Type> &data,int option=0);
	template <class Type> bool ReadXYVec(const char * finp, simplmat<Type> &data,int option=0);
    template <class Type> bool ClusterizeXY(const char * finp, int xsize,int ysize,float lc, simplmat<Type> &data, int readOption, int outOption);
	template <class Type> bool ReadTiff(const char * fname, simplmat<Type>& data);

};


//
// Genera un archivo del tipo "X Y Value", agrega la extencion ".map" al nombre
// y el numero total de puntos al pricipio del archivo.
//
//
template <class Type> bool RWFile::WriteMapXY(const char * fname, simplmat<Type>& data, int option, double factor, Type filter)
{
	int i,j;
	Type sp;
	long tot=0;

    string::size_type pos=0;

	string name = fname;
	if( (pos=name.find(".")) == string::npos )
    {
		name += ".map";
    }
	string dname = name + ".txt";
	
	ofstream xymap( dname.c_str() );

    DimX=data.getRows();
    DimY=data.getCols();
	
	for(i=0; i<DimY; i++)
		for( j=0;j<DimX;j++)
			if( data(j,i) != 0 )
				tot++;

	xymap << "DimX: " << DimX << endl;
	xymap << "DimY: " << DimY << endl;
	xymap << "No. Points != 0: " << tot << endl;
	xymap.close();
	
	xymap.open( name.c_str() );
	switch(option)
	{
	case 0:
		// File X Y Z
		//
    	for( i=0; i<DimY; i++)
    		for( j=0;j<DimX;j++)
    		{
	   			sp = data(j,i);
	   			if( filter == 0 )
	   			{
					if(sp!=0)
		    			xymap << j*factor << "\t" << i*factor << "\t" << sp << endl;
	   			}
		    	else
		    	{
					if(sp==filter)
		    			xymap << j*factor << "\t" << i*factor << "\t" << sp << endl;
		    	}
    		}
    	break;
	case 1:
		// File X Y prescence/auscence only
		//
    	for( i=0; i<DimY; i++)
    		for( j=0;j<DimX;j++)
    		{
    			sp = data(j,i);
	   			if( filter == 0 )
	   			{
					if( sp!=0 )
//	    				xymap.form("%7.2f%7.2f",j*factor,i*factor) << endl;
	    				xymap << j*factor << "\t" << i*factor << endl;
	   			}
	    		else
		    	{
					if(sp==filter)
		    			xymap << j*factor << "\t" << i*factor << "\t" << sp << endl;
		    	}
    		}
    	break;
	}
return true;
}

//
// Convierte una matriz tipo sed a un vector de coordenadas
//
template <class Type> bool RWFile::Conv2XY(simplmat<Type>& data,simplmat<Type>& out)
{
	int i,j;
	Type sp;
	long t=0;

    DimX=data.getRows();
    DimY=data.getCols();

	for(i=0; i<DimY; i++)
		for( j=0;j<DimX;j++)
			if( data(j,i) != 0)
				t++;
	out.resize(t,2);

	// File X Y prescence/auscence only
	//
	t=0;
	for( i=0; i<DimY; i++)
		for( j=0;j<DimX;j++)
		{
			sp = data(j,i);
			if( sp != 0 )
			{
				out(t,0)=j;
				out(t,1)=i;
				t++;
			}
		}
	return true;
}




template <class Type> bool RWFile::ReadSeed(const char * fname,simplmat<Type>& data, const char * dataType)
{
	ifstream in;
	string buff,dType;
	unsigned int dx,dy;
	Type spe;

	in.open(fname);
	if ( !in )
	{
		cerr << "Cannot open file: " << fname << endl;
		return false;
	}
	in >> dx;
	in >> dy;


	if( dx!=data.getRows() || dy!=data.getCols() )
		data.resize(dx,dy);
		
	DimX=dx;
	DimY=dy;

	dType = dataType;
	// if dataType is empty assing BI type
	if( dType.empty() )
		dType = "BI";
	
	while( !in.eof() )
	{
		getline(in,buff);
		if ( buff.find(dType)==string::npos )
                {
                    if(in.eof()){		// No encontro el tipo especificado
                        cerr << "Type " << dType << " not found.\n";
                        return false;
                    }
                
        	continue;
                }
        
            
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
				in >> spe;
				if( !in )
				{
					cerr << "Seed File invalid.\n";
					return  false;
				}

				data(dx,dy) = spe;

			}
      	break;
	}
	return true;
}

template <class Type> bool RWFile::ReadSeed(const char * fname,simplmat<Type>& data, const char * dataType,
			int xTL, int yTL, int xLen, int yLen)
{
	ifstream in;
	string buff,dType;
	int dx,dy,ddx,ddy;
	int tipo;

	Type spe;

	in.open(fname);
	if ( !in )
	{
		cerr << "Cannot open file: " << fname << endl;
		return false;
	}
	in >> DimX;
	in >> DimY;
	
	int xBR = xTL + xLen;
	int yBR = yTL + yLen;

	if( xTL<0 || yTL<0 || xLen<0  || yLen<0 || xBR>DimX  || yBR>DimY)
	{
		cerr << "Invalid window" << endl;
		return false;
	}

	if( xLen!=data.getRows() || yLen!=data.getCols() )
		data.resize(xLen,yLen);
		
	dType = dataType;
	if( dType.empty() )
		dType = "BI";

	while( !in.eof() )
	{
		getline(in,buff);
		if ( buff.find(dType)==string::npos )
        {
        	if(in.eof())		// No encontro el tipo especificado
            	return false;
                
        	continue;
        }
        
            
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
				in >> spe;
				if( !in )
				{
					cerr << "Seed File invalid.\n";
					return  false;
				}
				if(dx>=xTL && dx<xBR && dy>=yTL && dy<yBR)
				{
					ddx = dx-xTL;
					ddy = dy-yTL;
					data(ddx,ddy) = spe;
				}
			}
		break;
	}
	return true;
}



template <class Type> bool RWFile::WriteSeed(const char * fname,simplmat<Type>& data, const char * dataType)
{
	unsigned int i,j,dx,dy;
	bool privez=false;

	// Esta alreves simplmat= (row,col) y yo use siempre (col,row) ==> x=rows y=cols
	DimX = data.getRows();
	DimY = data.getCols();
	
	fstream sav;
	sav.open( fname, ios::in );
	if(!sav )
		privez=true;
	else
	{
		sav >> dx  >> dy;
		if( dx!=DimX || dy!=DimY )
		{
			cerr << "Appending to a file with different dimension" << endl;
    		return false;
		}
	}
	sav.close();
	sav.clear();

	sav.open( fname, ios::out | ios::app );
	if(!sav)
	{
		cerr << "Cannot open Seed file.\n";
		return 0;
	}

    string dType = dataType;
    if( dType.empty() )
    	dType = "BI";

	if(privez)
		sav << DimX << "\t" << DimY << endl;

	sav << dType << endl;
	for(i=0; i<DimY; i++)
	{
		for(j=0;j<DimX;j++)
		{
			sav<< data(j,i) << "\t";
		}
		sav << endl;
	}
	sav << endl;

	return true;
}


template <class Type> bool RWFile::WriteSeed(const char * fname,simplmat<Type>& data, const char * dataType,
			int xTL, int yTL, int xLen, int yLen)

{
	int i,j;
	bool privez=false;

	fstream sav;
	sav.open( fname, ios::in );
	if(!sav )
		privez=true;
	else
    {
		sav >> DimX  >> DimY;
		if( xLen!=DimX || yLen!=DimY )
		{
			cerr << "Appending to a file with different dimension" << endl;
    		return false;
		}
    }
	sav.close();
    
	sav.open( fname, ios::out | ios::app );
	if(!sav)
	{
		cerr << "Cannot open Seed file.\n";
		return false;
	}
	// Esta alreves simplmat= (row,col) y yo use siempre (col,row) ==> x=rows y=cols
	DimX = data.getRows();
	DimY = data.getCols();

    int xBR = xTL + xLen;
	int yBR = yTL + yLen;

	if( xTL<0 || yTL<0 || xBR>DimX  || yBR>DimY || xLen<0 || yLen<0)
	{
		cerr << "Invalid window" << endl;
		return false;
	}

    string dType = dataType;
    if( dType.empty() )
    	dType = "BI";

	if(privez)
		sav << xLen << "\t" << yLen << endl;
		
	sav << dType << endl;
	for(i=yTL; i<yBR; i++)
	{
		for(j=xTL;j<xBR;j++)
		{
			sav<< data(j,i) << "\t";
		}
		sav << endl;
	}
	sav << endl;

	return true;
};



template <class Type> bool RWFile::ReadIdrisi(const char * fname, simplmat<Type>& data)
{
	string dname,iname;
	iname = fname;

    string::size_type pos=0;

	if( (pos=iname.find(".rst")) == string::npos )
    {
		iname += ".rst";
		dname = fname;
		dname += ".rdc";
    }
	else
    {
    	dname = iname.substr(0,pos) + ".rdc";
    }
	
//	string buff;
	char buff[256],dType[10], *ptr;
	ifstream in;
	int dx=0,dy=0;

	in.open(dname.c_str());
	if( !in )
	{
		cerr << "Cannot open doc file." << endl;
		return 0;
	}

	while( !in.eof() )
	{
		in.getline(buff,255);
		
		ptr =strstr(buff,"columns");
		if( ptr!=NULL )
			dx = atoi(ptr+13);

		ptr=strstr(buff,"rows");
		if( ptr!=NULL )
			dy = atoi(ptr+13);

		ptr=strstr(buff,"data type");
		if( ptr!=NULL )
		{
			strncpy(dType,ptr+14,10);
			dType[9]='\0';
   		}
        ptr=strstr(buff,"file type");
		if( ptr!=NULL )
		{
	        ptr=strstr(buff,"ascii");
            if( ptr!=NULL )
            {
            	cerr << "ASCII File type not supported" << endl;
                return 0;
            }
		}

	}

	in.close();
	in.clear();

	if( dx!=data.getRows() || dy!=data.getCols() )
		data.resize(dx,dy);
		
	DimX=dx;
	DimY=dy;

	int spe=0;
	char echa=0;
//	float *eflo= data.pointer();
	float eflo=0;
   short int eint=0;
	
	in.open(iname.c_str(), ios::binary );
	if( !in )
	{
		cerr << "Cannot open rst file.\n";
		return 0;
	}

	if( strstr( dType,"byte")!=NULL)
	{
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
				in.read(&echa,1);
                if(echa<0)
                	spe = 256 + echa;
                else
					spe = echa;
				data(dx,dy) = spe;
			}
	}
	else if( strstr( dType,"integer")!=NULL )
	{
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
				in.read(reinterpret_cast<char *>(&eint),sizeof(short));
				spe = eint;
				data(dx,dy) = spe;
			}
	}
	else
	{
//		in.read(eflo,DimX*DimY);
		for(dy=0;dy<DimY; dy++)

			for(dx=0;dx<DimX; dx++)
			{
				in.read((char *)&eflo,sizeof(float));
				data(dx,dy) = eflo;
			}
	}

	return 1;
}


// Discretiza datos de precencia desde el formato x y a el formato de matriz.
//
// Parameters
// 		finp : Input File
//		xsize: x length of the total square 
//		xstep: Size of the step in x
//		ysize: y length of the total square 
//		ystep: Size of the step in y
//		data: ouput matrix
//		option: 0=Solamente coordenadas
//		        1=Coordenadas mas columna de datos
//				2=Coordenadas mas 2 columnas de datos: val val2 
//		Si hay coordenadas duplicadas debido a la discretizacion 
//      option=0 suma 1
//      option=1 suma el valor de la columna de datos
//      option=2 deja el val que corresponde al mayor val2
//      
//      OJO asume val val2 >= 0!!

template <class Type> bool RWFile::ReadMapXY(const char * finp, float xsize,float xstep,float ysize,float ystep, simplmat<Type> &data, int option)
{
	int x,y,xx,yy,numTies=0;

	ifstream in(finp);
	if(!in)
	{
		cerr << "Can not open " << finp << endl;
		return(0);
	}
	xx=static_cast<int>(xsize/xstep);
	yy=static_cast<int>(ysize/ystep);

	data.resize(xx,yy,0.0);
	simplmat<Type> ties;
	if(option==2) ties.resize(xx,yy,0.0);
	long numRecs=0;
	
   	float kx,ky;
	Type dat,dat2;
	while(!in.eof())
	{
		switch(option)
		{
		case 2:
			in >> kx >> ky >> dat >> dat2;
			break;
		case 1:
			in >> kx >> ky >> dat;
			break;
		case 0:
			in >> kx >> ky;
			dat = 1;
			break;
		default:
			kx=0; ky=0; dat=0;
		}
		
		if( kx<0 || ky<0 )
		{
			cerr << "Error negative coords " << kx << "\t" << ky << endl;
			exit(0);
		}
		
			
		if( in.eof() ) break;

		if( kx == xsize ) kx = xsize-1;
		if( ky == ysize ) ky = ysize-1;

		x = static_cast<int>(kx/xstep);
		y = static_cast<int>(ky/ystep);
		if( x>=xx || y>=yy )
		{
			cerr <<	"Coord " << x << ":" << kx << "\t" << y << ":" << ky << " out of defined range!" << endl;
			exit(0);
		}
		else
		{
			switch(option)
			{
			case 2:
				if(data(x,y)>0.0)
				{
					numTies++;
					if( dat2>ties(x,y) )
					{
						ties(x,y)=dat2;
						data(x,y)=dat;
					}
				}
				else
				{
					ties(x,y)=dat2;
					data(x,y)=dat;
				}
				break;
			default:
				if(data(x,y)>0.0)
					numTies++;
				data(x,y)+=dat;
			}			
			// Para calcular promedio sería (data * n + dataNew)/(n+1)
			numRecs++;
		}
	}
	cerr << "Number of ties: " << numTies << endl;
	cerr << "Number of recs: " << numRecs << endl;

	return(1);
}

// Build a lattice of dimension xsize,ysize from a vector {x,y,sp} or {x,y,sp,ba}, 
//
// where 	lc is the local scale to build the lattice in the units of x,y
//			sp is a species number
//       	ba is the basal area or importance of this species
// parameters
//			readOption = 0  reads x,y
// 			     	     1 	reads x,y,sp
//			        	 2 	reads x,y,sp,ba
//			outOption = 0  choose the most abundant species inside the site
//						anyNumber choose the species with that number 
//
// Example
//         xsize=1000,ysize=500,lc=1: the vector X,Y in file is measured from x=0 to 1000 and y=0 to 500
//		   xsize=1000,ysize=500,lc=10: the vector X,Y in file is measured from x=0 to 10000 and y=0 to 5000 and each site is 10x10
//         xsize=100,ysize=50,lc=10: the vector X,Y in file is measured from x=0 to 1000 and y=0 to 500 and each site is 10x10
//
template <class Type> bool RWFile::ClusterizeXY(const char * finp, int xsize,int ysize,float lc, simplmat<Type> &data, int readOption, int outOption)
{
   	struct lines{ double x; double y; int sp; double ba;};
   	double kx=0.0,ky=0.0,ba=0.0;
   	int sp=0;
   	vector<lines> li; //pag 48

	ifstream in(finp);
	if(!in)
	{
		cerr << "Can not open " << finp << endl;
		return(0);
	}
	while(!in.eof())
	{
		switch(readOption)
		{
		case 2:
			in >> kx >> ky >> sp >> ba;
			break;
		case 1:
			in >> kx >> ky >> sp;
			break;
		case 0:
			in >> kx >> ky;
			break;
		default:
			cerr << "Error invalid option: " << readOption << endl;
		}
		
		if( kx<0 || ky<0 || sp<0 || ba<0)
		{
			cerr << "Error negative values " << kx << "\t" << ky << endl;
			exit(0);
		} 
		if( in.eof() ) break;

		li.push_back({kx,ky,sp,ba});
	}


	data.resize(xsize,ysize,0.0);
    unsigned long sumSp=0,lostSp=0,sumPatch=0;
    
    typedef unordered_map<int,unsigned long> CounterMap;
    typedef unordered_map<int,unsigned long>::value_type  CounterMap_type;

	for(int i=1; i<=xsize; i++)
        for(int j=1; j<=ysize; j++) {

		    // this calculates the number of individuals for each species
		    //
		    CounterMap counts;
		    for (const auto& l : li)
		  		if (l.x<i && l.y<j && l.x>=i-1 && l.y>=j-1) {                        
		                sumSp++;
		                CounterMap::iterator ii(counts.find(l.sp));
		                if (ii != counts.end()){
		                    ii->second++;
		                } else {
		                    counts[l.sp] = 1;
		                    // numSp++; total number of species = counts.size()
		        		}
		        }

		//  Put the most abundant species into the lattice (lambda expression)
		//
            if( counts.size()>0)
            {
	   			switch(outOption)
				{
				case 0:
                    {
		            auto maxSp = max_element(counts.begin(), counts.end(), 
		                [](const CounterMap_type& p1, const CounterMap_type& p2) {return p1.second < p2.second; });

		            data(i-1,j-1) = maxSp->first;
                    lostSp += counts.size()-1;
                    sumPatch++;
                    //cout << i << "\t" << j << "\tSp: " << maxSp->first << " Lo:" << lostSp << " Pa:" << sumPatch << endl;
                    }
					break;

				default:
    				auto it = counts.find(outOption);
					if( it!=counts.end() )
			            data(i-1,j-1) = outOption;
                }
            }
        }
    
	cerr << "Number of species patches lost: " << lostSp << endl;
	cerr << "Total number of individuals: " << sumSp << endl;
	cerr << "Total number of patches: " << sumPatch << endl;
    

	return(true);
}


// Reads a text file with structure {x,y} or {x,y,value}
// only returns a matrix with {x,y}
//
template <class Type> bool RWFile::ReadXYVec(const char * finp, simplmat<Type> &data,int option)
{
	int x,y,xx,yy;

	ifstream in(finp);
	if(!in)
	{
		cerr << "Can not open " << finp << endl;
		return(0);
	}
   	float kx,ky;
   	long tot=0,t;
	Type dat;
	while(!in.eof())
	{
		switch(option)
		{
		case 1:
			in >> kx >> ky;
			dat = 1;
			break;
		case 0:
			in >> kx >> ky >> dat;
			break;
		default:
			kx=0; ky=0; dat=0;
		}
		tot++;		
		if( kx<0 || ky<0 )
		{
			cerr << "Error negative coords " << kx << "\t" << ky << endl;
			exit(0);
		}
		if( in.eof() ) break;
	}
	data.resize(tot,2);
	t=0;
	in.close();
	in.open(finp);
	
	while(!in.eof())
	{
		switch(option)
		{
		case 1:
			in >> kx >> ky;
			dat = 1;
			break;
		case 0:
			in >> kx >> ky >> dat;
			break;
		default:
			kx=0; ky=0; dat=0;
		}
		data(t,0)=kx;
		data(t,1)=ky;
		t++;
		if( in.eof() ) break;
	}
	if(t!=tot)
	{
		cerr << "Error t!=tot" << endl;
		exit(0);
	}
	
	return(1);
}

#ifdef __linux  // ReadTiff
template <class Type> bool RWFile::ReadTiff(const char * fname, simplmat<Type>& data)
{
	string iname;
	iname = fname;

    string::size_type pos=0;

	if( (pos=iname.find(".tif")) == string::npos )
    {
		iname += ".tif";
    }

	TIFF* in = TIFFOpen(iname.c_str(),"r");
	if (in == NULL)
	{
		cerr << "Cannot open " << iname.c_str() << "\n";
		return 0;
	}

	uint16 bitspersample = 1;
	uint16 samplesperpixel;
	uint16 config, photometric;
	uint32 dx=0,dy=0;

	TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &dx);
	TIFFGetField(in, TIFFTAG_IMAGELENGTH, &dy);
	TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	if (bitspersample != 8 && bitspersample != 16)
	{
		cerr << iname.c_str() <<  ": Image must have at least 8-bits/sample\n";
		return (0);
	}

	if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric) ||
				samplesperpixel>1 || (photometric!=PHOTOMETRIC_MINISWHITE
				&& photometric!=PHOTOMETRIC_MINISBLACK) )
	{
		cerr << iname.c_str() <<  ": Image must be monochromatic\n";
		return (0);
    }

	TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);
	if (config != PLANARCONFIG_CONTIG)
	{
		cerr << iname.c_str() <<  ": Can only handle contiguous data packing\n";
		return (0);
	}

	register unsigned char *inptr;
	register uint32 j, i;
	unsigned char *inputline = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(in));
	if (inputline == NULL)
	{
		cerr << iname.c_str() <<  ": No space for scanline buffer\n";
		return( 0 );
	}
	
	if( dx!=data.getRows() || dy!=data.getCols() )
		data.resize(dx,dy);
		
	DimX=dx;
	DimY=dy;

	for (i = 0; i < dy; i++)
	{
		if (TIFFReadScanline(in, inputline, i, 0) <= 0)
			break;
		inptr = inputline;
		for (j = dx; j-- > 0;)
		{
		   data(dx-j-1,i) = *inptr++;
		}
	}

	TIFFClose(in);
}
#else  // ReadTiff

template <class Type> bool RWFile::ReadTiff(const char * fname, simplmat<Type>& data)
{
	return 0;
}

#endif // ReadTiff


#endif  // RWFILE_H
