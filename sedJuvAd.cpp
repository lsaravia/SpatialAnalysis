#include "smattpl.h"
#include "RWFile.h"

int main(int argc, char * argv[])
{
	if( argc<4)
	{
		cerr << "Usage: sedJuvAd input output type"  << endl;
		cerr << " type  : A/J" << endl;
		cerr << " Reads the AJ section and output only adults or juvenils" << endl;
		cerr << " adults are odd and juvenils are even" << endl;
      
		exit(1);
	}
	simplmat<int> data;
	RWFile file;
	string iname = argv[1];
    string::size_type pos;


	if( !file.ReadSeed(argv[1],data,"AJ") )
            return 1;

    int dX=data.getRows();
    int dY=data.getCols();
	int i,j;
    
	if( *argv[3]=='A' )
	{
		for(i=0; i<dY; i++)
			for( j=0;j<dX;j++)
				if( (data(j,i)%2) == 0 )
					data(j,i) = 0;
	}
	else
	{
		for(i=0; i<dY; i++)
			for( j=0;j<dX;j++)
				if( (data(j,i)%2) == 1 )
					data(j,i) = 0;
	}
	file.WriteSeed(argv[2], data, "BI");
}
