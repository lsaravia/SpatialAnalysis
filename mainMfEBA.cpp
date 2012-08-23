#include "mf.h"
#include "RWFile.h"


int main(int argc, char * argv[])
{
	RWFile file;
	simplmat <double> data;	
	simplmat <double> q;	


	if( argc < 6)
	{
		cerr << "Usage: mfEBA inputFile qFile minBox maxBox deltaBox option" << endl
        	 << "       option N: Not normalize measure" << endl
        	 << "              S: Normalize measure" << endl
        	 << "              A: Normalize measure and save" << endl;
        
		exit(1);
	}
   	string fname = argv[1];
	
	if( fname.find(".rst") != string::npos )
    {
		if(!file.ReadIdrisi(argv[1], data))
			exit(1);
    }
	else if( fname.find(".tif")!=string::npos )
	{
		if(!file.ReadTiff(argv[1], data))
			exit(1);
	}
	else
    {
		if(!file.ReadSeed(argv[1], data))
			exit(1);
    }
   	
   
	if(!file.ReadSeed(argv[2], q))
		exit(1);

	int minBox = atoi(argv[3]);
	int maxBox = atoi(argv[4]);
	int deltaBox = atoi(argv[5]);
    
    char norm=argv[6][0];
    if(norm=='A') norm='S';
	
	MultifractalEBA(data, q,argv[1] ,minBox, maxBox, deltaBox, norm);

    if(argv[6][0]=='A')
    {
	    string::size_type pos=0;
		if( (pos=fname.find(".")) == string::npos )
	    {
			fname += "Norm.sed";
	    }
		else
		{
	    	fname = fname.substr(0,pos) + "Norm.sed";
		}

    	if(!file.WriteSeed(fname.c_str(),data))
			exit(1);
    }
}

