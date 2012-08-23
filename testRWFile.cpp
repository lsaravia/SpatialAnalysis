// TestRWFile
// 
#include "RWFile.h"

int main()
{
	RWFile file;
/*	simplmat <float> data;
    simplmat <int> dataSP;
    simplmat <int> dataAG;
    simplmat <int> dataAJ;
    
    float * p;
	if( !file.ReadSeed("testSed.sed", dataSP,"SP") )
    	cerr << "Error SP" << endl;
	if( !file.ReadSeed("testSed.sed", dataAG,"AG") )
    	cerr << "Error AG" << endl;
	if( !file.ReadSeed("testSed.sed", dataAJ,"AJ") )
    	cerr << "Error AJ" << endl;
    if( !file.WriteSeed("testSP.sed", dataSP, "SP") )
    	cerr << "Error write SP" << endl;
    if( !file.WriteSeed("testSP.sed", dataAG, "AG", 0,0,20,20) );
    if( !file.WriteSeed("testSP.sed", dataAJ, "AJ", 0,0,20,20) );
        

*/
	simplmat <float> data;
	simplmat <int> idata;
	data.resize(128,64,0);
	idata.resize(128,64,0);
//	p = data.pointer();
    
    float peri=16;
	data.fill(0.0);
	idata.fill(0);
    
    for(int i=0;i<data.getRows();i++)
	   	for(int j=0; j<data.getCols(); j++)
        {
           peri = 100+i;
           idata(i,j)=peri;
        }


//    file.WriteIdrisi("testRwfile.rst", idata);

//	data.fill(0.0);

	file.ReadIdrisi("a03w.rst", data);
	file.WriteSeed("a03w.sed", data);


}
