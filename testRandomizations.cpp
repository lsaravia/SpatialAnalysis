// TestRWFile
// 
#include "RWFile.h"
#include "Randomizations.h"

int main()
{
	RWFile file;
	simplmat <double> data;
	simplmat <double> rdata;

	file.ReadSeed("trnz.sed", data);
	Randomizations rz;

	rz.Randomize(data,rdata);

	int xd=rdata.getRows();
	int yd=rdata.getCols();

	int i,j;
	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%3.0f",rdata(i,j));
			
		cout << endl;
	}

	cout << endl;
	rz.Randomize(data);

	for(i=0; i<xd; i++)
	{
		for( j=0; j<yd; j++)
			cout.form("%3.0f",data(i,j));
			
		cout << endl;
	}

}
