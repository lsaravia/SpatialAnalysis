#include "RipleyK.h"

void RipleyK(char * outfile, simplmat <float>& ftr, int numAnn, float widAnn, int treeDensity, int xDist, int yDist )
{
	simplmat <float> ripknl(100);
	simplmat <float> pairs(100);
	simplmat <float> kripley(100);
	simplmat <float> lripley(100);
	simplmat <float> ce95(100);
	simplmat <float> l_d(100);
	simplmat <float> pdf(100);
	simplmat <float> frank(100);
	simplmat <float> kfrank(100);
	simplmat <float> lfrank(100);

	// cout << "Ripley's K for run #" << run << "Tree density of " << treeDensity << endl;

	float searchDist = numAnn * widAnn;
	
	const int xCoor = 0;
	const int yCoor = 1;

	ripknl.fill(0.0);
	frank.fill(0.0);
	kfrank.fill(0.0);
	l_d.fill(0.0);
	ce95.fill(0.0);
	pdf.fill(0.0);

	float dist=0;
	int distpt=0;
	int foc,tar;
	for(foc=0; foc<treeDensity; foc++)
	{
		int fx = ftr(foc,xCoor);
		int fy = ftr(foc,yCoor);
                
        int focb1=0,focb2=0;
		// Determine distance from point FOC to nearest 2 boundaries
		//
		focBound(focb1,focb2,fx,fy,xDist,yDist);

		for(tar=0; tar<treeDensity; tar++)
		{
			if(tar!=foc)
			{
				int tx=ftr(tar,xCoor);
				int ty=ftr(tar,yCoor);
				int xDiff = abs(tx-fx);
				int yDiff = abs(ty-fy);
				
				if( xDiff <= searchDist && yDiff <= searchDist )
				{
					dist = sqrt(xDiff*xDiff + yDiff*yDiff );
					distpt = int(dist/widAnn) + 1;
					if(distpt <= numAnn)
					{
						double wij=0;
						if( dist<=focb1) // No adjustement
							wij=1;
						else
						{
							double arccos1=0,arccos2=0;
							double x1=0,x2=0;
							if( (dist*dist) <= (focb1*focb1+focb2*focb2) ) // 1-border adjustement
							{
								x1 = focb1/dist;
								arccos1 = Pi/2 - atan(x1/sqrt(1-x1*x1));
								wij = 1 / (1 - arccos1/Pi);
							}
							else
							{
								x1 = focb1/dist;
								x2 = focb2/dist;
								arccos1 = Pi/2 - atan(x1/sqrt(1-x1*x1));
								arccos2 = Pi/2 - atan(x2/sqrt(1-x2*x2));
								wij = 1 / (1 - (arccos1+arccos2+Pi/2)/(2*Pi));
							}
						}
						if( wij < 1 )
							cerr << "Error";

						if( wij > 4 )
							wij = 4;
							
						ripknl(distpt) += wij;
						double c1 = (4/3/Pi)*(dist/xDist + dist/yDist);
						double c2 = ((11/3/Pi) - 1 )*(dist*dist/(yDist*xDist));
						double fc1c2 = 1 - c1-c2;
						if( fc1c2 < 0 )
						{
							cerr << "Bogus fc1c2 of " << fc1c2 << endl;
							exit(1);
						}
						frank(distpt) += fc1c2;
					
					}
				}
			}
		} // next tar
	} // next foc

	double area = xDist*yDist;
	double areaDens = area/(treeDensity*(treeDensity-1));
	int ppp=0;
	
	for(distpt=0;distpt<numAnn; distpt++)
	{
		for(ppp=0; ppp<distpt; ppp++)
		{
			kripley(distpt) += ripknl(ppp);
			kfrank(distpt) += frank(ppp);
		}
		kripley(distpt) = areaDens * kripley(distpt);
		kfrank(distpt) = areaDens * kfrank(distpt);
	}
		
	int ddd;
	for(ddd=0; ddd<numAnn; ddd++)
	{
		lripley(ddd) = sqrt(kripley(ddd)/Pi);
		lfrank(ddd)  = sqrt(kfrank(ddd)/Pi);
		double d = ddd*widAnn;
		double c1 = (4/(3*Pi))*(d/xDist + d/yDist);
		double c2 = (4- 35/(3*Pi))*(d*d/area);
		double corr = 1-c1-c2;
		ce95(ddd) = 1.95*sqrt(areaDens/(2*Pi*corr));
		l_d(ddd) = lripley(ddd) - d;
		if( ddd>0 )
			pdf(ddd) = lripley(ddd) - lripley(ddd-1) - widAnn;
	}

	cout << " Tree density was" << treeDensity << endl;
	for(ddd=0; ddd<numAnn; ddd++)
	{
		// cout << head << ddd*widAnn << l_d(ddd) << endl;
		cout << ddd*widAnn << "\t" << lripley(ddd) << "\t" << l_d(ddd) << "\t" << ce95(ddd);
		if( fabs(l_d(ddd)) >= ce95(ddd))
			cout << "\t*";
		else
			cout << "\t";
			
		cout << "\t" << pdf(ddd) << endl;
	}
}



void focBound(int & focb1, int & focb2, int fx, int fy, int xDist, int yDist)
{
	// Nearest boundary distance
	focb1 = min(min(fx,xDist-fx), min(fy, yDist-fy));

	// Second closest boundary
	if( focb1 == fx || focb1 == xDist-fx)
		focb2 = min(fy, yDist-fy);
	else
		focb2 = min(fx, xDist-fx);
}
