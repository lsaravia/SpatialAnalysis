#ifndef RINGSTAT_H
#define RINGSTAT_H

#include "smattpl.h"
//#include "r250.h"
#include <time.h>
#include "randlib.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

class RingStatAnalysis
{
	public:

    RingAnalysis(int rSeed=0) {
    					if(rSeed==0)
                        	rSeed=time(0);
                       	setall(rSeed,rSeed+1);
                        };


    unsigned Rand(unsigned r) { return ignuin(0,r-1); };


};


#endif  // SPECTRAL_H


