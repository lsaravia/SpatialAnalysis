#ifndef RIPLEYK_H
#define RIPLEYK_H
#include "smattpl.h"
#include <math.h>
#include <stdlib.h>

const double Pi = M_PI;

#ifndef min
#define min(x,y)    (((x) < (y)) ?  (x) : (y))
#endif

void RipleyK(char * outfile, simplmat <float>& ftr, int numAnn, float widAnn, int treeDensity, int xDist, int yDist );
void focBound(int & focb1, int & focb2, int fx, int fy, int xDist, int yDist);

#endif

