/*--------------------------------------------------------------------------*/
/*                                           						        */
/* Computation of Ripley's K-function		                                */
/* for a sampling window of any shape					                    */
/* included within a rectangular initial window                             */
/* F.G. &  R.P. 2/09/96 -- 3/8/99  --- MODIFICADO POR L.A.S. 20/03/2000     */
/*                                           						        */
/* OJO usa indices a partir de 1            						        */
/*                                           						        */
/*--------------------------------------------------------------------------*/

/******************************************************************************/
/*Global declarations	                                                      */
/******************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
//#include "r250.h"
#include "randlib.h"
#include "RipleyTriangle.h"
//#include "fortify.h"


/*****************************************************************************/
/* Loading routines  	     		                                         */
/*****************************************************************************/

/* Load point coordinates: (x,y) data file (2 columns)*/

int charge_semis(char *nomfic, int *point_nb, float x[], float y[])
{       FILE *fp;

        fp=fopen(nomfic,"r");
        if (fp==NULL)
        {       printf("!! error: this file doesn't exist");
                return -1;
        }
        else
        {
        	*point_nb=0;
                while (feof(fp)==0)
                {
                        *point_nb+=1;
                        fscanf(fp,"%f",&x[*point_nb]);
                        fscanf(fp,"%f",&y[*point_nb]);
                        if(feof(fp)!=0)
                        	*point_nb-=1;

                }
                fclose(fp);
                return 0;
		}
}

/* Load sampling window parameters: (x,y) data file (2 columns).
 * Rows #1 and #2: vertices of the rectangular window, (xmin,xmax) and (ymin,ymax).
 * From row #3: one triangle vertice per row, (ax,ay), (bx,by) and (cx,cy)          */

int charge_triangles(char *nomfic, float *xmi, float *xma, float *ymi,
		float *yma,	int *triangle_nb, float ax[], float ay[], float bx[],
		float by[], float cx[], float cy[])
{
	FILE *fp;

	fp=fopen(nomfic,"r");
	if (fp==NULL)
	{
		printf("!! error: this file doesn't exist");
		return -1;
	}
	else
	{
		fscanf(fp,"%f",xmi);
		fscanf(fp,"%f",xma);
		fscanf(fp,"%f",ymi);
		fscanf(fp,"%f",yma);
		*triangle_nb=0;
		while (feof(fp)==0)
		{
			*triangle_nb+=1;
			fscanf(fp,"%f",&ax[*triangle_nb]);
			fscanf(fp,"%f",&ay[*triangle_nb]);
			fscanf(fp,"%f",&bx[*triangle_nb]);
			fscanf(fp,"%f",&by[*triangle_nb]);
			fscanf(fp,"%f",&cx[*triangle_nb]);
			fscanf(fp,"%f",&cy[*triangle_nb]);
            if(feof(fp)!=0)
              	*triangle_nb-=1;
		}
		fclose(fp);
	}

    return 0;
}


/******************************************************************************/
/* Trigonometric routines	                                                  */
/******************************************************************************/

float bacos(float a)
{
	float b;

	if (a>=1)
		b=0;
	else if (a<=-1)
		b=Pi;
	else
		b=acos(a);

	return b;
}

/* return 1 if (x,y) is on the same side of the straight line (ab) than c */

int in_droite(float x, float y, float ax, float ay, float bx, float by,float cx, float cy)
{
	float vabx,vaby,vacx,vacy,vamx,vamy,pv1,pv2;

	vabx=bx-ax;
	vaby=by-ay;
	vacx=cx-ax;
	vacy=cy-ay;
	vamx=x-ax;
	vamy=y-ay;
	pv1=vabx*vacy-vaby*vacx;
	pv2=vabx*vamy-vaby*vamx;

	if (((pv1>0)&&(pv2>=0))||((pv1<0)&&(pv2<=0)))
		return 1;
	else
		return 0;
}

/* return 1 if (x,y) is inside triangle abc or on its boundaries */

int in_triangle(float x, float y, float ax, float ay, float bx, float by,float cx, float cy)
{	int res;

	res=0;
	if (in_droite(x, y, ax, ay, bx, by, cx, cy)==1)
		if (in_droite(x, y, bx, by, cx, cy, ax, ay)==1)
			if (in_droite(x, y, cx, cy, ax, ay, bx, by)==1)
				res=1;

	return res;
}

/*Edge effect correction: give the perimeter/ddd of the circle centred on (xxx,yyy) with radius ddd which is
 *included within the rectangle defined by (xmi,xma,ymi,yma)*/

float perim_in_rect(float xxx, float yyy, float ddd, float xmi, float xma,float ymi, float yma)
{
		float p;
        float d1,d2,d3,d4;

        if ((xxx>=xmi+ddd)&&(yyy>=ymi+ddd)&&(xxx<=xma-ddd)&&(yyy<=yma-ddd))
                return 2*Pi;
        else
        {       d1=(xxx-xmi)/ddd;
                d2=(yyy-ymi)/ddd;
                d3=(xma-xxx)/ddd;
                d4=(yma-yyy)/ddd;
                if (d1>=1)
                {       if (d2>=1)
                        {       if (d3>=1)
/*contact with 1 boundary in d4*/
                                        return (2*(Pi-acos(d4)));
                                else
                                {       if (d4>=1)
/*contact with 1 boundary in d3*/
                                                return (2*(Pi-acos(d3)));
                                        else
/*contact with 2 boundaries in d3,d4*/
                                        {       if (d3*d3+d4*d4<1)
                                                        return(1.5*Pi-acos(d3)-acos(d4));
                                                else
                                                        return(2*(Pi-acos(d3)-acos(d4)));
                                        }
                                }
                        }
                        else
                        {       if (d3>=1)
                                {       if (d4>=1)
/*contact with 1 boundary in d2*/
                                                return (2*(Pi-acos(d2)));
                                        else
/*contact with 2 boudaries in d2,d4*/
                                                return(2*(Pi-acos(d2)-acos(d4)));
                                }
                                else
                                {       if (d4>=1)
/*contact with 2 boundaries in d2,d3*/
                                        {       if (d2*d2+d3*d3<1)
                                                        return((1.5*Pi-acos(d2)-acos(d3)));
                                                else
                                                        return(2*(Pi-acos(d2)-acos(d3)));
                                        }
                                        else
/*contact with 3 boundaries in d2,d3,d4*/
                                        {       if (d2*d2+d3*d3<1)
                                                {       if (d3*d3+d4*d4<1)
                                                                return((Pi-acos(d2)-acos(d4)));
                                                        else
                                                                return((1.5*Pi-acos(d2)-acos(d3)-2*acos(d4)));
                                                }
                                                else
                                                {       if (d3*d3+d4*d4<1)
                                                                return((1.5*Pi-2*acos(d2)-acos(d3)-acos(d4)));
                                                        else
                                                                return(2*(Pi-acos(d2)-acos(d3)-acos(d4)));
                                                }
                                        }
                                }
                        }
                }
                else
                {       if (d2>=1)
                        {       if (d3>=1)
                                {       if (d4>=1)
/*contact with 1 boundary in d1*/
                                                return (2*(Pi-acos(d1)));
                                        else
/*contact with 2 boundaries in d1,d4*/
                                        {       if (d1*d1+d4*d4<1)
                                                        return((1.5*Pi-acos(d1)-acos(d4)));
                                                else
                                                        return(2*(Pi-acos(d1)-acos(d4)));
                                        }
                                }
                                else
                                {       if (d4>=1)
/*contact with 2 boundaries in d1,d3*/
                                                return(2*(Pi-acos(d1)-acos(d3)));
                                        else
/*contact with 3 boundaries in d1,d3,d4*/
                                        {       if (d3*d3+d4*d4<1)
                                                {       if (d4*d4+d1*d1<1)
                                                                return((Pi-acos(d3)-acos(d1)));
                                                        else
                                                                return((1.5*Pi-acos(d3)-acos(d4)-2*acos(d1)));
                                                }
                                                else
                                                {       if (d4*d4+d1*d1<1)
                                                                return((1.5*Pi-2*acos(d3)-acos(d4)-acos(d1)));
                                                        else
                                                                return(2*(Pi-acos(d3)-acos(d4)-acos(d1)));
                                                }
                                        }
                                }
                        }
                        else
                        {       if (d3>=1)
                                {       if (d4>=1)
/*contact with 2 boundaries in d1,d2*/
                                        {       if (d1*d1+d2*d2<1)
                                                        return((1.5*Pi-acos(d1)-acos(d2)));
                                                else
                                                        return(2*(Pi-acos(d1)-acos(d2)));
                                        }
                                        else
/*contact with 3 boudaries in d1,d2,d4*/
                                        {       if (d4*d4+d1*d1<1)
                                                {       if (d1*d1+d2*d2<1)
                                                                return((Pi-acos(d4)-acos(d2)));
                                                        else
                                                                return((1.5*Pi-acos(d4)-acos(d1)-2*acos(d2)));
                                                }
                                                else
                                                {       if (d1*d1+d2*d2<1)
                                                                return((1.5*Pi-2*acos(d4)-acos(d1)-acos(d2)));
                                                        else
                                                                return(2*(Pi-acos(d4)-acos(d1)-acos(d2)));
                                                }
                                        }
                                }
                                else

                                {       if (d4>=1)
/*contact with 3 boundaries in d1,d2,d3*/
                                        {       if (d1*d1+d2*d2<1)
                                                {       if (d2*d2+d3*d3<1)
                                                                return((Pi-acos(d1)-acos(d3)));
                                                        else
                                                                return((1.5*Pi-acos(d1)-acos(d2)-2*acos(d3)));
                                                }
                                                else
                                                {       if (d2*d2+d3*d3<1)
                                                                return((1.5*Pi-2*acos(d1)-acos(d2)-acos(d3)));
                                                        else
                                                                return(2*(Pi-acos(d1)-acos(d2)-acos(d3)));
                                                }
                                        }
                                        else
                                        {       printf("erreur");
                                                return 0;
                                        }
                                }
                        }
                }
        }
}

/* a outside ; b and c inside*/

float un_point(	float ax, float ay, float bx, float by, float cx, float cy,
				float x, float y, float d)
{
	float alpha, beta, gamma, delta, ttt, ang;
	float ex,ey,fx,fy;

	/*first intersection point*/

	alpha=(bx-ax)*(bx-ax)+(by-ay)*(by-ay);
	beta=(2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay));
	gamma=((ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	if (delta<=0)
		printf("erreur1");
	ttt=(-beta-sqrt(delta))/(2*alpha);
	if ((ttt<=0)||(ttt>=1))
		printf("erreur2");
	ex=ax+ttt*(bx-ax);
	ey=ay+ttt*(by-ay);

	/*second intersection point*/

	alpha=(cx-ax)*(cx-ax)+(cy-ay)*(cy-ay);
	beta=(2*(ax-x)*(cx-ax)+2*(ay-y)*(cy-ay));
	delta=beta*beta-4*alpha*gamma;
	if (delta<=0)
		printf("erreur3");
	ttt=(-beta-sqrt(delta))/(2*alpha);
	if ((ttt<=0)||(ttt>=1))
		printf("erreur4");
	fx=ax+ttt*(cx-ax);
	fy=ay+ttt*(cy-ay);

	/*angle computation*/

	ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
	return ang;
}


/* a outside, b inside, c on the boundary*/

float ununun_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d)
{
	float alpha, beta, gamma, delta, ttt, ang;
	float ex,ey,fx,fy;

	/*first intersection point ab*/

	alpha=(bx-ax)*(bx-ax)+(by-ay)*(by-ay);
	beta=(2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay));
	gamma=((ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	if (delta<=0)
		printf("erreur1b");
	ttt=(-beta-sqrt(delta))/(2*alpha);
	if ((ttt<=0)||(ttt>=1))
		printf("erreur2b");
	ex=ax+ttt*(bx-ax);
	ey=ay+ttt*(by-ay);

	/*second intersection point ac*/

	alpha=(cx-ax)*(cx-ax)+(cy-ay)*(cy-ay);
	beta=(2*(ax-x)*(cx-ax)+2*(ay-y)*(cy-ay));
	delta=beta*beta-4*alpha*gamma;
	ttt=1;
	if (delta>0)
	{
		ttt=(-beta-sqrt(delta))/(2*alpha);
		if ((ttt<=0)||(ttt>1))
			ttt=1;
		if (ttt<=0)
			printf("e3b");
	}
	fx=ax+ttt*(cx-ax);
	fy=ay+ttt*(cy-ay);

	/*angle computation*/

	ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
	return ang;
}


/* a outside, b and c on the boundary*/

float deuxbord_point(float ax, float ay, float bx, float by, float cx, float cy,
						float x, float y, float d)
{
	float alpha, beta, gamma, delta, te,tf, ang;
	float ex,ey,fx,fy;

	/*first segment ab*/

	alpha=(bx-ax)*(bx-ax)+(by-ay)*(by-ay);
	beta=(2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay));
	gamma=((ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	te=1;
	if (delta>0)
	{
		te=(-beta-sqrt(delta))/(2*alpha);
		if ((te<=0)||(te>=1))
			te=1;
		if (te<=0)
			printf("e1t ");
	}
	ex=ax+te*(bx-ax);
	ey=ay+te*(by-ay);

	/*second segment ac*/

	alpha=(cx-ax)*(cx-ax)+(cy-ay)*(cy-ay);
	beta=(2*(ax-x)*(cx-ax)+2*(ay-y)*(cy-ay));
	delta=beta*beta-4*alpha*gamma;
	tf=1;
	if (delta>0)
	{
		tf=(-beta-sqrt(delta))/(2*alpha);
		if ((tf<=0)||(tf>=1))
			tf=1;
		if (tf<=0)
			printf("e4t ");
	}
	fx=ax+tf*(cx-ax);
	fy=ay+tf*(cy-ay);

	/*angle computation*/

	ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
	return ang;
}


/* a inside, b and c outside*/

float deux_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d)
{
	float alpha, beta, gamma, delta, ttt, ang;
	float ex,ey,fx,fy,gx,gy,hx,hy;
	int cas;

	/*first intersection point*/

	alpha=((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
	beta=(2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay));
	gamma=((ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	if (delta<=0)
		printf("erreur6");
	ttt=(-beta+sqrt(delta))/(2*alpha);
	if ((ttt<=0)||(ttt>=1))
		printf("erreur7");
	ex=ax+ttt*(bx-ax);
	ey=ay+ttt*(by-ay);

	/*second intersection point*/

	alpha=((cx-ax)*(cx-ax)+(cy-ay)*(cy-ay));
	beta=(2*(ax-x)*(cx-ax)+2*(ay-y)*(cy-ay));
	delta=beta*beta-4*alpha*gamma;
	if (delta<=0)
		printf("erreur8");
	ttt=(-beta+sqrt(delta))/(2*alpha);
	if ((ttt<=0)||(ttt>=1))
		printf("erreur9");
	fx=ax+ttt*(cx-ax);
	fy=ay+ttt*(cy-ay);

	/*Are there two other intersection points ?*/

	cas=0;
	alpha=((cx-bx)*(cx-bx)+(cy-by)*(cy-by));
	beta=(2*(bx-x)*(cx-bx)+2*(by-y)*(cy-by));
	gamma=((bx-x)*(bx-x)+(by-y)*(by-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	if (delta>0)
	{
		ttt=(-beta-sqrt(delta))/(2*alpha);
		if ((ttt>=0)&&(ttt<=1))
		{
			gx=bx+ttt*(cx-bx);
			gy=by+ttt*(cy-by);
			ttt=(-beta+sqrt(delta))/(2*alpha);
			if ((ttt>=0)&&(ttt<=1))
			{
				cas=1;
				hx=bx+ttt*(cx-bx);
				hy=by+ttt*(cy-by);
			}
			else
				printf("erreur9bis");
		}
	}

	/*angle computation*/

	if (cas==0)
	{
		ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
	}
	else
	{
		ang=bacos(((ex-x)*(gx-x)+(ey-y)*(gy-y))/(d*d));
		ang+=bacos(((fx-x)*(hx-x)+(fy-y)*(hy-y))/(d*d));
	}

	return ang;
}


/* a on the boundary , b and c outside*/

float deuxun_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d)
{
	float alpha, beta, gamma, delta, te,tf,tg,th, ang;
	float ex,ey,fx,fy,gx,gy,hx,hy;
	int cas;

	/*first intersection point*/

	alpha=((bx-ax)*(bx-ax)+(by-ay)*(by-ay));
	beta=(2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay));
	gamma=((ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	te=0;
	if (delta>0)
	{
		te=(-beta+sqrt(delta))/(2*alpha);
		if ((te<0)||(te>=1))
			te=0;
		if (te>=1)
			printf("e15 ");
	}
	ex=ax+te*(bx-ax);
	ey=ay+te*(by-ay);

	/*second intersection point*/

	alpha=((cx-ax)*(cx-ax)+(cy-ay)*(cy-ay));
	beta=(2*(ax-x)*(cx-ax)+2*(ay-y)*(cy-ay));
	delta=beta*beta-4*alpha*gamma;
	tf=0;
	if (delta>0)
	{
		tf=(-beta+sqrt(delta))/(2*alpha);
		if ((tf<0)||(tf>=1))
			tf=0;
		if (tf>=1)
			printf("e15 ");
	}
	fx=ax+tf*(cx-ax);
	fy=ay+tf*(cy-ay);

	/*Are there two other intersection points ?*/

	cas=0;
	alpha=((cx-bx)*(cx-bx)+(cy-by)*(cy-by));
	beta=(2*(bx-x)*(cx-bx)+2*(by-y)*(cy-by));
	gamma=((bx-x)*(bx-x)+(by-y)*(by-y)-d*d);
	delta=beta*beta-4*alpha*gamma;
	if (delta>0)
	{
		tg=(-beta-sqrt(delta))/(2*alpha);
		if ((tg>=0)&&(tg<=1))
		{
			gx=bx+tg*(cx-bx);
			gy=by+tg*(cy-by);
			th=(-beta+sqrt(delta))/(2*alpha);
			if ((th>=0)&&(th<=1))
			{
				cas=1;
				hx=bx+th*(cx-bx);
				hy=by+th*(cy-by);
			}
			else
				printf("erreur9ter");
		}
	}

	/*angle computation*/

	if (cas==0)
	{
		if ((te==0)&&(tf==0))
			ang=0;
		else
			ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
	}
	else
	{
		ang=bacos(((ex-x)*(gx-x)+(ey-y)*(gy-y))/(d*d));
		ang+=bacos(((fx-x)*(hx-x)+(fy-y)*(hy-y))/(d*d));
	}

	return ang;
}


/*a,b and c outside*/

float trois_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d)
{
	float alpha, beta, gamma, delta, te,tf,tg,th,ti,tj, ang;
	float ex,ey,fx,fy,gx,gy,hx,hy,ix,iy,jx,jy;
	int cas;

	/*first segment ab*/

	alpha=(bx-ax)*(bx-ax)+(by-ay)*(by-ay);
	beta=2*(ax-x)*(bx-ax)+2*(ay-y)*(by-ay);
	gamma=(ax-x)*(ax-x)+(ay-y)*(ay-y)-d*d;
	delta=beta*beta-4*alpha*gamma;
	if (delta<0)
	{
		te=-1;
		tf=-1;
	}
	else
	{
		te=(-beta-sqrt(delta))/(2*alpha);
		tf=(-beta+sqrt(delta))/(2*alpha);
		if ((te<0)||(te>=1)||(tf==0))
		{
			te=-1;
			tf=-1;
		}
		else
		{
			ex=ax+te*(bx-ax);
			ey=ay+te*(by-ay);
			fx=ax+tf*(bx-ax);
			fy=ay+tf*(by-ay);
			if ((tf<=0)||(tf>1))
				printf("\npb te %f tf %f",te,tf);
		}
	}

	/*second segment bc*/

	alpha=(cx-bx)*(cx-bx)+(cy-by)*(cy-by);
	beta=2*(bx-x)*(cx-bx)+2*(by-y)*(cy-by);
	gamma=(bx-x)*(bx-x)+(by-y)*(by-y)-d*d;
	delta=beta*beta-4*alpha*gamma;
	if (delta<0)
	{
		tg=-1;
		th=-1;
	}
	else
	{
		tg=(-beta-sqrt(delta))/(2*alpha);
		th=(-beta+sqrt(delta))/(2*alpha);
		if ((tg<0)||(tg>=1)||(th==0))
		{
			tg=-1;
			th=-1;
		}
		else
		{
			gx=bx+tg*(cx-bx);
			gy=by+tg*(cy-by);
			hx=bx+th*(cx-bx);
			hy=by+th*(cy-by);
			if ((th<=0)||(th>1))
				printf("\npb tg %f th %f",tg,th);
		}
	}

	/*third segment ca*/

	alpha=(ax-cx)*(ax-cx)+(ay-cy)*(ay-cy);
	beta=2*(cx-x)*(ax-cx)+2*(cy-y)*(ay-cy);
	gamma=(cx-x)*(cx-x)+(cy-y)*(cy-y)-d*d;
	delta=beta*beta-4*alpha*gamma;
	if (delta<0)
	{
		ti=-1;
		tj=-1;
	}
	else
	{
		ti=(-beta-sqrt(delta))/(2*alpha);
		tj=(-beta+sqrt(delta))/(2*alpha);
		if ((ti<0)||(ti>=1)||(tj==0))
		{
			ti=-1;
			tj=-1;
		}
		else
		{
			ix=cx+ti*(ax-cx);
			iy=cy+ti*(ay-cy);
			jx=cx+tj*(ax-cx);
			jy=cy+tj*(ay-cy);
			if ((tj<=0)||(tj>1))
				printf("\npb ti %f tj %f",ti,tj);
		}
	}

	/*which case ?*/

	if (te<0)
	{
		if (tg<0)
		{
			if (ti<0)
			{	/*no intersection*/
				ang=0;
			}
			else
			{	/*segment ca cut the circle in i,j*/
				ang=bacos(((ix-x)*(jx-x)+(iy-y)*(jy-y))/(d*d));
			}
		}
		else
		{	if (ti<0)
			{	/*segment bc cut the circle in g,h*/
				ang=bacos(((gx-x)*(hx-x)+(gy-y)*(hy-y))/(d*d));
			}
			else
			{
				/*segments bc and ca cut the circle in g,h,i,j*/
				ang=bacos(((gx-x)*(jx-x)+(gy-y)*(jy-y))/(d*d));
				ang+=bacos(((hx-x)*(ix-x)+(hy-y)*(iy-y))/(d*d));
			}
		}
	}
	else
	{
		if (tg<0)
		{
			if (ti<0)
			{
				/*segment ab cut the circle in e,f*/
				ang=bacos(((ex-x)*(fx-x)+(ey-y)*(fy-y))/(d*d));
			}
			else
			{	/*segments ab and ca cut the circle in e,f,i,j*/
				ang=bacos(((ex-x)*(jx-x)+(ey-y)*(jy-y))/(d*d));
				ang+=bacos(((fx-x)*(ix-x)+(fy-y)*(iy-y))/(d*d));
			}
		}
		else
		{	if (ti<0)
			{
				/*segments ab and bc cut the circle in e,f,g,h*/
				ang=bacos(((ex-x)*(hx-x)+(ey-y)*(hy-y))/(d*d));
				ang+=bacos(((fx-x)*(gx-x)+(fy-y)*(gy-y))/(d*d));
			}
			else
			{
				/*three segments cut the circle*/
				ang=bacos(((ex-x)*(jx-x)+(ey-y)*(jy-y))/(d*d));
				ang+=bacos(((hx-x)*(ix-x)+(hy-y)*(iy-y))/(d*d));
				ang+=bacos(((fx-x)*(gx-x)+(fy-y)*(gy-y))/(d*d));
			}
		}
	}
	if ((ang<0)||(ang>Pi))
		printf("\nerreur12 : ang=%3.2f, %d %d %d %d %d %d",ang,te,tf,tg,th,ti,tj);

	return ang;
}


/*return the sum of angles included within the triangles*/

float perim_triangle(float x,float y, float d, int triangle_nb, float ax[], float ay[],
					float bx[], float by[], float cx[], float cy[])
{
	float angle;
	float doa,dob,doc;
	int h,i;

	angle=0;
	for(h=1;h<=triangle_nb;h++)
	{
		doa=sqrt((x-ax[h])*(x-ax[h])+(y-ay[h])*(y-ay[h]));
		dob=sqrt((x-bx[h])*(x-bx[h])+(y-by[h])*(y-by[h]));
		doc=sqrt((x-cx[h])*(x-cx[h])+(y-cy[h])*(y-cy[h]));
		if (doa-d<-epsilon)
		{
			if (dob-d<-epsilon)
			{
				if (doc-d<-epsilon)
					i=1;		 /*triangle inside the circle, TVB*/
				else if (doc-d>epsilon)
					angle+=un_point(cx[h],cy[h],ax[h],ay[h],bx[h],by[h],x,y,d);
				else
					i=1;		 /*triangle inside the circle, TVB*/
			}
			else if (dob-d>epsilon)
			{
				if (doc-d<-epsilon)
				angle+=un_point(bx[h],by[h],ax[h],ay[h],cx[h],cy[h],x,y,d);
				else if (doc-d>epsilon)
					angle+=deux_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
				else
					angle+=ununun_point(bx[h],by[h],ax[h],ay[h],cx[h],cy[h],x,y,d);
			}
			else /* b on the boundary*/
			{
				if (doc-d<-epsilon)
					i=1;		 /*triangle inside the circle, TVB*/
				else if (doc-d>epsilon)
					angle+=ununun_point(cx[h],cy[h],ax[h],ay[h],bx[h],by[h],x,y,d);
				else
					i=1;		 /*triangle inside the circle, TVB*/
			}
		}
		else if (doa-d>epsilon)
		{
			if (dob-d<-epsilon)
			{
				if (doc-d<-epsilon)
					angle+=un_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
				else if (doc-d>epsilon)
					angle+=deux_point(bx[h],by[h],ax[h],ay[h],cx[h],cy[h],x,y,d);
				else
					angle+=ununun_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
			}
			else if (dob-d>epsilon)
			{
				if (doc-d<-epsilon)
					angle+=deux_point(cx[h],cy[h],ax[h],ay[h],bx[h],by[h],x,y,d);
				else if (doc-d>epsilon)
					angle+=trois_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
				else
					angle+=deuxun_point(cx[h],cy[h],ax[h],ay[h],bx[h],by[h],x,y,d);
			}
			else /* b on the boundary*/
			{
				if (doc-d<-epsilon)
					angle+=ununun_point(ax[h],ay[h],cx[h],cy[h],bx[h],by[h],x,y,d);
				else if (doc-d>epsilon)
					angle+=deuxun_point(bx[h],by[h],ax[h],ay[h],cx[h],cy[h],x,y,d);
				else
					angle+=deuxbord_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
			}
		}
		else /* a on the boundary*/
		{
			if (dob-d<-epsilon)
			{
				if (doc-d<-epsilon)
					i=1;		 /*triangle inside the circle, TVB*/
				else if (doc-d>epsilon)
					angle+=ununun_point(cx[h],cy[h],bx[h],by[h],ax[h],ay[h],x,y,d);
				else
					i=1;		 /*triangle inside the circle, TVB*/
			}
			else if (dob-d>epsilon)
			{
				if (doc-d<-epsilon)
					angle+=ununun_point(bx[h],by[h],cx[h],cy[h],ax[h],ay[h],x,y,d);
				else if (doc-d>epsilon)
					angle+=deuxun_point(ax[h],ay[h],bx[h],by[h],cx[h],cy[h],x,y,d);
				else
					angle+=deuxbord_point(bx[h],by[h],ax[h],ay[h],cx[h],cy[h],x,y,d);
			}
			else /*b on the boundary*/
			{
				if (doc-d<-epsilon)
					i=1;		 /*triangle inside the circle, TVB*/
				else if (doc-d>epsilon)
					angle+=deuxbord_point(cx[h],cy[h],ax[h],ay[h],bx[h],by[h],x,y,d);
				else
					i=1;		 /*triangle inside the circle, TVB*/
			}
		}
	}
	return angle;
}

/**********************************************************************************/
/* This routine computes Ripley's K-function for a point pattern given in tables  */
/* x and y. The points are included within the rectangular area (xmi,xma,ymi,yma) */
/* from which some triangles are excluded.						                  */
/* The edge effects are corrected following Ripley's method,i.e. using the inverse*/
/* rectangualr window and outside the triangles.				                  */
/* The computation is made for t2 distance intervals of length dt.                */
/* The routine computes also the pair density function g,the local density		  */
/* function n, and the linearized L-function.                              	      */
/**********************************************************************************/

void ripley_tr_rect(int point_nb, float x[],float y[],float xmi,float xma,
					float ymi,float yma, int triangle_nb,float ax[],float ay[],
					float bx[],float by[],float cx[],float cy[], int t2, float dt,
					float g[],float k[], float densite)
{
	int i,j,tt;
	float d,cin;

	/*compute the number of point couples per distance g*/
	for(tt=1;tt<=t2;tt=tt+1)
		g[tt]=0;
	for(i=2;i<=point_nb;i=i+1)
		for(j=1;j<i;j=j+1)
		{
			d=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
			if (d<t2*dt)
			{
				tt=d/dt+1;
				/*for point i*/
				cin=perim_in_rect(x[i],y[i],d,xmi,xma,ymi,yma);

				cin=cin-perim_triangle(x[i],y[i],d,triangle_nb,ax,ay,bx,by,cx,cy);
				if (cin<0)
					printf("\ncin1<0");
				g[tt]+=2*Pi/cin;
				/*for point j*/
				cin=perim_in_rect(x[j],y[j],d,xmi,xma,ymi,yma);

				cin=cin-perim_triangle(x[j],y[j],d,triangle_nb,ax,ay,bx,by,cx,cy);
				if (cin<0)
					printf("\ncin2<0");
				g[tt]+=2*Pi/cin;
			}
		}

	for(tt=1;tt<=t2;tt=tt+1)
		g[tt]=g[tt]/point_nb;	/*mean -> density*/

	k[1]=g[1];
	for(tt=2;tt<=t2;tt=tt+1)	/*integration*/
		k[tt]=k[tt-1]+g[tt];
}

/******************************************************************************/
/* Routines for CI computation                                                */
/******************************************************************************/

void initalea(void)
{
	int ii,cro;
    double chi;

    ii=1;
    chi=difftime(time(0),0);
    while (chi>INT_MAX)
    {
    	chi=chi/INT_MAX;
        ii=ii+1;
    }
    while (ii>1)
    {
    	cro=chi;
        chi=chi-cro;
        chi=INT_MAX*chi;
        ii=ii-1;
    }
    cro=chi;
    srand(cro);
//    r250_init(cro);
	setall(rand(),rand());
    
}


/*return a long int between 0 and z*/
long aleaS(long z)	/* v4 */
{
/*
	int r1,r2;
	double f;
	unsigned long r;
	if (((RAND_MAX+1)*(RAND_MAX+1)-1)>LONG_MAX)
	{	printf("ea4");
		return -1;
	}
	else
	{	r=(long) rand()*float(RAND_MAX+1)+rand();
		f=(double) r/(RAND_MAX+1);
      	f=(double) f/(RAND_MAX+1);
		f=(double) f*(z+1);
      	r=(long) f;
		return r;
	}
*/
//	return r250n(z+1);
//	double Rand() { return ranf(); };
   	return ignuin(0,z);
    
}

void s_alea_tr_rect(int point_nb,float x[], float y[], float xmi, float xma,
				float ymi, float yma, int triangle_nb, float ax[], float ay[],
				float bx[], float by[], float cx[], float cy[], float p)
{	int i,j;

	for(i=1;i<=point_nb;i=i+1)
	{	x[i]=xmi+aleaS((xma-xmi)/p)*p;
		y[i]=ymi+aleaS((yma-ymi)/p)*p;
		for (j=1;j<=triangle_nb;j=j+1)
		if (in_triangle(x[i],y[i],ax[j],ay[j],bx[j],by[j],cx[j],cy[j]))
		{	j=triangle_nb+1;
			i=i-1;
		}
	}
}

void swap(float *tab,int a, int b)
{       float f;

        f=tab[a];
        tab[a]=tab[b];
        tab[b]=f;
}

void trirapide(float *tab,int g, int d)
{       int i,de;

        if (g>=d)
                return;
        swap(tab,g,(g+d)/2);
        de=g;
        for (i=g+1;i<=d;i++)
                if (tab[i]<tab[g])
                        swap(tab,++de,i);
        swap(tab,g,de);
        trirapide(tab,g,de-1);
        trirapide(tab,de+1,d);
}



