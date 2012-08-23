#ifndef RIPLEY_TRIANGLE_H
#define RIPLEY_TRIANGLE_H

#define ver "RipleyTriangle FG & RP 3/8/99"

#define pointmax 500
#define rmax 61
#define trimax 21
#define epsilon 0.0001

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

const float Pi=M_PI;



/*
 * Load point coordinates: (x,y) data file (2 columns)
 */
int charge_semis(char *nomfic, int *point_nb, float x[], float y[]);

/*
 * Load sampling window parameters: (x,y) data file (2 columns).
 * Rows #1 and #2: vertices of the rectangular window, (xmin,xmax) and (ymin,ymax).
 * From row #3: one triangle vertice per row, (ax,ay), (bx,by) and (cx,cy)
 */

int charge_triangles(char *nomfic, float *xmi, float *xma, float *ymi,
		float *yma,	int *triangle_nb, float ax[], float ay[], float bx[],
		float by[], float cx[], float cy[]);

/******************************************************************************/
/* Trigonometric routines	                                                  */
/******************************************************************************/

float bacos(float a);

/* return 1 if (x,y) is on the same side of the straight line (ab) than c */
int in_droite(float x, float y, float ax, float ay, float bx, float by,float cx, float cy);


/* return 1 if (x,y) is inside triangle abc or on its boundaries */
int in_triangle(float x, float y, float ax, float ay, float bx, float by,float cx, float cy);

/*Edge effect correction: give the perimeter/ddd of the circle centred on (xxx,yyy) with radius ddd which is
 *included within the rectangle defined by (xmi,xma,ymi,yma)*/
float perim_in_rect(float xxx, float yyy, float ddd, float xmi, float xma,float ymi, float yma);

/* a outside ; b and c inside*/
float un_point(	float ax, float ay, float bx, float by, float cx, float cy,
				float x, float y, float d);


/* a outside, b inside, c on the boundary*/
float ununun_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d);

/* a outside, b and c on the boundary*/
float deuxbord_point(float ax, float ay, float bx, float by, float cx, float cy,
						float x, float y, float d);

/* a inside, b and c outside*/
float deux_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d);

/* a on the boundary , b and c outside*/
float deuxun_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d);
					
/*a,b and c outside*/
float trois_point(float ax, float ay, float bx, float by, float cx, float cy,
					float x, float y, float d);

/*return the sum of angles included within the triangles*/
float perim_triangle(float x,float y, float d, int triangle_nb, float ax[], float ay[],
					float bx[], float by[], float cx[], float cy[]);


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
					float g[],float k[], float densite);


/******************************************************************************/
/* Routines for CI computation                                                */
/******************************************************************************/
void initalea(void);

/*return a long int between 0 and z*/
long aleaS(long z);

void s_alea_tr_rect(int point_nb,float x[], float y[], float xmi, float xma,
				float ymi, float yma, int triangle_nb, float ax[], float ay[],
				float bx[], float by[], float cx[], float cy[], float p);

void swap(float *tab,int a, int b);

void trirapide(float *tab,int g, int d);

#endif  // RIPLEY_TRIANGLE_H
