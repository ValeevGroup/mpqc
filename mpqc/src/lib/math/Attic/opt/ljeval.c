#ifndef HEADER
#define HEADER 
#include <math.h>
#include <stdio.h>

/* Define DEBUG to print out more info than you want to see */
#define NDEBUG

/* This dimensions a single 1d array (ie this can be set very large) */
#define MAX_ATOMS 1000

/* define a 3d data type */
typedef struct {
  double x;
  double y;
  double z;
} coord_struct;

/* function to calculate internal energy given a list of cartesian coords */

void   angle_to_cart(double *, coord_struct *, int);

#endif

double calc_energy(coords,n_atoms)
     coord_struct coords[];
     int n_atoms;
{
  double energy=0.0,distance();
  int i,j;
  double r, rbar, width, y, f, rsub;	
  for (i=0;i<n_atoms; i++) 
    for (j=0;j<i;j++) {
      /* Here is where the energy is calced, so modify it at will */
      r     = distance(coords[i],coords[j]);
      rbar  = 0.05;
      width = 0.10;
      y     = (r-rbar) / width;
      f     = 10*y*y*y - 15.0 *y*y*y*y+6*y*y*y*y*y; 
      if(y<0) f = 0.0;
      if(y>1) f = 1.0;
      rsub  = r + (rbar+width)*(1.0-f);
      energy+=1./rsub;
    }
  return(energy);
}

/* function to calculate the distance between two 3d points */
double distance(coord1,coord2)
     coord_struct coord1,coord2;
{
  double value;
  value=sqrt((coord1.x-coord2.x)*(coord1.x-coord2.x)+
	     (coord1.y-coord2.y)*(coord1.y-coord2.y)+
	     (coord1.z-coord2.z)*(coord1.z-coord2.z));
  return(value);
}

/* This is a function that converts a list of
   angles into a list of cartesian coordinates */

void angle_to_cart(double angle[], coord_struct coords[], int n_atoms)
{

  int i;
  double angle_v1;
  coord_struct v1,v2;

  /* set up initial points and direction vectors */
  coords[0].x=coords[0].y=coords[0].z=0.0;
  coords[1].x=coords[1].z=0.0;
  coords[1].y = -1;

  /* v1 points from atom i-2 to atom i-1 */
  v1.x=0.0;
  v1.y= -1.0;
  v1.z=0.0;

  /* v2 is the new v1, but z must be set since it's never used */
  v2.z=0.0;

  /* now loop over the rest of the natoms, calculating their 
     cartesian positions */
  for (i=2; i<n_atoms; i++)
    {
      /* calculate next vector relative to the x-axis */
      v2.x=cos(angle[i-2]);
      v2.y=sin(angle[i-2]);

      /* rotate this vector to v1 origin */
      angle_v1= -atan2(-v1.y,-v1.x);

      v1.x=cos(angle_v1)*v2.x+sin(angle_v1)*v2.y;
      v1.y= -sin(angle_v1)*v2.x+cos(angle_v1)*v2.y;
      
      coords[i].x=coords[i-1].x+v1.x;
      coords[i].y=coords[i-1].y+v1.y;
      coords[i].z=coords[i-1].z+v1.z;

#ifdef DEBUG
      printf("vector %d: %lf %lf %lf\n",i,v1.x,v1.y,v1.z);
#endif
    }
}

#define MBIG 1000000000
#define MSEED 161803123
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(idum)
int *idum;
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
