
// some inline functions for dealing with 3 dimensional vectors

#ifndef _localdef_h
#define _localdef_h

#include <math.h>

static const double pi=3.14159265358979323846;
static const double pih=1.57079632679489661923;
static double tpi=2.0*pi;

static const double bohr = 0.52917706;

///////////////////////////////////////////////////////////

static inline void
delta(Point& u, Point& a, Point& b)
{
  u[0]=a[0]-b[0];
  u[1]=a[1]-b[1];
  u[2]=a[2]-b[2];
  }

static inline void
delta(double u[], const double a[], const double b[])
{
  u[0]=a[0]-b[0];
  u[1]=a[1]-b[1];
  u[2]=a[2]-b[2];
  }

///////////////////////////////////////////////////////////

// returns the distance between two points
static inline double
dist(Point& a, Point& b)
{
  double x,y,z;
  return (sqrt((x=a[0]-b[0])*x + (y=a[1]-b[1])*y + (z=a[2]-b[2])*z));
}

static inline double
dist(const double a[], const double b[])
{
  double x,y,z;
  return (sqrt((x=a[0]-b[0])*x + (y=a[1]-b[1])*y + (z=a[2]-b[2])*z));
}

///////////////////////////////////////////////////////////

// given sin(x) returns cos(x) 
static inline double
s2(double x)
{
  return (sqrt(1.0-x*x));
  }

///////////////////////////////////////////////////////////

// returns the dot product for two vectors
static inline double
scalar(Point& a, Point& b)
{
  double x = a[0]*b[0];
  double x1 = a[1]*b[1];
  x += a[2]*b[2];
  return x+x1;
  }

static inline double
scalar(const double a[], const double b[])
{
  double x = a[0]*b[0];
  double x1 = a[1]*b[1];
  x += a[2]*b[2];
  return x+x1;
  }

///////////////////////////////////////////////////////////

// given vectors a and b, returns a unit vector directed along the difference
// of the two vectors
static inline void
norm(Point& u, Point& a, Point& b)
{
  delta(u,a,b);
  double x = 1.0/sqrt(scalar(u,u));
  u[0] *= x; u[1] *= x; u[2] *= x;
  }

static inline void
norm(double u[], const double a[], const double b[])
{
  delta(u,a,b);
  double x = 1.0/sqrt(scalar(u,u));
  u[0] *= x; u[1] *= x; u[2] *= x;
  }

///////////////////////////////////////////////////////////

// given two vectors, returns the normalized cross product of those vectors
static inline void
normal(Point& a, Point& b, Point& w)
{
  w[0] = a[1]*b[2]-a[2]*b[1];
  w[1] = a[2]*b[0]-a[0]*b[2];
  w[2] = a[0]*b[1]-a[1]*b[0];
  double x = 1.0/sqrt(scalar(w,w));
  w[0] *= x; w[1] *= x; w[2] *= x;
  }

static inline void
normal(const double a[], const double b[], double w[])
{
  w[0] = a[1]*b[2]-a[2]*b[1];
  w[1] = a[2]*b[0]-a[0]*b[2];
  w[2] = a[0]*b[1]-a[1]*b[0];
  double x = 1.0/sqrt(scalar(w,w));
  w[0] *= x; w[1] *= x; w[2] *= x;
  }

#endif
