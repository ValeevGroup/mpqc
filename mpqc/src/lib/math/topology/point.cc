
#include <stdlib.h>
#include "point.h"
#include <math/scmat/matrix.h>
#include <util/keyval/keyval.h>

cart_point::cart_point()
{
}
cart_point::cart_point(const cart_point&p)
{
  r[0]=p[0];
  r[1]=p[1];
  r[2]=p[2];
}
cart_point::cart_point(const double* p)
{
  r[0]=p[0];
  r[1]=p[1];
  r[2]=p[2];
}
cart_point::~cart_point() {};
double& cart_point::operator[](int i) { return r[i]; };
const double& cart_point::operator[](int i) const { return r[i]; };
double& cart_point::x() { return r[0]; };
double& cart_point::y() { return r[1]; };
double& cart_point::z() { return r[2]; };
const double& cart_point::x() const { return r[0]; };
const double& cart_point::y() const { return r[1]; };
const double& cart_point::z() const { return r[2]; };


#define CLASSNAME Point
#define PARENTS virtual public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Point::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

// Constructors
Point::Point(int in_dim) {dim=in_dim; x=new double[dim];}

Point::Point(const double *y, int in_dim=3)
{
    dim = in_dim;
    if (dim > 0) {
        x=new double[dim];
        if (!memcpy (x,y,dim*sizeof(double)))
          fprintf(stderr,"Bad copy in Point cnstrctr\n");
      }
}
    
Point::Point(const Point &in_point)
{
    dim=in_point.dim;
    if (dim > 0) {
        x=new double[dim];
        if (!memcpy (x,in_point.x,dim*sizeof(double)))
          fprintf(stderr,"Bad copy in Point copy ctor\n");
      }
}
    
Point::Point(RefSCVector &in)
{
    dim=in.dim().n();
    if (dim > 0) {
        x=new double[dim];
        in.convert(x);
      }
}

Point::Point(double a)
{
  dim = 1;
  x = new double[1];
  x[0] = a;
}

Point::Point(double a,double b)
{
  dim = 2;
  x = new double[2];
  x[0] = a;
  x[1] = b;
}

Point::Point(double a,double b,double c)
{
  dim = 3;
  x = new double[3];
  x[0] = a;
  x[1] = b;
  x[2] = c;
}

Point::Point(KeyVal&k)
{
  dim = k.count();
  if (dim) x = new double[dim];
  for (int i=0; i<dim; i++) {
      x[i] = k.doublevalue(i);
    }
}

Point::Point(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(dim);
  if (dim > 0) s.get(x);
}

void
Point::save_data_state(StateOut& so)
{
  so.put(dim);
  if (dim > 0) so.put(x,dim);
}

// Destructor
Point::~Point(void) { if (dim > 0) delete[] x;};

// Copy to double *
double *Point::copy(void) const
{
  if (dim > 0) {
      double *x_out = new double[dim];
      if (!memcpy (x_out,x,dim*sizeof(double)))
        fprintf(stderr,"Bad copy in point copy function\n");
      return x_out;
    }
  return 0;
}

Point& Point::operator=(const Point&p)
{
  if (this != &p) {
      if (dim > 0) delete[] x;
      dim = p.dim;
      if (dim > 0) {
          x = new double[dim];
          for (int i=0; i<dim; i++) x[i] = p.x[i];
        }
    }
  return *this;
}

    // Set an element to a particular value
double &Point::operator[](int i)
{
    if (i < 0 || i >= dim)
    {
	fprintf(stderr," Dimension out of bounds in Point[]\n");
	exit(1);
    }
    return x[i];
}

    // Set an element to a particular value
const double &Point::operator[](int i) const
{
    if (i < 0 || i >= dim)
    {
	fprintf(stderr," Dimension out of bounds in Point[]\n");
	exit(1);
    }
    return x[i];
}

// Print out a Point
void Point::print(FILE *fp)
{
    int i;
    for (i=0;i<dim;i++)	fprintf(fp," %12.8f",x[i]);
    fprintf(fp,"\n");
}

void
Point::resize(int d)
{
  if (dim > 0) delete[] x;
  dim = d;
  if (dim > 0) x = new double[dim];
}

int Point::dimension() const
{
  return dim;
}

DescribedClass_REF_def(Point);
ARRAY_def(RefPoint);
SET_def(RefPoint);
ARRAYSET_def(RefPoint);
