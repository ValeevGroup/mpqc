
extern "C" {
# include <stdio.h>
# include <math.h>
};

#include "volume.h"

#define CLASSNAME Volume
#define PARENTS public NLP2
#include <util/state/statei.h>
#include <util/class/classia.h>

void *
Volume::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP2::_castdown(cd) };
  return do_castdowns(casts,cd);
}

Volume::Volume(int dim):
  NLP2(dim),
  _value(fvalue),
  _gradient(grad),
  _hessian(Hessian),
  _interp_acc(1.0e-6)
{
}

Volume::~Volume()
{
}

// interpolate using the bisection algorithm
RefPoint
Volume::interpolate(RefPoint p1,RefPoint p2,double val)
{
  int dim = GetDim();
  
  DVector A(p1);
  DVector B(p2);

  SetX(A);
  double value0 = value() - val;

  SetX(B);
  double value1 = value() - val;

  if (value0*value1 > 0.0) {
      failure("interpolate(): values at endpoints don't bracket val");
    }
  else if (value0 == 0.0) {
      return p1;
    }
  else if (value1 == 0.0) {
      return p2;
    }

  DVector BA = B - A;

  double length = sqrt(BA.dot(BA));
  int niter = (int) (log(length/_interp_acc)/M_LN2);

  double f0 = 0.0;
  double f1 = 1.0;
  double fnext = 0.5;

  for (int i=0; i<niter; i++) {
      DVector X = A + fnext*BA;
      SetX(X);
      double valuenext = value() - val;

      if (valuenext*value0 <= 0.0) {
          value1 = valuenext;
          f1 = fnext;
          fnext = (f0 + fnext)*0.5;
        }
      else {
          value0 = valuenext;
          f0 = fnext;
          fnext = (fnext + f1)*0.5;
        }
    }

  DVector X = A + fnext*BA;
  RefPoint result = new Point(X);
  return result;
}

RefPoint
Volume::solve(RefPoint start,DVector& grad,double val)
{
  double direction;
  double startvalue = value(start);
  if (startvalue == val) return start;
  else if (startvalue < val) direction = 1.0;
  else direction = -1.0;
  DVector startv(start);
  int i=0;
  RefPoint next;
  double trialvalue;
  do {
      if (i>10) {
          fprintf(stderr,"Volume::solve: couldn't find end points\n");
          abort();
        }
      i++;
      DVector nextv = startv + (direction*i)*grad;
      next = new Point(nextv);
      trialvalue = value(next);
    } while ((startvalue-val)*(trialvalue-val)>0.0);
  return interpolate(start,next,val);
}

void
Volume::failure(const char * msg)
{
  fprintf(stderr,"Volume::failure: \"%s\"\n",msg);
  abort();
}

void
Volume::SetX(const Point& x)
{
  int dim = GetDim();
  ColumnVector tmp(dim);
  for (int i=0; i<dim; i++) tmp.element(i) = x[i];
  SetX(tmp);
}

void
Volume::gradient(DVector& g)
{
  int dim = GetDim();
  const ColumnVector grad = gradient();
  for (int i=0; i<dim; i++) g[i] = grad.element(i);
}

// the default point set is a uniform lattice
// others can be used by overriding this member
void
Volume::pointset(double resolution,
                 double valuemin, double valuemax,
                 SetRefPoint& points)
{
  int dim = GetDim();
  Point p1(dim),p2(dim);
  boundingbox(valuemin, valuemax, p1, p2);
  double* incr = new double [dim];

  int i;
  for (i=0; i<dim; i++) {
      incr[i] = resolution;
      if ((p2[i]-p1[i])/incr[i] < 3) incr[i] = (p2[i]-p1[i])/3;
    }

  if (dim!=3) {
      fprintf(stderr,"Volume::pointset: dim is != 3\n");
      abort();
    }
  
  Point p(dim);
  for (p[0] = p1[0]; p[0] < p2[0]+incr[0]; p[0] += incr[0]) {
      for (p[1] = p1[1]; p[1] < p2[1]+incr[1]; p[1] += incr[1]) {
          for (p[2] = p1[2]; p[2] < p2[2]+incr[2]; p[2] += incr[2]) {
              RefPoint rp(new Point(p));
              points.add(rp);
            }
        }
    }

  delete[] incr;
}

void
Volume::Eval()
{
  hessian(); gradient(); value();
}

double
Volume::EvalF()
{
  return value();
}

ColumnVector
Volume::EvalG()
{
  ColumnVector result = gradient();
  return result;
}

SymmetricMatrix
Volume::EvalH()
{
  SymmetricMatrix result = hessian();
  return result;
}

void
Volume::set_value(double e)
{
  _value = e;
  _value.computed() = 1;
}

double
Volume::value()
{
  return _value;
}

void
Volume::set_gradient(ColumnVector&g)
{
  _gradient = g;
  _gradient.computed() = 1;
}

void
Volume::set_gradient(DVector&g)
{
  ColumnVector& v = _gradient.result_noupdate();
  v.ReDimension(g.dim());
  for (int i=0; i<g.dim(); i++) v.element(i) = g[i];
  _gradient.computed() = 1;
}


const ColumnVector&
Volume::gradient()
{
  return _gradient;
}

void
Volume::set_hessian(SymmetricMatrix&h)
{
  _hessian = h;
  _hessian.computed() = 1;
}

const SymmetricMatrix&
Volume::hessian()
{
  return _hessian;
}

int
Volume::do_value()
{
  return _value.compute();
}

int
Volume::do_gradient()
{
  return _gradient.compute();
}

int
Volume::do_hessian()
{
  return _hessian.compute();
}

int
Volume::do_value(int f)
{
  return _value.compute(f);
}

int
Volume::do_gradient(int f)
{
  return _gradient.compute(f);
}

int
Volume::do_hessian(int f)
{
  return _hessian.compute(f);
}

SavableState_REF_def(Volume);
ARRAY_def(RefVolume);
SET_def(RefVolume);
ARRAYSET_def(RefVolume);
