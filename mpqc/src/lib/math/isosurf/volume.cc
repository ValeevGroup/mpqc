
extern "C" {
# include <stdio.h>
# include <math.h>
};

#include <math/scmat/vector3.h>
#include "volume.h"

#define CLASSNAME Volume
#define PARENTS public NLP2
#include <util/state/statei.h>
#include <util/class/classia.h>

void *
Volume::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = NLP2::_castdown(cd);
  return do_castdowns(casts,cd);
}

Volume::Volume(RefSCDimension& dim):
  NLP2(dim),
  _interp_acc(1.0e-6)
{
}

Volume::~Volume()
{
}

// interpolate using the bisection algorithm
RefSCVector
Volume::interpolate(RefSCVector& p1,RefSCVector& p2,double val)
{
  RefSCVector A(p1);
  RefSCVector B(p2);

  set_x(A);
  double value0 = value() - val;

  set_x(B);
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

  RefSCVector BA = B - A;

  double length = sqrt(BA.dot(BA));
  int niter = (int) (log(length/_interp_acc)/M_LN2);

  double f0 = 0.0;
  double f1 = 1.0;
  double fnext = 0.5;

  for (int i=0; i<niter; i++) {
      RefSCVector X = A + fnext*BA;
      set_x(X);
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

  RefSCVector X = A + fnext*BA;
  return X;
}

RefSCVector
Volume::solve(RefSCVector& start,RefSCVector& grad,double val)
{
  double direction;
  double startvalue = (set_x(start),value());
  if (startvalue == val) return start;
  else if (startvalue < val) direction = 1.0;
  else direction = -1.0;
  int i=0;
  RefSCVector next;
  double trialvalue;
  do {
      if (i>10) {
          fprintf(stderr,"Volume::solve: couldn't find end points\n");
          abort();
        }
      i++;
      next = start + (direction*i)*grad;
      trialvalue = (set_x(next),value());
    } while ((startvalue-val)*(trialvalue-val)>0.0);
  return interpolate(start,next,val);
}

void
Volume::failure(const char * msg)
{
  fprintf(stderr,"Volume::failure: \"%s\"\n",msg);
  abort();
}

// the default point set is a uniform lattice
// others can be used by overriding this member
void
Volume::pointset(double resolution,
                 double valuemin, double valuemax,
                 SetRefSCVector& points)
{
  if (dimension().n()!=3) {
      fprintf(stderr,"Volume::pointset: dim is != 3\n");
      abort();
    }

  int dim = dimension().n();
  RefSCVector p1(dimension()),p2(dimension());
  boundingbox(valuemin, valuemax, p1, p2);
  double* incr = new double [dim];

  int i;
  for (i=0; i<dim; i++) {
      incr[i] = resolution;
      if ((p2(i)-p1(i))/incr[i] < 3) incr[i] = (p2(i)-p1(i))/3;
    }
  
  SCVector3 p;
  SCVector3 p1_(p1);
  SCVector3 p2_(p2);
  for (p(0) = p1_(0); p(0) < p2_(0)+incr[0]; p(0) += incr[0]) {
      for (p(1) = p1_(1); p(1) < p2_(1)+incr[1]; p(1) += incr[1]) {
          for (p(2) = p1_(2); p(2) < p2_(2)+incr[2]; p(2) += incr[2]) {
              RefSCVector rp(dimension());
              rp.assign(p.data());
              points.add(rp);
            }
        }
    }

  delete[] incr;
}

SavableState_REF_def(Volume);
ARRAY_def(RefVolume);
SET_def(RefVolume);
ARRAYSET_def(RefVolume);
