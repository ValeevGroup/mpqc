
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
# include <stdio.h>
# include <math.h>
};

#include <util/keyval/keyval.h>
#include <math/scmat/vector3.h>
#include <math/scmat/local.h>
#include <math/isosurf/volume.h>

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

Volume::Volume():
  NLP2(new LocalSCDimension(3)),
  _interp_acc(1.0e-6)
{
}

Volume::Volume(const RefKeyVal&keyval):
  NLP2(keyval)
{
  _interp_acc = keyval->doublevalue("interpolation_accuracy");
  if (keyval->error() != KeyVal::OK) _interp_acc = 1.0e-6;
}

Volume::~Volume()
{
}

// interpolate using the bisection algorithm
void
Volume::interpolate(const SCVector3& A,
                    const SCVector3& B,
                    double val,
                    SCVector3& result)
{
  set_x(A);
  double value0 = value() - val;

  set_x(B);
  double value1 = value() - val;

  if (value0*value1 > 0.0) {
      failure("interpolate(): values at endpoints don't bracket val");
    }
  else if (value0 == 0.0) {
      result = A;
      return;
    }
  else if (value1 == 0.0) {
      result = B;
      return;
    }

  SCVector3 BA = B - A;

  double length = sqrt(BA.dot(BA));
  int niter = (int) (log(length/_interp_acc)/M_LN2);

  double f0 = 0.0;
  double f1 = 1.0;
  double fnext = 0.5;

  for (int i=0; i<niter; i++) {
      SCVector3 X = A + fnext*BA;
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

  result = A + fnext*BA;
}

void
Volume::solve(const SCVector3& start,
              const SCVector3& grad,
              double val,
              SCVector3& result)
{
  double direction;
  set_x(start);
  double startvalue = value();
  if (startvalue == val) {
      result = start;
      return;
    }
  else if (startvalue < val) direction = 1.0;
  else direction = -1.0;
  int i=0;
  SCVector3 next;
  double trialvalue;
  do {
      if (i>10) {
          fprintf(stderr,"Volume::solve: couldn't find end points\n");
          abort();
        }
      i++;
      next = start + (direction*i)*grad;
      set_x(next);
      trialvalue = value();
    } while ((startvalue-val)*(trialvalue-val)>0.0);
  interpolate(start,next,val,result);
}

void
Volume::failure(const char * msg)
{
  fprintf(stderr,"Volume::failure: \"%s\"\n",msg);
  abort();
}

void
Volume::set_gradient(const SCVector3& g)
{
  RefSCVector& grad = _gradient.result_noupdate();
  grad.set_element(0, g[0]);
  grad.set_element(1, g[1]);
  grad.set_element(2, g[2]);
  _gradient.computed() = 1;
}

void
Volume::get_gradient(SCVector3& g)
{
  const RefSCVector v = gradient();
  g[0] = v.get_element(0);
  g[1] = v.get_element(1);
  g[2] = v.get_element(2);
}

void
Volume::set_x(const SCVector3& x)
{
  _x.set_element(0, x[0]);
  _x.set_element(1, x[1]);
  _x.set_element(2, x[2]);
  obsolete();
}

void
Volume::get_x(SCVector3& x)
{
  const RefSCVector& v = get_x_no_copy();
  x[0] = v.get_element(0);
  x[1] = v.get_element(1);
  x[2] = v.get_element(2);
}

SavableState_REF_def(Volume);
ARRAY_def(RefVolume);
SET_def(RefVolume);
ARRAYSET_def(RefVolume);
