//
// volume.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/math.h>
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/vector3.h>
#include <math/scmat/local.h>
#include <math/isosurf/volume.h>

using namespace std;
using namespace sc;

static ClassDesc Volume_cd(
  typeid(Volume),"Volume",1,"public Function",
  0, 0, 0);

Volume::Volume():
  _interp_acc(1.0e-6)
{
  set_dimension(new SCDimension(3));
}

Volume::Volume(const Ref<KeyVal>&keyval):
  Function(keyval)
{
  set_dimension(new SCDimension(3));
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
          ExEnv::errn() << "Volume::solve: couldn't find end points" << endl;
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
  ExEnv::errn() << scprintf("Volume::failure: \"%s\"\n",msg);
  abort();
}

void
Volume::set_gradient(const SCVector3& g)
{
  RefSCVector grad(dimension(), matrixkit());
  grad.set_element(0, g[0]);
  grad.set_element(1, g[1]);
  grad.set_element(2, g[2]);
  set_gradient(grad);
}

void
Volume::set_gradient(RefSCVector& g)
{
  Function::set_gradient(g);
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
  RefSCVector xx(dimension(), matrixkit());
  xx.set_element(0, x[0]);
  xx.set_element(1, x[1]);
  xx.set_element(2, x[2]);
  set_x(xx);
}

void
Volume::set_x(const RefSCVector& x)
{
  Function::set_x(x);
}

void
Volume::get_x(SCVector3& x)
{
  const RefSCVector& v = get_x_no_copy();
  x[0] = v.get_element(0);
  x[1] = v.get_element(1);
  x[2] = v.get_element(2);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
