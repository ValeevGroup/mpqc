//
// volume.h
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

#ifndef _math_isosurf_volume_h
#define _math_isosurf_volume_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/optimize/function.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>

namespace sc {

/** A Volume is a Function of three variables. */
class Volume: public Function {
  private:
    double _interp_acc;
  protected:
    double& interpolation_accuracy();

    virtual void compute() = 0;

    virtual void failure(const char*);
  public:
    Volume();
    Volume(const Ref<KeyVal>&);
    ~Volume();

    void set_gradient(const SCVector3& g);
    void set_gradient(RefSCVector& g);
    void get_gradient(SCVector3& g);
    void set_x(const SCVector3& x);
    void set_x(const RefSCVector& x);
    void get_x(SCVector3& x);

    // find the corners of a bounding box which approximately
    // contains all points with a value between valuemin and valuemax
    // the result must satisfy p1[i] < p2[i]
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2) = 0;

    virtual void interpolate(const SCVector3& p1,
                             const SCVector3& p2,
                             double value,
                             SCVector3& result);
    virtual void solve(const SCVector3& p,
                       const SCVector3& grad,
                       double value,
                       SCVector3& result);
};

}

#ifdef INLINE_FUNCTIONS
#include <math/isosurf/volume_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
