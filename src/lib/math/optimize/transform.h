//
// transform.h
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

#ifndef _math_optimize_transform_h
#define _math_optimize_transform_h

#include <math/scmat/matrix.h>

namespace sc {

/** The NonlinearTransform class transforms between
    two nonlinear coordinate systems.  It is needed when a change
    of coordinates occurs in the middle of an optimization. */
class NonlinearTransform: public RefCount {
  protected:
    // The linear part of the nonlinear transform.  This must
    // be initialized by derived classes in their
    // transform_coordinates routine (or the transform
    // members must be overridden so it is ignored).
    RefSCMatrix linear_transform_;
  public:
    ~NonlinearTransform();

    /// Transform the coordinates.
    virtual void transform_coordinates(const RefSCVector& x) = 0;
    /** Transform the gradient at a point in the new
        coordinate system.  transform_coordinates must be
        called first to give the point. */
    virtual void transform_gradient(const RefSCVector& g);
    /** Transform the hessian to the new coordinate system.
        transform_gradient must be called first to
        initialize this routine. */
    virtual void transform_hessian(const RefSymmSCMatrix& h);
    /** Transform the inverse of the hessian.
        transform_gradient must be called first to
        initialize this routine. */
    virtual void transform_ihessian(const RefSymmSCMatrix &ih);
};



/** The IdentityTransform is a special case of
    NonlinearTransform were no transformation takes place.
*/
class IdentityTransform: public NonlinearTransform {
  public:
    ~IdentityTransform();

    /** These override the tranformation members of
        NonlinearTransform and do nothing. */
    void transform_coordinates(const RefSCVector& x);
    void transform_gradient(const RefSCVector& g);
    void transform_hessian(const RefSymmSCMatrix& h);
    void transform_ihessian(const RefSymmSCMatrix &ih);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
