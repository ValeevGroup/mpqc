
#ifndef _math_optimize_transform_h
#define _math_optimize_transform_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/matrix.h>

//texi The @code{NonlinearTransform} class transforms between
// two nonlinear coordinate systems.  It is needed when a change
// of coordinates occurs in the middle of an optimization.
class NonlinearTransform: public VRefCount {
  protected:
    // The linear part of the nonlinear transform.  This must
    // be initialized by derived classes in their
    // transform_coordinates routine (or the transform
    // members must be overridden so it is ignored).
    RefSCMatrix linear_transform_;
  public:
    ~NonlinearTransform();

    //texi Transform the coordinates.
    virtual void transform_coordinates(const RefSCVector& x) = 0;
    //texi Transform the gradient at a point in the new
    // coordinate system.  @code{transform\_coordinates} must be
    // called first to give the point.
    virtual void transform_gradient(const RefSCVector& g);
    //texi Transform the hessian to the new coordinate system.
    // @code{transform\_gradient} must be called first to
    // initialize this routine.
    virtual void transform_hessian(const RefSymmSCMatrix& h);
    //texi Transform the inverse of the hessian.
    // @code{transform\_gradient} must be called first to
    // initialize this routine.
    virtual void transform_ihessian(const RefSymmSCMatrix &ih);
};

REF_dec(NonlinearTransform);

//texi The @code{IdentityTransform} is a special case of
// @code{NonlinearTransform} were no transformation takes place.
class IdentityTransform: public NonlinearTransform {
  public:
    ~IdentityTransform();

    //texi These override the tranformation members of
    // @code{NonlinearTransform} and do nothing.
    void transform_coordinates(const RefSCVector& x);
    void transform_gradient(const RefSCVector& g);
    void transform_hessian(const RefSymmSCMatrix& h);
    void transform_ihessian(const RefSymmSCMatrix &ih);
};

#endif
