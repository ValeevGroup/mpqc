
#ifndef _math_isosurf_volume_h
#define _math_isosurf_volume_h

#include <math/optimize/nlp.h>
#include <util/container/ref.h>
#include <math/scmat/matrix.h>

class Volume: public NLP2 {
#   define CLASSNAME Volume
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double _interp_acc;
  protected:
    void set_value(double);
    void set_gradient(RefSCVector&);
    void set_hessian(RefSymmSCMatrix&);

    double& interpolation_accuracy();

    virtual void compute() = 0;

    virtual void failure(const char*);
  public:
    Volume(RefSCDimension&);
    ~Volume();

    // find the corners of a bounding box which approximately
    // contains all points with a value between valuemin and valuemax
    // the result must satisfy p1[i] < p2[i]
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             RefSCVector& p1, RefSCVector& p2) = 0;
    virtual void pointset(double resolution,
                          double valuemin,
                          double valuemax,
                          SetRefSCVector& points);

    virtual RefSCVector interpolate(RefSCVector&,RefSCVector&,double value);
    virtual RefSCVector solve(RefSCVector&,RefSCVector&,double value);
    
};

SavableState_REF_dec(Volume);
ARRAY_dec(RefVolume);
SET_dec(RefVolume);
ARRAYSET_dec(RefVolume);

#ifdef INLINE_FUNCTIONS
#include <math/isosurf/volume_i.h>
#endif

#endif
