
#ifndef _math_isosurf_volume_h
#define _math_isosurf_volume_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/optimize/nlp.h>
#include <util/container/ref.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>

class Volume: public NLP2 {
#   define CLASSNAME Volume
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    double _interp_acc;
  protected:
    double& interpolation_accuracy();

    virtual void compute() = 0;

    virtual void failure(const char*);
  public:
    Volume();
    Volume(const RefKeyVal&);
    ~Volume();

    void set_gradient(const SCVector3& g);
    void get_gradient(SCVector3& g);
    void set_x(const SCVector3& x);
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

SavableState_REF_dec(Volume);
ARRAY_dec(RefVolume);
SET_dec(RefVolume);
ARRAYSET_dec(RefVolume);

#ifdef INLINE_FUNCTIONS
#include <math/isosurf/volume_i.h>
#endif

#endif
