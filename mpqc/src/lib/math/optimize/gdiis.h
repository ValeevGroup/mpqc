
#ifndef _math_optimize_gdiis_h
#define _math_optimize_gdiis_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/function.h>
#include <math/optimize/opt.h>
#include <math/optimize/update.h>

////////////////////////////////////////////////////////////////////////
// gdiis

class GDIISOpt: public Optimize {
#   define CLASSNAME GDIISOpt
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int nsave;
    int diis_iter;
  
    double maxabs_gradient;
    double convergence_;
    double accuracy_;

    RefSCVector *coords_;
    RefSCVector *grad_;
    RefSCVector *error_;

    RefSymmSCMatrix ihessian_;
    RefHessianUpdate update_;

  public:
    GDIISOpt(const RefKeyVal&);
    GDIISOpt(StateIn&);
    ~GDIISOpt();
    void save_data_state(StateOut&);

    void init();
    int update();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
