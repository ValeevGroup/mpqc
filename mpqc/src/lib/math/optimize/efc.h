
#ifndef _math_optimize_efc_h
#define _math_optimize_efc_h

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
// eigenvector following a la Baker (JCC Vol 7, No 4, 385-395, 1986)

class EFCOpt: public Optimize {
#   define CLASSNAME EFCOpt
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int tstate;
    int modef;
  
    double maxabs_gradient;
    double convergence_;
    double accuracy_;

    RefSymmSCMatrix hessian_;
    RefHessianUpdate update_;
    RefSCVector last_mode_;

  public:
    EFCOpt(const RefKeyVal&);
    EFCOpt(StateIn&);
    ~EFCOpt();
    void save_data_state(StateOut&);

    void init();
    int update();
};

#endif
