
#ifndef _math_optimize_qnewton_h
#define _math_optimize_qnewton_h

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
// newton and related methods

class QNewtonOpt: public Optimize {
#   define CLASSNAME QNewtonOpt
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double maxabs_gradient;
    double convergence_;
    double accuracy_;

    RefSymmSCMatrix ihessian_;
    RefHessianUpdate update_;
    RefLineOpt lineopt_;

    int take_newton_step_;
  public:
    QNewtonOpt(const RefKeyVal&);
    QNewtonOpt(StateIn&);
    ~QNewtonOpt();
    void save_data_state(StateOut&);

    void init();
    int update();
};

#endif
