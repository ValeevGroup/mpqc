
#ifndef _math_optimize_opt_h
#define _math_optimize_opt_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/class/class.h>
#include <math/scmat/matrix.h>
#include <math/optimize/nlp.h>

////////////////////////////////////////////////////////////////////////

class Optimize: virtual public SavableState {
#   define CLASSNAME Optimize
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int max_iterations_;
    int n_iterations_;
  public:
    Optimize();
    Optimize(StateIn&);
    Optimize(const RefKeyVal&);
    virtual ~Optimize();
    void save_data_state(StateOut&);

    virtual int optimize();

    // initialize the optimizer
    virtual void init();
    // take a step
    // returns 1 if the optimization has converged, otherwise 0
    virtual int update() = 0;

    // returns the problem being optimized
    virtual RefNLP0 nlp() = 0;
};
SavableState_REF_dec(Optimize);

class LineOpt: public Optimize {
#   define CLASSNAME LineOpt
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefSCVector search_direction_;
  public:
    LineOpt();
    LineOpt(StateIn&);
    LineOpt(const RefKeyVal&);
    ~LineOpt();
    void save_data_state(StateOut&);

    void set_search_direction(RefSCVector&);
};
SavableState_REF_dec(LineOpt);

////////////////////////////////////////////////////////////////////////
//  hessian update classes.  based on the value of inverse_hessian_
//  x and g may be reversed (see Schlegel, ab initio Methods in Quantum
//  Chemistry I, 1987, p 10

class HessianUpdate: virtual public SavableState {
#   define CLASSNAME HessianUpdate
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int inverse_hessian_;
  public:
    HessianUpdate();
    HessianUpdate(StateIn&);
    HessianUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    virtual ~HessianUpdate();
    virtual void update(RefSymmSCMatrix&hessian,RefNLP2&nlp,
                        RefSCVector&xnew,RefSCVector&gnew) = 0;
    void set_inverse();
};
SavableState_REF_dec(HessianUpdate);

class DFPUpdate: public HessianUpdate {
#   define CLASSNAME DFPUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    DFPUpdate();
    DFPUpdate(StateIn&);
    DFPUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~DFPUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
};

class BFGSUpdate: public DFPUpdate {
#   define CLASSNAME BFGSUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    BFGSUpdate();
    BFGSUpdate(StateIn&);
    BFGSUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~BFGSUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
};

class PowellUpdate: public HessianUpdate {
#   define CLASSNAME PowellUpdate
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCVector xprev;
    RefSCVector gprev;
  public:
    PowellUpdate();
    PowellUpdate(StateIn&);
    PowellUpdate(const RefKeyVal&);
    void save_data_state(StateOut&);
    ~PowellUpdate();
    void update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                RefSCVector&xnew,RefSCVector&gnew);
};

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

    RefNLP2 nlp_;
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
    RefNLP0 nlp();
};

#endif
