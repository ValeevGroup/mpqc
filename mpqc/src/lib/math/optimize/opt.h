
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

#endif
