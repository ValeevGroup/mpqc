
#ifndef _math_optimize_nlp_h
#define _math_optimize_nlp_h

#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/result.h>

class NLP0: virtual public SavableState, public Compute {
#   define CLASSNAME NLP0
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefSCDimension _dim;
    RefSCVector    _x;    // variables
    Resultdouble   _value;// value of function at _x
    double       _acc_value;    // value should be evaluated with this accuracy
    double       _maxacc_value; // maximum value for _acc_value
    virtual void set_value(double);
  public:
    NLP0(RefSCDimension&);
    NLP0(StateIn&);
    NLP0(KeyVal&);
    virtual ~NLP0();

    virtual RefSCDimension dimension();

    virtual void save_data_state(StateOut&);

    virtual double value();
    int do_value(int);
    int do_value();

    virtual void set_x(int i, double);
    virtual void set_x(RefSCVector&);
    virtual RefSCVector get_x();
};
SavableState_REF_dec(NLP0);

//------------------------------------------------------------------------
// Derived Classes from NonLinear Programming Problem
// NLP0: NLP  + No derivative info
// NLP1: NLP0 + First derivatives
// NLP2: NLP1 + Second derivatives
//------------------------------------------------------------------------

class NLP1: public NLP0 {
#   define CLASSNAME NLP1
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    ResultRefSCVector _gradient;     // gradient at _x
    virtual void set_gradient(RefSCVector&);
  public:
    NLP1(RefSCDimension&);
    NLP1(StateIn&);
    NLP1(KeyVal&);
    virtual ~NLP1();

    virtual void save_data_state(StateOut&);

    virtual RefSCVector gradient();
    int do_gradient(int);
    int do_gradient();

    // gradients by values at finite displacements
    // virtual RefSCVector fd0_gradient();
};
SavableState_REF_dec(NLP1);

class NLP2: public NLP1 {
#   define CLASSNAME NLP2
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    ResultRefSymmSCMatrix _hessian;
    virtual void set_hessian(RefSymmSCMatrix&);
  public:
    NLP2(RefSCDimension&);
    NLP2(StateIn&);
    NLP2(KeyVal&);
    virtual ~NLP2();

    void save_data_state(StateOut&);

    virtual RefSymmSCMatrix hessian();
    int do_hessian(int);
    int do_hessian();

    // hessian by gradients at finite displacements
    // virtual RefSCMatrix fd1_hessian();
};
SavableState_REF_dec(NLP2);

#endif
