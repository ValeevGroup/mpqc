
#ifndef _math_optimize_nlp_h
#define _math_optimize_nlp_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/scostream.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/result.h>

class NLP0: virtual_base public SavableState, public Compute {
#   define CLASSNAME NLP0
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit _matrixkit;
  protected:
    RefSCDimension _dim;
    RefSCVector    _x;    // variables
    AccResultdouble   _value;// value of function at _x
    virtual void set_value(double);

    virtual void set_matrixkit(const RefSCMatrixKit&);
    virtual void set_dimension(const RefSCDimension&);
  public:
    NLP0();
    NLP0(StateIn&);
    NLP0(const RefKeyVal&);
    virtual ~NLP0();

    RefSCMatrixKit matrixkit();
    RefSCDimension dimension();

    virtual void save_data_state(StateOut&);

    virtual double value();
    int do_value(int);
    int do_value();

    virtual void set_x(const RefSCVector&);
    RefSCVector get_x() const;
    const RefSCVector& get_x_no_copy() const;

    virtual void set_desired_value_accuracy(double);
    virtual double actual_value_accuracy();
    virtual double desired_value_accuracy();

    virtual void print(SCostream& =SCostream::cout);
};
SavableState_REF_dec(NLP0);

inline RefSCVector
NLP0::get_x() const
{
  return _x.copy();
}

inline const RefSCVector&
NLP0::get_x_no_copy() const
{
  return _x;
}

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
    AccResultRefSCVector _gradient; // gradient at _x
    virtual void set_gradient(RefSCVector&);
  public:
    NLP1();
    NLP1(StateIn&);
    NLP1(const RefKeyVal&);
    virtual ~NLP1();

    virtual void save_data_state(StateOut&);

    virtual RefSCVector gradient();
    int do_gradient(int);
    int do_gradient();

    virtual void set_desired_gradient_accuracy(double);
    virtual double actual_gradient_accuracy();
    virtual double desired_gradient_accuracy();

    // gradients by values at finite displacements
    // virtual RefSCVector fd0_gradient();

    void set_dimension(const RefSCDimension&);

    virtual void print(SCostream& =SCostream::cout);
};
SavableState_REF_dec(NLP1);

class NLP2: public NLP1 {
#   define CLASSNAME NLP2
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    AccResultRefSymmSCMatrix _hessian;
    virtual void set_hessian(RefSymmSCMatrix&);
  public:
    NLP2();
    NLP2(StateIn&);
    NLP2(const RefKeyVal&);
    virtual ~NLP2();

    void save_data_state(StateOut&);

    virtual RefSymmSCMatrix hessian();
    int do_hessian(int);
    int do_hessian();

    virtual void set_desired_hessian_accuracy(double);
    virtual double actual_hessian_accuracy();
    virtual double desired_hessian_accuracy();

    // hessian by gradients at finite displacements
    // virtual RefSCMatrix fd1_hessian();

    // quick, approximate hessian
    virtual void guess_hessian(RefSymmSCMatrix&);
    virtual RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    virtual void print(SCostream& =SCostream::cout);

    void set_dimension(const RefSCDimension&);

    // information about the availability of values, gradients, and hessians
    virtual int value_implemented();
    virtual int gradient_implemented();
    virtual int hessian_implemented();
};
SavableState_REF_dec(NLP2);

#endif
