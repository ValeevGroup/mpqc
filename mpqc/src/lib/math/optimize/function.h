
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_optimize_function_h
#define _math_optimize_function_h

#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/result.h>

class Function: virtual_base public SavableState, public Compute {
#   define CLASSNAME Function
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit _matrixkit;
  protected:
    RefSCDimension _dim;
    RefSCVector    _x;    // variables
    AccResultdouble   _value;// value of function at _x
    virtual void set_value(double);

    AccResultRefSCVector _gradient; // gradient at _x
    virtual void set_gradient(RefSCVector&);

    AccResultRefSymmSCMatrix _hessian;
    virtual void set_hessian(RefSymmSCMatrix&);

    virtual void set_matrixkit(const RefSCMatrixKit&);
    virtual void set_dimension(const RefSCDimension&);
  public:
    Function();
    Function(StateIn&);
    Function(const Function&);
    Function(const RefKeyVal&);
    virtual ~Function();

    Function & operator=(const Function&);
    
    RefSCMatrixKit matrixkit();
    RefSCDimension dimension();

    virtual void save_data_state(StateOut&);

    virtual double value();
    int do_value(int);
    int do_value();

    virtual RefSCVector gradient();
    int do_gradient(int);
    int do_gradient();

    virtual void set_desired_gradient_accuracy(double);
    virtual double actual_gradient_accuracy();
    virtual double desired_gradient_accuracy();

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

    // information about the availability of values, gradients, and hessians
    virtual int value_implemented();
    virtual int gradient_implemented();
    virtual int hessian_implemented();

    virtual void set_x(const RefSCVector&);
    RefSCVector get_x() const { return _x.copy(); }
    const RefSCVector& get_x_no_copy() const { return _x; }

    virtual void set_desired_value_accuracy(double);
    virtual double actual_value_accuracy();
    virtual double desired_value_accuracy();

    virtual void print(ostream& = cout);
};
SavableState_REF_dec(Function);

#endif
