
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <util/keyval/keyval.h>

#include <math/optimize/function.h>

////////////////////////////////////////////////////////////////////////

SavableState_REF_def(Function);

#define CLASSNAME Function
#define PARENTS virtual_base public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

void *
Function::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Function::Function():
  _value(this),
  _gradient(this),
  _hessian(this)
{
  _matrixkit = SCMatrixKit::default_matrixkit();
  _value.set_desired_accuracy(DBL_EPSILON);
  _gradient.set_desired_accuracy(DBL_EPSILON);
  _hessian.set_desired_accuracy(DBL_EPSILON);
}

Function::Function(const Function& func):
  _value(func._value,this),
  _gradient(func._gradient,this),
  _hessian(func._hessian,this)
{
  _matrixkit = func._matrixkit;
  _dim = func._dim;
  _x = func._x;
}

Function::Function(const RefKeyVal&kv):
  _value(this),
  _gradient(this),
  _hessian(this)
{
  _matrixkit = kv->describedclassvalue("matrixkit");

  if (_matrixkit.null()) _matrixkit = SCMatrixKit::default_matrixkit();

  _value.set_desired_accuracy(kv->doublevalue("value_accuracy"));
  if (_value.desired_accuracy() < DBL_EPSILON)
    _value.set_desired_accuracy(DBL_EPSILON);

  _gradient.set_desired_accuracy(kv->doublevalue("gradient_accuracy"));
  if (_gradient.desired_accuracy() < DBL_EPSILON)
    _gradient.set_desired_accuracy(DBL_EPSILON);

  _hessian.set_desired_accuracy(kv->doublevalue("hessian_accuracy"));
  if (_hessian.desired_accuracy() < DBL_EPSILON)
    _hessian.set_desired_accuracy(DBL_EPSILON);
}

Function::Function(StateIn&s):
  SavableState(s),
  _value(s,this),
  _gradient(s,this),
  _hessian(s,this)
{
  _matrixkit = SCMatrixKit::default_matrixkit();
  _dim.restore_state(s);
  _x = _matrixkit->vector(_dim);
  _x.restore(s);

  _gradient.result_noupdate() = matrixkit()->vector(_dim);
  _gradient.result_noupdate()->restore(s);

  _hessian.result_noupdate() = matrixkit()->symmmatrix(_dim);
  _hessian.result_noupdate()->restore(s);
}

Function::~Function()
{
}

Function &
Function::operator=(const Function& func)
{
  _matrixkit = func._matrixkit;
  _dim = func._dim;
  _x = func._x;
  _value = func._value;
  _gradient = func._gradient;
  _hessian = func._hessian;
  return *this;
}

void
Function::save_data_state(StateOut&s)
{
  _value.save_data_state(s);
  _dim.save_state(s);
  _x.save(s);
  _gradient.save_data_state(s);
  _gradient.result_noupdate()->save(s);
  _hessian.save_data_state(s);
  _hessian.result_noupdate()->save(s);
}

RefSCMatrixKit
Function::matrixkit()
{
  return _matrixkit;
}

RefSCDimension
Function::dimension()
{
  return _dim;
}

void
Function::set_x(const RefSCVector&v)
{
  _x.assign(v);
  obsolete();
}

double
Function::value()
{
  return _value;
}

int
Function::do_value()
{
  return _value.compute();
}

int
Function::do_value(int f)
{
  return _value.compute(f);
}

void
Function::set_value(double e)
{
  _value.result_noupdate() = e;
  _value.computed() = 1;
}

void
Function::set_desired_value_accuracy(double a)
{
  _value.set_desired_accuracy(a);
}

double
Function::desired_value_accuracy()
{
  return _value.desired_accuracy();
}

double
Function::actual_value_accuracy()
{
  return _value.actual_accuracy();
}

RefSCVector
Function::gradient()
{
  RefSCVector ret = _gradient.result();
  return ret;
}

int
Function::do_gradient()
{
  return _gradient.compute();
}

int
Function::do_gradient(int f)
{
  return _gradient.compute(f);
}

void
Function::set_gradient(RefSCVector&g)
{
  _gradient.result_noupdate() = g;
  _gradient.computed() = 1;
}

void
Function::set_desired_gradient_accuracy(double a)
{
  _gradient.set_desired_accuracy(a);
}

double
Function::actual_gradient_accuracy()
{
  return _gradient.actual_accuracy();
}

double
Function::desired_gradient_accuracy()
{
  return _gradient.desired_accuracy();
}

RefSymmSCMatrix
Function::hessian()
{
  return _hessian;
}

int
Function::do_hessian()
{
  return _hessian.compute();
}

int
Function::do_hessian(int f)
{
  return _hessian.compute(f);
}

void
Function::set_hessian(RefSymmSCMatrix&h)
{
  _hessian.result_noupdate() = h;
  _hessian.computed() = 1;
}

// the default guess hessian is the unit diagonal
void
Function::guess_hessian(RefSymmSCMatrix&hessian)
{
  RefSCElementOp op(new SCElementShiftDiagonal(1.0));
  hessian.assign(0.0);
  hessian.element_op(op);
}

RefSymmSCMatrix
Function::inverse_hessian(RefSymmSCMatrix&hessian)
{
  return hessian.gi();
}

void
Function::set_desired_hessian_accuracy(double a)
{
  _hessian.set_desired_accuracy(a);
}

double
Function::desired_hessian_accuracy()
{
  return _hessian.desired_accuracy();
}

double
Function::actual_hessian_accuracy()
{
  return _hessian.actual_accuracy();
}

void
Function::print(ostream&o)
{
  o << indent << "value_accuracy = " << desired_value_accuracy() << endl;
  o << indent << "gradient.accuracy = "
              << desired_gradient_accuracy()
              << endl;
  o << indent << "hessian_accuracy = " << desired_hessian_accuracy() << endl;
}

void
Function::set_matrixkit(const RefSCMatrixKit& kit)
{
  _matrixkit = kit;
}

void
Function::set_dimension(const RefSCDimension& dim)
{
  _dim = dim;
  _x = _matrixkit->vector(dim);
  _gradient = matrixkit()->vector(dim);
  _hessian = matrixkit()->symmmatrix(dim);
}

int
Function::value_implemented()
{
  return 0;
}

int
Function::gradient_implemented()
{
  return 0;
}

int
Function::hessian_implemented()
{
  return 0;
}
