
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
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
  value_(this),
  gradient_(this),
  hessian_(this)
{
  matrixkit_ = SCMatrixKit::default_matrixkit();
  value_.set_desired_accuracy(DBL_EPSILON);
  gradient_.set_desired_accuracy(DBL_EPSILON);
  hessian_.set_desired_accuracy(DBL_EPSILON);
}

Function::Function(const Function& func):
  value_(func.value_,this),
  gradient_(func.gradient_,this),
  hessian_(func.hessian_,this)
{
  matrixkit_ = func.matrixkit_;
  dim_ = func.dim_;
  x_ = func.x_;
}

Function::Function(const RefKeyVal&kv):
  value_(this),
  gradient_(this),
  hessian_(this)
{
  matrixkit_ = kv->describedclassvalue("matrixkit");

  if (matrixkit_.null()) matrixkit_ = SCMatrixKit::default_matrixkit();

  value_.set_desired_accuracy(kv->doublevalue("value_accuracy"));
  if (value_.desired_accuracy() < DBL_EPSILON)
    value_.set_desired_accuracy(DBL_EPSILON);

  gradient_.set_desired_accuracy(kv->doublevalue("gradient_accuracy"));
  if (gradient_.desired_accuracy() < DBL_EPSILON)
    gradient_.set_desired_accuracy(DBL_EPSILON);

  hessian_.set_desired_accuracy(kv->doublevalue("hessian_accuracy"));
  if (hessian_.desired_accuracy() < DBL_EPSILON)
    hessian_.set_desired_accuracy(DBL_EPSILON);
}

Function::Function(StateIn&s):
  SavableState(s),
  value_(s,this),
  gradient_(this),
  hessian_(this)
{
  matrixkit_ = SCMatrixKit::default_matrixkit();
  dim_.restore_state(s);
  x_ = matrixkit_->vector(dim_);
  x_.restore(s);

  gradient_.result_noupdate() = matrixkit()->vector(dim_);
  gradient_.restore_state(s);
  gradient_.result_noupdate().restore(s);

  hessian_.result_noupdate() = matrixkit()->symmmatrix(dim_);
  hessian_.restore_state(s);
  hessian_.result_noupdate().restore(s);
}

Function::~Function()
{
}

Function &
Function::operator=(const Function& func)
{
  matrixkit_ = func.matrixkit_;
  dim_ = func.dim_;
  x_ = func.x_;
  value_ = func.value_;
  gradient_ = func.gradient_;
  hessian_ = func.hessian_;
  return *this;
}

void
Function::save_data_state(StateOut&s)
{
  value_.save_data_state(s);
  dim_.save_state(s);
  x_.save(s);
  gradient_.save_data_state(s);
  gradient_.result_noupdate().save(s);
  hessian_.save_data_state(s);
  hessian_.result_noupdate().save(s);
}

RefSCMatrixKit
Function::matrixkit()
{
  return matrixkit_;
}

RefSCDimension
Function::dimension()
{
  return dim_;
}

void
Function::set_x(const RefSCVector&v)
{
  x_.assign(v);
  obsolete();
}

double
Function::value()
{
  return value_;
}

int
Function::value_needed()
{
  return value_.needed();
}

int
Function::do_value()
{
  return value_.compute();
}

int
Function::do_value(int f)
{
  return value_.compute(f);
}

void
Function::set_value(double e)
{
  value_.result_noupdate() = e;
  value_.computed() = 1;
}

void
Function::set_desired_value_accuracy(double a)
{
  value_.set_desired_accuracy(a);
}

void
Function::set_actual_value_accuracy(double a)
{
  value_.set_actual_accuracy(a);
}

double
Function::desired_value_accuracy()
{
  return value_.desired_accuracy();
}

double
Function::actual_value_accuracy()
{
  return value_.actual_accuracy();
}

RefSCVector
Function::gradient()
{
  RefSCVector ret = gradient_.result();
  return ret;
}

int
Function::gradient_needed()
{
  return gradient_.needed();
}

int
Function::do_gradient()
{
  return gradient_.compute();
}

int
Function::do_gradient(int f)
{
  return gradient_.compute(f);
}

void
Function::set_gradient(RefSCVector&g)
{
  gradient_.result_noupdate() = g;
  gradient_.computed() = 1;
}

void
Function::set_desired_gradient_accuracy(double a)
{
  gradient_.set_desired_accuracy(a);
}

void
Function::set_actual_gradient_accuracy(double a)
{
  gradient_.set_actual_accuracy(a);
}

double
Function::actual_gradient_accuracy()
{
  return gradient_.actual_accuracy();
}

double
Function::desired_gradient_accuracy()
{
  return gradient_.desired_accuracy();
}

RefSymmSCMatrix
Function::hessian()
{
  return hessian_.result();
}

int
Function::hessian_needed()
{
  return hessian_.needed();
}

int
Function::do_hessian()
{
  return hessian_.compute();
}

int
Function::do_hessian(int f)
{
  return hessian_.compute(f);
}

void
Function::set_hessian(RefSymmSCMatrix&h)
{
  hessian_.result_noupdate() = h;
  hessian_.computed() = 1;
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
  hessian_.set_desired_accuracy(a);
}

void
Function::set_actual_hessian_accuracy(double a)
{
  hessian_.set_actual_accuracy(a);
}

double
Function::desired_hessian_accuracy()
{
  return hessian_.desired_accuracy();
}

double
Function::actual_hessian_accuracy()
{
  return hessian_.actual_accuracy();
}

void
Function::print(ostream&o)
{
  if (matrixkit_->messagegrp()->me()==0) {
    o << indent << "Function Parameters:\n" << incindent;
    o << indent << "value_accuracy = " << desired_value_accuracy() << endl;
    o << indent << "gradient.accuracy = "
      << desired_gradient_accuracy()
      << endl;
    o << indent << "hessian_accuracy = " << desired_hessian_accuracy() << endl;
    o << decindent << endl;
  }
}

void
Function::set_matrixkit(const RefSCMatrixKit& kit)
{
  matrixkit_ = kit;
}

void
Function::set_dimension(const RefSCDimension& dim)
{
  dim_ = dim;
  x_ = matrixkit_->vector(dim);
  gradient_ = matrixkit()->vector(dim);
  hessian_ = matrixkit()->symmmatrix(dim);
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
