
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math/scmat/matrix.h>
#include <util/keyval/keyval.h>

#include "nlp.h"

////////////////////////////////////////////////////////////////////////

SavableState_REF_def(NLP0);

#define CLASSNAME NLP0
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
NLP0::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

NLP0::NLP0(RefSCDimension&dim):
  _value(this),
  _dim(dim),
  _x(dim)
{
  _value.set_desired_accuracy(DBL_EPSILON);
}

NLP0::NLP0(KeyVal&kv):
  _value(this)
{
  _dim = kv.describedclassvalue("dimension");

  _value.set_desired_accuracy(kv.doublevalue("value_accuracy"));
  if (_value.desired_accuracy() < DBL_EPSILON)
    _value.set_desired_accuracy(DBL_EPSILON);

  RefSCVector x(_dim);
  _x = x;
}

NLP0::NLP0(StateIn&s):
  SavableState(s,class_desc_),
  _value(s,this)
{
  _dim.restore_state(s);
}  

NLP0::~NLP0()
{
}

void
NLP0::save_data_state(StateOut&s)
{
  _value.save_data_state(s);
  _dim.save_state(s);
}


RefSCDimension
NLP0::dimension()
{
  return _dim;
}

double
NLP0::value()
{
  return _value;
}

int
NLP0::do_value()
{
  return _value.compute();
}

int
NLP0::do_value(int f)
{
  return _value.compute(f);
}

void
NLP0::set_x(RefSCVector&v)
{
  _x.assign(v);
  obsolete();
}

RefSCVector
NLP0::get_x()
{
  return _x.copy();
}

void
NLP0::set_value(double e)
{
  _value.result_noupdate() = e;
  _value.computed() = 1;
}

void
NLP0::set_desired_value_accuracy(double a)
{
  _value.set_desired_accuracy(a);
}

double
NLP0::desired_value_accuracy()
{
  return _value.desired_accuracy();
}

double
NLP0::actual_value_accuracy()
{
  return _value.actual_accuracy();
}

void
NLP0::print(SCostream&o)
{
  o.indent(); o << "value_accuracy = " << desired_value_accuracy() << endl;
}

///////////////////////////////////////////////////////////////////////////

SavableState_REF_def(NLP1);

#define CLASSNAME NLP1
#define PARENTS public NLP0
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
NLP1::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = NLP0::_castdown(cd);
  return do_castdowns(casts,cd);
}

NLP1::NLP1(RefSCDimension&dim):
  NLP0(dim),
  _gradient(this)
{
  RefSCVector g(dim);
  _gradient.result_noupdate() = g;
  _gradient.set_desired_accuracy(DBL_EPSILON);
}

NLP1::NLP1(KeyVal&kv):
  NLP0(kv),
  _gradient(this)
{
  RefSCVector gradient(_dim);

  _gradient.set_desired_accuracy(kv.doublevalue("gradient_accuracy"));
  if (_gradient.desired_accuracy() < DBL_EPSILON)
    _gradient.set_desired_accuracy(DBL_EPSILON);

  _gradient = gradient;
}

NLP1::NLP1(StateIn&s):
  SavableState(s,class_desc_),
  NLP0(s),
  _gradient(s,this)
{
}

NLP1::~NLP1()
{
}

void
NLP1::save_data_state(StateOut&s)
{
  NLP0::save_data_state(s);
  _gradient.save_data_state(s);
}

RefSCVector
NLP1::gradient()
{
  return _gradient;
}

int
NLP1::do_gradient()
{
  return _gradient.compute();
}

int
NLP1::do_gradient(int f)
{
  return _gradient.compute(f);
}

void
NLP1::set_gradient(RefSCVector&g)
{
  _gradient.result_noupdate() = g;
  _gradient.computed() = 1;
}

void
NLP1::set_desired_gradient_accuracy(double a)
{
  _gradient.set_desired_accuracy(a);
}

double
NLP1::actual_gradient_accuracy()
{
  return _gradient.actual_accuracy();
}

double
NLP1::desired_gradient_accuracy()
{
  return _gradient.desired_accuracy();
}

void
NLP1::print(SCostream&o)
{
  NLP0::print(o);
  o.indent(); o << "gradient.accuracy = "
                << desired_gradient_accuracy()
                << endl;
}

///////////////////////////////////////////////////////////////////////////

SavableState_REF_def(NLP2);

#define CLASSNAME NLP2
#define PARENTS public NLP1
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
NLP2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = NLP1::_castdown(cd);
  return do_castdowns(casts,cd);
}

NLP2::NLP2(RefSCDimension&dim):
  NLP1(dim),
  _hessian(this)
{
  RefSymmSCMatrix h(dim);
  _hessian.result_noupdate() = h;
  _hessian.set_desired_accuracy(DBL_EPSILON);
}

NLP2::NLP2(KeyVal&kv):
  NLP1(kv),
  _hessian(this)
{
  RefSymmSCMatrix hessian(_dim);

  _hessian.set_desired_accuracy(kv.doublevalue("hessian_accuracy"));
  if (_hessian.desired_accuracy() < DBL_EPSILON)
    _hessian.set_desired_accuracy(DBL_EPSILON);

  _hessian = hessian;
}

NLP2::NLP2(StateIn&s):
  SavableState(s,class_desc_),
  NLP1(s),
  _hessian(s,this)
{
}

NLP2::~NLP2()
{
}

void
NLP2::save_data_state(StateOut&s)
{
  NLP1::save_data_state(s);
  _hessian.save_data_state(s);
}

RefSymmSCMatrix
NLP2::hessian()
{
  return _hessian;
}

int
NLP2::do_hessian()
{
  return _hessian.compute();
}

int
NLP2::do_hessian(int f)
{
  return _hessian.compute(f);
}

void
NLP2::set_hessian(RefSymmSCMatrix&h)
{
  _hessian.result_noupdate() = h;
  _hessian.computed() = 1;
}

// the default guess hessian is the unit diagonal
void
NLP2::guess_hessian(RefSymmSCMatrix&hessian)
{
  RefSCElementOp op(new SCElementShiftDiagonal(1.0));
  hessian.assign(0.0);
  hessian.element_op(op);
}

void
NLP2::set_desired_hessian_accuracy(double a)
{
  _hessian.set_desired_accuracy(a);
}

double
NLP2::desired_hessian_accuracy()
{
  return _hessian.desired_accuracy();
}

double
NLP2::actual_hessian_accuracy()
{
  return _hessian.actual_accuracy();
}

void
NLP2::print(SCostream&o)
{
  NLP1::print(o);
  o.indent(); o << "hessian_accuracy = " << desired_hessian_accuracy() << endl;
}

///////////////////////////////////////////////////////////////////////////
