
static char rcsid[] = "$Header$";

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
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

NLP0::NLP0(RefSCDimension&dim):
  _value(this),
  _dim(dim),
  _x(dim),
  _acc_value(DBL_EPSILON),
  _maxacc_value(DBL_EPSILON)
{
}

NLP0::NLP0(KeyVal&kv):
  _value(this),
  _acc_value(DBL_EPSILON),
  _maxacc_value(DBL_EPSILON)
{
  _dim = kv.describedclassvalue("dimension");
  RefSCVector x(_dim);
  _x = x;
}

NLP0::NLP0(StateIn&s):
  SavableState(s,class_desc_),
  _value(this)
{
  _dim.restore_state(s);
  s.get(_value.computed());
  s.get(_value.compute());
  s.get(_acc_value);
  s.get(_maxacc_value);
  if (_value.computed()) {
      s.get(_value.result_noupdate());
    }
}  

NLP0::~NLP0()
{
}

void
NLP0::save_data_state(StateOut&s)
{
  _dim.save_state(s);
  s.put(_value.computed());
  s.put(_value.compute());
  s.put(_acc_value);
  s.put(_maxacc_value);
  if (_value.computed()) {
      s.put(_value.result_noupdate());
    }
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
NLP0::set_x(int i,double val)
{
  _x(i) = val;
  obsolete();
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

///////////////////////////////////////////////////////////////////////////

SavableState_REF_def(NLP1);

#define CLASSNAME NLP1
#define PARENTS public NLP0
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
NLP1::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP0::_castdown(cd) };
  return do_castdowns(casts,cd);
}

NLP1::NLP1(RefSCDimension&dim):
  NLP0(dim),
  _gradient(this)
{
  RefSCMatrix g(dim);
  _gradient.result_noupdate() = g;
}

NLP1::NLP1(KeyVal&kv):
  NLP0(kv),
  _gradient(this)
{
  RefSCVector gradient(_dim);
  _gradient = gradient;
}

NLP1::NLP1(StateIn&s):
  SavableState(s,class_desc_),
  NLP0(s),
  _gradient(this)
{
  s.get(_gradient.computed());
  s.get(_gradient.compute());
  if (_gradient.computed()) {
      _gradient.result_noupdate().restore_state(s);
    }
}

NLP1::~NLP1()
{
}

void
NLP1::save_data_state(StateOut&s)
{
  NLP0::save_data_state(s);
  s.put(_gradient.computed());
  s.put(_gradient.compute());
  if (_gradient.computed()) {
      _gradient.result_noupdate().save_state(s);
    }
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

///////////////////////////////////////////////////////////////////////////

SavableState_REF_def(NLP2);

#define CLASSNAME NLP2
#define PARENTS public NLP1
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
NLP2::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { NLP1::_castdown(cd) };
  return do_castdowns(casts,cd);
}

NLP2::NLP2(RefSCDimension&dim):
  NLP1(dim),
  _hessian(this)
{
  RefSymmSCMatrix h(dim);
  _hessian.result_noupdate() = h;
}

NLP2::NLP2(KeyVal&kv):
  NLP1(kv),
  _hessian(this)
{
  RefSymmSCMatrix hessian(_dim);
  _hessian = hessian;
}

NLP2::NLP2(StateIn&s):
  SavableState(s,class_desc_),
  NLP1(s),
  _hessian(this)
{
  s.get(_hessian.computed());
  s.get(_hessian.compute());
  if (_hessian.computed()) {
      _hessian.result_noupdate().restore_state(s);
    }
}

NLP2::~NLP2()
{
}

void
NLP2::save_data_state(StateOut&s)
{
  NLP1::save_data_state(s);
  s.put(_hessian.computed());
  s.put(_hessian.compute());
  if (_hessian.computed()) {
      _hessian.result_noupdate().save_state(s);
    }
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

///////////////////////////////////////////////////////////////////////////
