
#include "compute.h"
#include <util/state/state.h>

Result_def_nc(int);
Result_def_nc(float);
Result_def_nc(double);
SSAccResult_def_nc(double);

ARRAY_def(ResultP);
SET_def(ResultP);

Compute::Compute()
{
}

Compute::~Compute()
{
}

void
Compute::add(Result*r)
{
  _results.add(r);
}

void
Compute::obsolete()
{
  // go thru all of the results and mark them as obsolete
  for (Pix i = _results.first(); i; _results.next(i)) {
      _results(i)->computed() = 0;
    }
}

void
Result::update() {
  if (!computed()) {
      int oldcompute = compute(1);
      _c->compute();
      compute() = oldcompute;
      if (!computed()) {
          fprintf(stderr,"Result::compute: nothing was computed\n");
          abort();
        }
    }
}

Result::~Result()
{
}

Result::Result(StateIn&s,Compute*c):
  _c(c)
{
  s.get(_compute);
  s.get(_computed);

  c->add(this);
}

void
Result::save_data_state(StateOut&s)
{
  Result::save_data_state(s);
  s.put(_compute);
  s.put(_computed);
}

AccResult::AccResult(Compute*c):
  Result(c),
  _actual_accuracy(0.0),
  _desired_accuracy(0.01)
{
}

AccResult::~AccResult()
{
}

double
AccResult::actual_accuracy()
{
  return _actual_accuracy;
}

double
AccResult::desired_accuracy()
{
  return _desired_accuracy;
}

void
AccResult::set_desired_accuracy(double a)
{
  _desired_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      computed() = 0;
    }
}

void
AccResult::set_actual_accuracy(double a)
{
  _actual_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      fprintf(stderr,"AccResult: setting actual accuracy greater than"
              " desired accuracy\n");
      abort();
    }
}

AccResult::AccResult(StateIn&s,Compute*c):
  Result(s,c)
{
  s.get(_actual_accuracy);
  s.get(_desired_accuracy);
}

void
AccResult::save_data_state(StateOut&s)
{
  Result::save_data_state(s);
  s.put(_actual_accuracy);
  s.put(_desired_accuracy);
}
