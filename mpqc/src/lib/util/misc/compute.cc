
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/compute.h>
#include <util/state/state.h>

#ifdef __GNUC__
template class Result<int>;
template class Result<double>;
template class NCAccResult<double>;
#endif

ARRAY_def(ResultBaseP);
SET_def(ResultBaseP);

Compute::Compute()
{
}

Compute::~Compute()
{
}

void
Compute::add(ResultBase*r)
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

////////////////////////////////////////////////////////////////////////

ResultBase::ResultBase(Compute*c):
  _compute(0),_computed(0),_c(c)
{
  c->add(this);
}

void
ResultBase::update() {
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

ResultBase::~ResultBase()
{
}

ResultBase::ResultBase(StateIn&s,Compute*c):
  _c(c)
{
  s.get(_compute);
  s.get(_computed);

  c->add(this);
}

ResultBase::ResultBase(const ResultBase&r, Compute*c) :
  _c(c)
{
  _compute=r._compute;
  _computed=r._computed;
  
  c->add(this);
}

void
ResultBase::save_data_state(StateOut&s)
{
  s.put(_compute);
  s.put(_computed);
}

ResultBase&
ResultBase::operator=(const ResultBase&r)
{
  _compute=r._compute;
  _computed=r._computed;
  return *this;
}

/////////////////////////////////////////////////////////////////////////

AccResultBase::AccResultBase(Compute*c):
  ResultBase(c),
  _actual_accuracy(0.0),
  _desired_accuracy(0.01)
{
}

AccResultBase::~AccResultBase()
{
}

double
AccResultBase::actual_accuracy() const
{
  return _actual_accuracy;
}

double
AccResultBase::desired_accuracy() const
{
  return _desired_accuracy;
}

void
AccResultBase::set_desired_accuracy(double a)
{
  _desired_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      computed() = 0;
    }
}

void
AccResultBase::set_actual_accuracy(double a)
{
  _actual_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      fprintf(stderr,"AccResult: setting actual accuracy greater than"
              " desired accuracy\n");
      abort();
    }
  computed() = 1;
}

AccResultBase::AccResultBase(StateIn&s,Compute*c):
  ResultBase(s,c)
{
  s.get(_actual_accuracy);
  s.get(_desired_accuracy);
}

AccResultBase::AccResultBase(const AccResultBase&a, Compute*c) :
  ResultBase(a,c)
{
  _actual_accuracy=a._actual_accuracy;
  _desired_accuracy=a._desired_accuracy;
}

void
AccResultBase::save_data_state(StateOut&s)
{
  ResultBase::save_data_state(s);
  s.put(_actual_accuracy);
  s.put(_desired_accuracy);
}

AccResultBase&
AccResultBase::operator=(const AccResultBase&a)
{
  ResultBase::operator=(a);
  _actual_accuracy=a._actual_accuracy;
  _desired_accuracy=a._desired_accuracy;
  return *this;
}
