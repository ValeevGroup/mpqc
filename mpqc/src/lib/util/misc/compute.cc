
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

ARRAY_def(ResultInfoP);
SET_def(ResultInfoP);

Compute::Compute()
{
}

Compute::~Compute()
{
}

void
Compute::add(ResultInfo*r)
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

ResultInfo::ResultInfo(Compute*c):
  _compute(0),_computed(0),_c(c)
{
  c->add(this);
}

void
ResultInfo::update() {
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

ResultInfo::~ResultInfo()
{
}

ResultInfo::ResultInfo(StateIn&s,Compute*c):
  _c(c)
{
  s.get(_compute);
  s.get(_computed);

  c->add(this);
}

ResultInfo::ResultInfo(const ResultInfo&r, Compute*c) :
  _c(c)
{
  _compute=r._compute;
  _computed=r._computed;
  
  c->add(this);
}

void
ResultInfo::save_data_state(StateOut&s)
{
  s.put(_compute);
  s.put(_computed);
}

ResultInfo&
ResultInfo::operator=(const ResultInfo&r)
{
  _compute=r._compute;
  _computed=r._computed;
  return *this;
}

/////////////////////////////////////////////////////////////////////////

AccResultInfo::AccResultInfo(Compute*c):
  ResultInfo(c),
  _actual_accuracy(0.0),
  _desired_accuracy(0.01)
{
}

AccResultInfo::~AccResultInfo()
{
}

double
AccResultInfo::actual_accuracy() const
{
  return _actual_accuracy;
}

double
AccResultInfo::desired_accuracy() const
{
  return _desired_accuracy;
}

void
AccResultInfo::set_desired_accuracy(double a)
{
  _desired_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      computed() = 0;
    }
}

void
AccResultInfo::set_actual_accuracy(double a)
{
  _actual_accuracy = a;
  if (_desired_accuracy < _actual_accuracy) {
      fprintf(stderr,"AccResult: setting actual accuracy greater than"
              " desired accuracy\n");
      abort();
    }
  computed() = 1;
}

AccResultInfo::AccResultInfo(StateIn&s,Compute*c):
  ResultInfo(s,c)
{
  s.get(_actual_accuracy);
  s.get(_desired_accuracy);
}

AccResultInfo::AccResultInfo(const AccResultInfo&a, Compute*c) :
  ResultInfo(a,c)
{
  _actual_accuracy=a._actual_accuracy;
  _desired_accuracy=a._desired_accuracy;
}

void
AccResultInfo::save_data_state(StateOut&s)
{
  ResultInfo::save_data_state(s);
  s.put(_actual_accuracy);
  s.put(_desired_accuracy);
}

AccResultInfo&
AccResultInfo::operator=(const AccResultInfo&a)
{
  ResultInfo::operator=(a);
  _actual_accuracy=a._actual_accuracy;
  _desired_accuracy=a._desired_accuracy;
  return *this;
}
