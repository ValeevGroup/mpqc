
#include "compute.h"

Result_def_nc(int);
Result_def_nc(float);
Result_def_nc(double);

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
