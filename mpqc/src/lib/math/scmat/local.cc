
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// LocalSCMatrixKit member functions

#define CLASSNAME LocalSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LocalSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCMatrixKit::LocalSCMatrixKit()
{
}

LocalSCMatrixKit::LocalSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
}

LocalSCMatrixKit::LocalSCMatrixKit(StateIn&s):
  SCMatrixKit(s)
{
}

LocalSCMatrixKit::~LocalSCMatrixKit()
{
}

void
LocalSCMatrixKit::save_data_state(StateOut&s)
{
  SCMatrixKit::save_data_state(s);
}

SCDimension*
LocalSCMatrixKit::dimension(int n, const char* name)
{
  return new LocalSCDimension(n, name);
}

SCMatrixKit*
SCMatrixKit::default_matrixkit()
{
  return new LocalSCMatrixKit;
}

SCMatrix*
LocalSCMatrixKit::restore_matrix(StateIn& s,
                                const RefSCDimension& d1,
                                const RefSCDimension& d2)
{
  abort();
}

SymmSCMatrix*
LocalSCMatrixKit::restore_symmmatrix(StateIn& s, const RefSCDimension& d)
{
  abort();
}

DiagSCMatrix*
LocalSCMatrixKit::restore_diagmatrix(StateIn& s, const RefSCDimension& d)
{
  abort();
}

SCVector*
LocalSCMatrixKit::restore_vector(StateIn& s, const RefSCDimension& d)
{
  abort();
}

/////////////////////////////////////////////////////////////////////////////
// LocalSCDimension member functions

#define CLASSNAME LocalSCDimension
#define PARENTS public SCDimension
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LocalSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCDimension::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCDimension::LocalSCDimension(int n, const char* name):
  SCDimension(name),
  n_(n)
{
}

LocalSCDimension::LocalSCDimension(StateIn&s):
  SCDimension(s)
{
  s.get(n_);
}

LocalSCDimension::LocalSCDimension(const RefKeyVal&keyval)
{
  n_ = keyval->intvalue("n");
}

void
LocalSCDimension::save_data_state(StateOut&s)
{
  SCDimension::save_data_state(s);
  s.put(n_);
}

LocalSCDimension::~LocalSCDimension()
{
}

int
LocalSCDimension::equiv(SCDimension*a) const
{
  // for now require another LocalSCDimension
  LocalSCDimension *la = LocalSCDimension::castdown(a);

  if (!la || n_ != la->n_)
    return 0;

  return 1;
}

int
LocalSCDimension::n()
{
  return n_;
}
SCMatrix*
LocalSCDimension::create_matrix(SCDimension*a)
{
  LocalSCDimension*coldim
    = LocalSCDimension::require_castdown(a,"LocalSCDimension::create_matrix");
  return new LocalSCMatrix(this,coldim);
}
SymmSCMatrix*
LocalSCDimension::create_symmmatrix()
{
  return new LocalSymmSCMatrix(this);
}
DiagSCMatrix*
LocalSCDimension::create_diagmatrix()
{
  return new LocalDiagSCMatrix(this);
}
SCVector*
LocalSCDimension::create_vector()
{
  return new LocalSCVector(this);
}

SavableState_REF_def(LocalSCDimension);
