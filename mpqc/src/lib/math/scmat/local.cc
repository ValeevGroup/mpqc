
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

LocalSCMatrixKit::~LocalSCMatrixKit()
{
}

SCMatrix*
LocalSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return new LocalSCMatrix(d1,d2,this);
}

SymmSCMatrix*
LocalSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return new LocalSymmSCMatrix(d,this);
}

DiagSCMatrix*
LocalSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return new LocalDiagSCMatrix(d,this);
}

SCVector*
LocalSCMatrixKit::vector(const RefSCDimension&d)
{
  return new LocalSCVector(d,this);

}

///////////////////////////////////////////////////////////////////////
// The static SCMatrixKit members.

static RefSCMatrixKit defaultmatrixkit;

SCMatrixKit*
SCMatrixKit::default_matrixkit()
{
  if (defaultmatrixkit.null()) defaultmatrixkit = new LocalSCMatrixKit;
  return defaultmatrixkit.pointer();
}

void
SCMatrixKit::set_default_matrixkit(const RefSCMatrixKit &k)
{
  defaultmatrixkit = k;
}
