
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrixKit member functions

#define CLASSNAME ReplSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

ReplSCMatrixKit::ReplSCMatrixKit()
{
  grp_ = MessageGrp::get_default_messagegrp();
}

ReplSCMatrixKit::ReplSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
  grp_ = keyval->describedclassvalue("messagegrp");
  if (grp_.null()) grp_ = MessageGrp::get_default_messagegrp();
}

ReplSCMatrixKit::~ReplSCMatrixKit()
{
}

SCMatrix*
ReplSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return new ReplSCMatrix(d1,d2,this);
}

SymmSCMatrix*
ReplSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return new ReplSymmSCMatrix(d,this);
}

DiagSCMatrix*
ReplSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return new ReplDiagSCMatrix(d,this);
}

SCVector*
ReplSCMatrixKit::vector(const RefSCDimension&d)
{
  return new ReplSCVector(d,this);

}

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrixKit member functions

ReplSCMatrixListSubblockIter::ReplSCMatrixListSubblockIter(
    Access access,
    const RefSCMatrixBlockList &list,
    const RefMessageGrp &grp,
    double *data,
    int ndata
    ):
  SCMatrixListSubblockIter(access, list),
  grp_(grp),
  data_(data),
  ndata_(ndata)
{
  if (access == Write) {
      for (int i=0; i<ndata; i++) data[i] = 0.0;
    }
}

ReplSCMatrixListSubblockIter::~ReplSCMatrixListSubblockIter()
{
  if (access() == Write || access() == Accum) {
      grp_->sum(data_,ndata_);
    }
}
