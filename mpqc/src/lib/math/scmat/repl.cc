
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
  grp_ = keyval->describedclassvalue("MessageGrp");
  if (grp_.null()) grp_ = MessageGrp::get_default_messagegrp();
}

ReplSCMatrixKit::~ReplSCMatrixKit()
{
}

SCDimension*
ReplSCMatrixKit::dimension(int n, const char* name)
{
  return new ReplSCDimension(n, grp_, name);
}

/////////////////////////////////////////////////////////////////////////////
// ReplSCDimension member functions

#define CLASSNAME ReplSCDimension
#define PARENTS public SCDimension
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ReplSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCDimension::_castdown(cd);
  return do_castdowns(casts,cd);
}

ReplSCDimension::ReplSCDimension(int n, const RefMessageGrp& grp,
                                 const char* name):
  SCDimension(name),
  grp_(grp),
  n_(n)
{
  int nproc = messagegrp()->n();
  int nperproc = n_/nproc;
  int nleftover = n_%nproc;

  if (nperproc) nblocks_ = nproc;
  else nblocks_ = nleftover;

  blocks_ = new int[nblocks_+1];
  int i;
  blocks_[0] = 0;
  for (i=1; i<=nproc; i++) {
      if (nperproc) blocks_[i] = blocks_[i-1] + nperproc;
      if (i<=nleftover) blocks_[i]++;
    }
}

ReplSCDimension::~ReplSCDimension()
{
  delete[] blocks_;
}

int
ReplSCDimension::equiv(SCDimension *a) const
{
  ReplSCDimension *ra = ReplSCDimension::castdown(a);

  if (!ra)
    return 0;

  if (n_ != ra->n_ || nblocks_ != ra->nblocks_)
    return 0;

  for (int i=0; i < nblocks_; i++)
    if (blocks_[i] != ra->blocks_[i])
      return 0;

  return grp_==ra->grp_;
}

int
ReplSCDimension::n()
{
  return n_;
}
SCMatrix*
ReplSCDimension::create_matrix(SCDimension*a)
{
  ReplSCDimension*coldim
    = ReplSCDimension::require_castdown(a,"ReplSCDimension::create_matrix");
  return new ReplSCMatrix(this,coldim);
}
SymmSCMatrix*
ReplSCDimension::create_symmmatrix()
{
  return new ReplSymmSCMatrix(this);
}
DiagSCMatrix*
ReplSCDimension::create_diagmatrix()
{
  return new ReplDiagSCMatrix(this);
}
SCVector*
ReplSCDimension::create_vector()
{
  return new ReplSCVector(this);
}

SavableState_REF_def(ReplSCDimension);

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrixKit member functions

ReplSCMatrixListSubblockIter::ReplSCMatrixListSubblockIter(
    Access access,
    const RefSCMatrixBlockList &list,
    const RefMessageGrp &grp,
    double *data,
    int ndata
    ):
  SCMatrixListSubblockIter(access_, list),
  grp_(grp),
  data_(data),
  ndata_(ndata)
{
  if (access == Write || (access == Read && grp->me() != 0)) {
      for (int i=0; i<ndata; i++) data[i] = 0.0;
    }
}

ReplSCMatrixListSubblockIter::~ReplSCMatrixListSubblockIter()
{
  if (access() == Write || access() == Accum) {
      grp_->sum(data_,ndata_);
    }
}
