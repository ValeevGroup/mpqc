
#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// DistSCMatrixKit member functions

#define CLASSNAME DistSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
DistSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

DistSCMatrixKit::DistSCMatrixKit(const RefMessageGrp &grp)
{
  // if grp is nonnull, then reset grp_ (it gets set to the default in the
  // default SCMatrixKit constructor
  if (grp.nonnull())
    grp_ = grp;
}

DistSCMatrixKit::DistSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
}

DistSCMatrixKit::~DistSCMatrixKit()
{
}

SCMatrix*
DistSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  return new DistSCMatrix(d1,d2,this);
}

SymmSCMatrix*
DistSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  return new DistSymmSCMatrix(d,this);
}

DiagSCMatrix*
DistSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  return new DistDiagSCMatrix(d,this);
}

SCVector*
DistSCMatrixKit::vector(const RefSCDimension&d)
{
  return new DistSCVector(d,this);

}

/////////////////////////////////////////////////////////////////////////////
// DistSCMatrixKit member functions

DistSCMatrixListSubblockIter::DistSCMatrixListSubblockIter(
    Access access,
    const RefSCMatrixBlockList &list,
    const RefMessageGrp &grp
    ):
  SCMatrixListSubblockIter(access, list->deepcopy()),
  locallist_(list),
  grp_(grp),
  step_(0),
  out_(grp),
  in_(grp)
{
  if (access == Write) {
      cerr << indent
           << "DistSCMatrixListSubblockIter: write access not allowed"
           << endl;
      abort();
    }

  if (grp->n() == 1) return;

  out_.target(grp->me() == grp->n()-1 ? 0: grp->me()+1);
  in_.source(grp->me() == 0 ? grp->n()-1: grp->me()-1);

  out_.copy_references();
  in_.copy_references();
}

void
DistSCMatrixListSubblockIter::begin()
{
  if (step_ == grp_->n()) step_ = 0;
  else if (step_ != 0) {
      cerr << indent << "DistSCMatrixListSubblockIter::begin(): "
           << "step != 0: tried to begin in middle of iteration"
           << endl;
      abort();
    }
  SCMatrixListSubblockIter::begin();
  maybe_advance_list();
}

void
DistSCMatrixListSubblockIter::maybe_advance_list()
{
  while (!ready() && grp_->n() > 1 && step_ < grp_->n() - 1) {
      advance_list();
    }
}

void
DistSCMatrixListSubblockIter::advance_list()
{
  list_.save_state(out_);
  out_.flush();
  list_.restore_state(in_);
  SCMatrixListSubblockIter::begin();
  step_++;
}

void
DistSCMatrixListSubblockIter::next()
{
  SCMatrixListSubblockIter::next();
  maybe_advance_list();
}

DistSCMatrixListSubblockIter::~DistSCMatrixListSubblockIter()
{
  if (access_ == Accum) {
      while (step_%grp_->n() != 0) {
          advance_list();
        }
      SCMatrixBlockListIter i1, i2;
      for (i1=list_->begin(),i2=locallist_->begin();
           i1!=list_->end() && i2!=locallist_->end();
           i1++,i2++) {
          int n = i1.block()->ndat();
          if (n != i2.block()->ndat()) {
              cerr << indent
                   << "DistSCMatrixListSubblockIter: block mismatch: "
                   << "internal error" << endl;
              abort();
            }
          double *dat1 = i1.block()->dat();
          double *dat2 = i2.block()->dat();
          for (int i=0; i<n; i++) {
              dat2[i] += dat1[i];
            }
        }
    }
}
