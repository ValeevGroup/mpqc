
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// DistSCMatrixKit member functions

#define CLASSNAME DistSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
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
  if (grp.null()) grp_ = MessageGrp::get_default_messagegrp();
  else grp_ = grp;
}

DistSCMatrixKit::DistSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
  grp_ = keyval->describedclassvalue("MessageGrp");
  if (grp_.null()) grp_ = MessageGrp::get_default_messagegrp();
}

DistSCMatrixKit::~DistSCMatrixKit()
{
}

SCDimension*
DistSCMatrixKit::dimension(int n, const char* name)
{
  return new DistSCDimension(n, this, name);
}

SCDimension*
DistSCMatrixKit::dimension(int n, int nblocks,
                           const int *blocksizes,
                           const char* name)
{
  return new DistSCDimension(n, this, nblocks, blocksizes, name);
}

/////////////////////////////////////////////////////////////////////////////
// DistSCDimension member functions

#define CLASSNAME DistSCDimension
#define PARENTS public SCDimension
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DistSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCDimension::_castdown(cd);
  return do_castdowns(casts,cd);
}

DistSCDimension::DistSCDimension(int n, const RefDistSCMatrixKit& kit,
                                 const char* name):
  SCDimension(name),
  kit_(kit),
  n_(n)
{
  blocks_ = new SCBlockInfo(n, (int)sqrt(2.0*kit->messagegrp()->n()));
}

DistSCDimension::DistSCDimension(int n, const RefDistSCMatrixKit& kit,
                                 int nblocks, const int *blocksizes,
                                 const char* name):
  SCDimension(name),
  kit_(kit),
  n_(n)
{
  blocks_ = new SCBlockInfo(n, nblocks, blocksizes);
}

DistSCDimension::~DistSCDimension()
{
}

int
DistSCDimension::equiv(SCDimension *a) const
{
  DistSCDimension *ra = DistSCDimension::castdown(a);

  if (n_ != ra->n_) return 0;

  if (!blocks_->equiv(ra->blocks_)) return 0;

  return kit_==ra->kit_;
}

int
DistSCDimension::n()
{
  return n_;
}
SCMatrix*
DistSCDimension::create_matrix(SCDimension*a)
{
  DistSCDimension*coldim
    = DistSCDimension::require_castdown(a,"DistSCDimension::create_matrix");
  return new DistSCMatrix(this,coldim);
}
SymmSCMatrix*
DistSCDimension::create_symmmatrix()
{
  return new DistSymmSCMatrix(this);
}
DiagSCMatrix*
DistSCDimension::create_diagmatrix()
{
  return new DistDiagSCMatrix(this);
}
SCVector*
DistSCDimension::create_vector()
{
  return new DistSCVector(this);
}

SavableState_REF_def(DistSCDimension);

/////////////////////////////////////////////////////////////////////////////
// SCBlockInfo member functions

SCBlockInfo::SCBlockInfo(int n, int nblocks, const int *blocksizes)
{
  n_ = n;
  nblocks_ = nblocks;

  if (n_ == 0) nblocks_ = 0;

  if (n_ != 0 && nblocks_ == 0) {
      nblocks_ = 1;
      start_ = new int[1];
      size_ = new int[1];
      start_[0] = 0;
      size_[0] = n;
    }
  else if (nblocks_ == 0) {
      start_ = 0;
      size_ = 0;
    }
  else {
      int i;
      start_ = new int[nblocks_];
      size_ = new int[nblocks_];
      start_[0] = 0;
      if (blocksizes) {
          for (i=0; i<nblocks_; i++) {
              size_[i] = blocksizes[i];
              if (i+1<nblocks_) start_[i+1] = start_[i] + size_[i];
            }
        }
      else {
          int nper = n/nblocks;
          int nleft = n%nblocks;
          for (i=0; i<nblocks_; i++) {
              size_[i] = nper;
              if (i<nleft) size_[i]++;
              if (i+1<nblocks_) start_[i+1] = start_[i] + size_[i];
            }
        }
    }
}

int
SCBlockInfo::equiv(SCBlockInfo *bi)
{
  if (bi == 0) return 0;
  if (n_ != bi->n_) return 0;
  if (nblocks_ != bi->nblocks_) return 0;
  int i;
  for (i=0; i<nblocks_; i++) {
      if (start_[i] != bi->start_[i]) return 0;
    }
  return 1;
}

void
SCBlockInfo::elem_to_block(int elem, int &block, int &offset)
{
  for (int i=0; i<nblocks_; i++) {
      if (start_[i] <= elem) {
          block = i;
          offset = elem-start_[i];
          return;
        }
    }
  cerr << "SCBlockInfo::elem_to_block: couldn't find block" << endl;
}

SCBlockInfo::~SCBlockInfo()
{
  delete[] start_;
  delete[] size_;
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
      cerr << "DistSCMatrixListSubblockIter: write access not allowed"
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
      cerr << "DistSCMatrixListSubblockIter::begin(): "
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
  while (!ready() && grp_->n() > 1 && step_ < grp_->n()) {
      list_.save_state(out_);
      out_.flush();
      list_.restore_state(in_);
      SCMatrixListSubblockIter::begin();
      step_++;
    }
}

void
DistSCMatrixListSubblockIter::next()
{
  SCMatrixListSubblockIter::next();
  maybe_advance_list();
}

DistSCMatrixListSubblockIter::~DistSCMatrixListSubblockIter()
{
  if (step_%grp_->n() != 0) {
      cerr << "DistSCMatrixListSubblockIter: DTOR: "
           << "step != 0: tried to end in middle of iteration"
           << endl;
      abort();
    }
  if (access_ == Accum) {
      SCMatrixBlockListIter i1, i2;
      for (i1=list_->begin(),i2=locallist_->begin();
           i1!=list_->end() && i2!=locallist_->end();
           i1++,i2++) {
          int n = i1.block()->ndat();
          if (n != i2.block()->ndat()) {
              cerr << "DistSCMatrixListSubblockIter: block mismatch: "
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
