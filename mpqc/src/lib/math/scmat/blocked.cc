
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <math.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedSCMatrixKit member functions

#define CLASSNAME BlockedSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCMatrixKit::BlockedSCMatrixKit()
{
}

BlockedSCMatrixKit::BlockedSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
}

BlockedSCMatrixKit::BlockedSCMatrixKit(StateIn&s):
  SCMatrixKit(s)
{
}

BlockedSCMatrixKit::~BlockedSCMatrixKit()
{
}

void
BlockedSCMatrixKit::save_data_state(StateOut&s)
{
  SCMatrixKit::save_data_state(s);
}

SCDimension*
BlockedSCMatrixKit::dimension(int n, const char* name)
{
  return new BlockedSCDimension(default_matrixkit(), n, name);
}

SCDimension*
BlockedSCMatrixKit::dimension(int n, int *nelem, const char* name)
{
  return new BlockedSCDimension(default_matrixkit(), n, nelem, name);
}

SCDimension*
BlockedSCMatrixKit::dimension(const RefSCMatrixKit& mk,
                              int n, const char* name)
{
  return new BlockedSCDimension(mk, n, name);
}

SCDimension*
BlockedSCMatrixKit::dimension(const RefSCMatrixKit& mk,
                              int n, int *nelem, const char* name)
{
  return new BlockedSCDimension(mk, n, nelem, name);
}


SCMatrix*
BlockedSCMatrixKit::restore_matrix(StateIn& s,
                                const RefSCDimension& d1,
                                const RefSCDimension& d2)
{
  abort();
}

SymmSCMatrix*
BlockedSCMatrixKit::restore_symmmatrix(StateIn& s, const RefSCDimension& d)
{
  abort();
}

DiagSCMatrix*
BlockedSCMatrixKit::restore_diagmatrix(StateIn& s, const RefSCDimension& d)
{
  abort();
}

SCVector*
BlockedSCMatrixKit::restore_vector(StateIn& s, const RefSCDimension& d)
{
  abort();
}

SavableState_REF_def(BlockedSCMatrixKit);

/////////////////////////////////////////////////////////////////////////////
// BlockedSCDimension member functions

#define CLASSNAME BlockedSCDimension
#define PARENTS public SCDimension
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BlockedSCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCDimension::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSCDimension::init(const RefSCMatrixKit& mk, int *nb, const char *name)
{
  if (dims_) {
    delete[] dims_;
    dims_=0;
  }
  if (first_) {
    delete[] first_;
    first_=0;
  }
  if (last_) {
    delete[] last_;
    last_=0;
  }

  if (!nblocks_)
    return;

  if (!nb) {
    if (nblocks_ > 1) {
      fprintf(stderr,"BlockedSCDimension::init: no block info\n");
      abort();
    }

    int *newnb = new int[1];
    newnb[0] = n_;
    init(mk,newnb,name);
    delete[] newnb;
    return;
  }

  dims_ = new RefSCDimension[nblocks_];
  first_ = new int[nblocks_];
  last_ = new int[nblocks_];

  n_ = nb[0];
  first_[0] = 0;
  last_[0] = nb[0];
  dims_[0] = mk->dimension(nb[0], name);
  
  for (int i=1; i < nblocks_; i++) {
    n_ += nb[0];
    first_[i] = last_[i-1];
    last_[i] = first_[i]+nb[i];
    dims_[i] = mk->dimension(nb[i], name);
  }
}

BlockedSCDimension::BlockedSCDimension(const RefSCMatrixKit& mk,
                                       int n, int *nelem, const char* name) :
  SCDimension(name),
  nblocks_(n),
  dims_(0),
  first_(0),
  last_(0)
{
  init(mk, nelem, name);
}

BlockedSCDimension::BlockedSCDimension(const RefSCMatrixKit& mk,
                                       int n, const char* name) :
  SCDimension(name),
  n_(n),
  dims_(0),
  first_(0),
  last_(0)
{
  nblocks_=1;
  init(mk, 0, name);
}

BlockedSCDimension::BlockedSCDimension(StateIn&s):
  SCDimension(s),
  first_(0),
  last_(0)
{
  s.get(n_);
  s.get(nblocks_);
  s.get(first_);
  s.get(last_);

  dims_ = new RefSCDimension[nblocks_];
  for (int i=0; i < nblocks_; i++)
    dims_[i].restore_state(s);
}

BlockedSCDimension::BlockedSCDimension(const RefKeyVal&keyval)
{
  abort();
}

void
BlockedSCDimension::save_data_state(StateOut&s)
{
  SCDimension::save_data_state(s);
  s.put(n_);
  s.put(nblocks_);
  s.put(first_,nblocks_);
  s.put(last_,nblocks_);
  for (int i=0; i < nblocks_; i++)
    dims_[i].save_state(s);
}

BlockedSCDimension::~BlockedSCDimension()
{
  n_=nblocks_=0;
  init(0,0,0);
}

int
BlockedSCDimension::equiv(SCDimension *a) const
{
  BlockedSCDimension *ba = BlockedSCDimension::castdown(a);

  if (!ba || ba->n_ != n_ || ba->nblocks_ != nblocks_)
    return 0;

  for (int i=0; i < nblocks_; i++) {
    if (first_[i] != ba->first_[i] ||
        last_[i] != ba->last_[i] ||
        !dims_[i]->equiv(ba->dims_[i].pointer()))
      return 0;
  }

  return 1;
}

int
BlockedSCDimension::n()
{
  return n_;
}

int
BlockedSCDimension::n(int i) const
{
  return dims_[i]->n();
}

int
BlockedSCDimension::nblocks() const
{
  return nblocks_;
}

int
BlockedSCDimension::first(int i) const
{
  return first_[i];
}

int
BlockedSCDimension::last(int i) const
{
  return last_[i];
}

int
BlockedSCDimension::block(int i) const
{
  for (int j=0; j < nblocks_; j++)
    if (i < last_[j])
      return j;
}

RefSCDimension
BlockedSCDimension::dim(int i) const
{
  return dims_[i];
}

SCMatrix*
BlockedSCDimension::create_matrix(SCDimension*a)
{
  BlockedSCDimension*coldim =
    BlockedSCDimension::require_castdown(a,
                                         "BlockedSCDimension::create_matrix");

  return new BlockedSCMatrix(this,coldim);
}

SymmSCMatrix*
BlockedSCDimension::create_symmmatrix()
{
  return new BlockedSymmSCMatrix(this);
}

DiagSCMatrix*
BlockedSCDimension::create_diagmatrix()
{
  return new BlockedDiagSCMatrix(this);
}

SCVector*
BlockedSCDimension::create_vector()
{
  return new BlockedSCVector(this);
}

SavableState_REF_def(BlockedSCDimension);
