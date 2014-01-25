//
// dist.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <iostream>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// DistSCMatrixKit member functions

static ClassDesc DistSCMatrixKit_cd(
  typeid(DistSCMatrixKit),"DistSCMatrixKit",1,"public SCMatrixKit",
  0, create<DistSCMatrixKit>, 0);

DistSCMatrixKit::DistSCMatrixKit(const Ref<MessageGrp> &grp)
{
  // if grp is nonnull, then reset grp_ (it gets set to the default in the
  // default SCMatrixKit constructor
  if (grp.nonnull())
    grp_ = grp;
}

DistSCMatrixKit::DistSCMatrixKit(const Ref<KeyVal>& keyval):
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
    const Ref<SCMatrixBlockList> &list,
    const Ref<MessageGrp> &grp
    ):
  SCMatrixListSubblockIter(access, list->deepcopy()),
  grp_(grp),
  out_(grp),
  in_(grp),
  step_(0),
  locallist_(list)
{
  if (access == Write) {
      ExEnv::errn() << indent
           << "DistSCMatrixListSubblockIter: write access not allowed"
           << endl;
      abort();
    }

  if (grp->n() == 1) return;

  out_.target(grp->me() == grp->n()-1 ? 0: grp->me()+1);
  in_.source(grp->me() == 0 ? grp->n()-1: grp->me()-1);

  out_.copy_references();
}

void
DistSCMatrixListSubblockIter::begin()
{
  if (step_ == grp_->n()) step_ = 0;
  else if (step_ != 0) {
      ExEnv::errn() << indent << "DistSCMatrixListSubblockIter::begin(): "
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
  SavableState::save_state(list_.pointer(), out_);
  out_.flush();
  list_ << SavableState::restore_state(in_);
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
              ExEnv::errn() << indent
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
