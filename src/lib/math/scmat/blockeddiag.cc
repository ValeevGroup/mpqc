//
// blockeddiag.cc
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

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// BlockedDiagSCMatrix member functions

static ClassDesc BlockedDiagSCMatrix_cd(
  typeid(BlockedDiagSCMatrix),"BlockedDiagSCMatrix",1,"public DiagSCMatrix",
  0, 0, 0);

void
BlockedDiagSCMatrix::resize(SCDimension *a)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }

  d = a;

  mats_ = new RefDiagSCMatrix[d->blocks()->nblock()];
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (d->blocks()->size(i))
      mats_[i] = subkit->diagmatrix(d->blocks()->subdim(i));
}

BlockedDiagSCMatrix::BlockedDiagSCMatrix(const RefSCDimension&a,
                                         BlockedSCMatrixKit*k):
  DiagSCMatrix(a,k),
  subkit(k->subkit()),
  mats_(0)
{
  resize(a);
}

BlockedDiagSCMatrix::~BlockedDiagSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
}

double
BlockedDiagSCMatrix::get_element(int i) const
{
  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  return mats_[bi]->get_element(bo);
}

void
BlockedDiagSCMatrix::set_element(int i,double a)
{
  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  mats_[bi]->set_element(bo,a);
}

void
BlockedDiagSCMatrix::accumulate_element(int i,double a)
{
  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  mats_[bi]->accumulate_element(bo,a);
}

void
BlockedDiagSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const BlockedDiagSCMatrix* la = require_dynamic_cast<const BlockedDiagSCMatrix*>(a,
                               "BlockedDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedDiagSCMatrix:: accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i])
      mats_[i]->accumulate(la->mats_[i].pointer());
}

double
BlockedDiagSCMatrix::invert_this()
{
  double det = 1.0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i])
      det *= mats_[i]->invert_this();

  return det;
}

double
BlockedDiagSCMatrix::determ_this()
{
  double det = 1.0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i])
      det *= mats_[i]->determ_this();

  return det;
}

double
BlockedDiagSCMatrix::trace()
{
  double det = 0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i])
      det += mats_[i]->trace();

  return det;
}

void
BlockedDiagSCMatrix::gen_invert_this(double condition_number_threshold)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i])
      mats_[i]->gen_invert_this(condition_number_threshold);
}

void
BlockedDiagSCMatrix::convert_accumulate(DiagSCMatrix*a)
{
  // ok, are we converting a non-blocked matrix to a blocked one, or are
  // we converting a blocked matrix from one specialization to another?

  if (dynamic_cast<BlockedDiagSCMatrix*>(a)) {
    BlockedDiagSCMatrix *ba = dynamic_cast<BlockedDiagSCMatrix*>(a);
    if (ba->nblocks() == this->nblocks()) {
      for (int i=0; i < nblocks(); i++)
        mats_[i]->convert_accumulate(ba->mats_[i]);
    } else {
      ExEnv::errn() << indent
           << "BlockedDiagSCMatrix::convert_accumulate: "
           << "can't convert from BlockedDiagSCMatrix with different nblock"
           << endl;
      abort();
    }
  }
  else {
    if (nblocks()==1) {
      mats_[0]->convert_accumulate(a);
    } else {
      ExEnv::errn() << indent
           << "BlockedDiagSCMatrix::convert_accumulate: "
           << "can't convert from DiagSCMatrix when nblocks != 1"
           << endl;
      abort();
    }
  }
}

void
BlockedDiagSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  BlockedSCElementOp *bop = dynamic_cast<BlockedSCElementOp*>(op.pointer());

  int nb = d->blocks()->nblock();

  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedDiagSCMatrix::element_op(const Ref<SCElementOp2>& op,
                              DiagSCMatrix* m)
{
  BlockedDiagSCMatrix *lm = require_dynamic_cast<BlockedDiagSCMatrix*>(m,
                                    "BlockedDiagSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    ExEnv::errn() << indent << "BlockedDiagSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = dynamic_cast<BlockedSCElementOp2*>(op.pointer());

  int nb = d->blocks()->nblock();

  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedDiagSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  BlockedDiagSCMatrix *lm = require_dynamic_cast<BlockedDiagSCMatrix*>(m,
                                      "BlockedDiagSCMatrix::element_op");
  BlockedDiagSCMatrix *ln = require_dynamic_cast<BlockedDiagSCMatrix*>(n,
                                      "BlockedDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    ExEnv::errn() << indent << "BlockedDiagSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = dynamic_cast<BlockedSCElementOp3*>(op.pointer());

  int nb = d->blocks()->nblock();

  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op,lm->mats_[i].pointer(),ln->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedDiagSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->blocks()->nblock(); i++) {
    if (mats_[i] == 0)
      continue;

    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCDimension
BlockedDiagSCMatrix::dim(int i) const
{
  return d->blocks()->subdim(i);
}

int
BlockedDiagSCMatrix::nblocks() const
{
  return d->blocks()->nblock();
}

RefDiagSCMatrix
BlockedDiagSCMatrix::block(int i)
{
  return mats_[i];
}

Ref<SCMatrixSubblockIter>
BlockedDiagSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i) == 0)
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->local_blocks(access));
    }
  Ref<SCMatrixSubblockIter> ret = iter.pointer();
  return ret;
}

Ref<SCMatrixSubblockIter>
BlockedDiagSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i) == 0)
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->all_blocks(access));
    }
  Ref<SCMatrixSubblockIter> ret = iter.pointer();
  return ret;
}

void
BlockedDiagSCMatrix::save(StateOut&s)
{
  int ndim = n();
  s.put(ndim);
  int has_subblocks = 1;
  s.put(has_subblocks);
  s.put(nblocks());
  for (int i=0; i<nblocks(); i++) {
      block(i).save(s);
    }
}

void
BlockedDiagSCMatrix::restore(StateIn& s)
{
  int ndimt, ndim = n();
  s.get(ndimt);
  if (ndimt != ndim) {
      ExEnv::errn() << indent
           << "BlockedDiagSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          ExEnv::errn() << indent
               << "BlockedDiagSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      ExEnv::errn() << indent
           << "BlockedDiagSCMatrix::restore(): no subblocks--cannot restore"
           << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
