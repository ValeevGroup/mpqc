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
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedDiagSCMatrix member functions

#define CLASSNAME BlockedDiagSCMatrix
#define PARENTS public DiagSCMatrix
#include <util/class/classi.h>
void *
BlockedDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DiagSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

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
BlockedDiagSCMatrix::get_element(int i)
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
BlockedDiagSCMatrix::accumulate(DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedDiagSCMatrix* la = BlockedDiagSCMatrix::require_castdown(a,
                               "BlockedDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedDiagSCMatrix:: accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->mats_[i].pointer());
}

double
BlockedDiagSCMatrix::invert_this()
{
  double det = 1.0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      det *= mats_[i]->invert_this();
  
  return det;
}

double
BlockedDiagSCMatrix::determ_this()
{
  double det = 1.0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      det *= mats_[i]->determ_this();
  
  return det;
}

double
BlockedDiagSCMatrix::trace()
{
  double det = 0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      det += mats_[i]->trace();
  
  return det;
}

void
BlockedDiagSCMatrix::gen_invert_this()
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->gen_invert_this();
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  for (int i=0; i < nb; i++) {
    if (i < nb-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);

    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op);
  }
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp2& op,
                              DiagSCMatrix* m)
{
  BlockedDiagSCMatrix *lm = BlockedDiagSCMatrix::require_castdown(m,
                                    "BlockedDiagSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    cerr << indent << "BlockedDiagSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = BlockedSCElementOp2::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  for (int i=0; i < nb; i++) {
    if (i < nb-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);

    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
}

void
BlockedDiagSCMatrix::element_op(const RefSCElementOp3& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  BlockedDiagSCMatrix *lm = BlockedDiagSCMatrix::require_castdown(m,
                                      "BlockedDiagSCMatrix::element_op");
  BlockedDiagSCMatrix *ln = BlockedDiagSCMatrix::require_castdown(n,
                                      "BlockedDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    cerr << indent << "BlockedDiagSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = BlockedSCElementOp3::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  for (int i=0; i < nb; i++) {
    if (i < nb-1)
      op->defer_collect(1);
    else
      op->defer_collect(0);

    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op,lm->mats_[i].pointer(),ln->mats_[i].pointer());
  }
}

void
BlockedDiagSCMatrix::vprint(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->blocks()->nblock(); i++) {
    if (mats_[i].null())
      continue;
    
    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCDimension
BlockedDiagSCMatrix::dim(int i)
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

RefSCMatrixSubblockIter
BlockedDiagSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixCompositeSubblockIter iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->local_blocks(access));
    }
  RefSCMatrixSubblockIter ret = iter.pointer();
  return ret;
}

RefSCMatrixSubblockIter
BlockedDiagSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  RefSCMatrixCompositeSubblockIter iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->all_blocks(access));
    }
  RefSCMatrixSubblockIter ret = iter.pointer();
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
      cerr << indent
           << "BlockedDiagSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          cerr << indent
               << "BlockedDiagSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      cerr << indent
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
