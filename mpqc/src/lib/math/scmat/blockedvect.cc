//
// blockedvect.cc
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
// BlockedSCVector member functions

#define CLASSNAME BlockedSCVector
#define PARENTS public SCVector
#include <util/class/classi.h>
void *
BlockedSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSCVector::resize(SCDimension *bsd)
{
  if (vecs_) {
    delete[] vecs_;
    vecs_=0;
  }

  d = bsd;
  
  if (!bsd || !bsd->blocks()->nblock())
    return;
  
  vecs_ = new RefSCVector[d->blocks()->nblock()];

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (d->blocks()->size(i))
      vecs_[i] = subkit->vector(d->blocks()->subdim(i));
}

BlockedSCVector::BlockedSCVector(const RefSCDimension&a,
                                 BlockedSCMatrixKit*k):
  SCVector(a,k),
  subkit(k->subkit()),
  vecs_(0)
{
  resize(a);
}

BlockedSCVector::~BlockedSCVector()
{
  if (vecs_) {
    delete[] vecs_;
    vecs_=0;
  }
}

void
BlockedSCVector::assign(double a)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(a);
}

void
BlockedSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedSCVector::assign(SCVector*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(la->vecs_[i]);
}

void
BlockedSCVector::assign(const double*a)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->assign(a+d->blocks()->start(i));
}

double
BlockedSCVector::get_element(int i)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    cerr << indent << "BlockedSCVector::get_element: bad index\n";
    abort();
  }

  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  return vecs_[bi].get_element(bo);
}

void
BlockedSCVector::set_element(int i,double a)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    cerr << indent << "BlockedSCVector::set_element: bad index\n";
    abort();
  }

  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  vecs_[bi].set_element(bo,a);
}

void
BlockedSCVector::accumulate_element(int i,double a)
{
  int size = d->n();
  if (i < 0 || i >= size) {
    cerr << indent << "BlockedSCVector::accumulate_element: bad index\n";
    abort();
  }

  int bi, bo;
  d->blocks()->elem_to_block(i,bi,bo);
  vecs_[bi].accumulate_element(bo,a);
}

void
BlockedSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "BlockedSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = BlockedSCMatrix::require_castdown(a,name);
  BlockedSCVector* lb = BlockedSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
    cerr << indent
         << "BlockedSCVector::accumulate_product(SCMatrix*a,SCVector*b): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate_product(la->mats_[i], lb->vecs_[i]);
}

void
BlockedSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "BlockedSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  BlockedSymmSCMatrix* la = BlockedSymmSCMatrix::require_castdown(a,name);
  BlockedSCVector* lb = BlockedSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
    cerr << indent
         << "BlockedSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate_product(la->mats_[i], lb->vecs_[i]);
}

void
BlockedSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedSCVector::accumulate(SCVector*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate(la->vecs_[i]);
}

void
BlockedSCVector::accumulate(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSCVector::accumulate");

  // make sure that the dimensions match
  if (!((la->rowdim()->equiv(dim()) && la->coldim()->n() == 1)
        || (la->coldim()->equiv(dim()) && la->rowdim()->n() == 1))) {
    cerr << indent << "BlockedSCVector::accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      vecs_[i]->accumulate(la->mats_[i]);
}

double
BlockedSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedSCVector::scale_product(SCVector*a): "
         << "dimensions don't match\n";
    abort();
  }

  double result=0;

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (vecs_[i].nonnull())
      result += vecs_[i]->scalar_product(la->vecs_[i]);
  
  return result;
}

void
BlockedSCVector::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  BlockedSCVector *lm
      = BlockedSCVector::require_castdown(m, "BlockedSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
    cerr << indent << "BlockedSCVector: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = BlockedSCElementOp2::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op, lm->vecs_[i]);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  BlockedSCVector *lm
      = BlockedSCVector::require_castdown(m, "BlockedSCVector::element_op");
  BlockedSCVector *ln
      = BlockedSCVector::require_castdown(n, "BlockedSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    cerr << indent << "BlockedSCVector: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = BlockedSCElementOp3::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (vecs_[i].nonnull())
      vecs_[i]->element_op(op, lm->vecs_[i], ln->vecs_[i]);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCVector::vprint(const char *title, ostream& os, int prec)
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < d->blocks()->nblock(); i++) {
    if (vecs_[i].null())
      continue;
    
    sprintf(newtitle,"%s:  block %d",title,i+1);
    vecs_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCDimension
BlockedSCVector::dim(int i)
{
  return d->blocks()->subdim(i);
}

int
BlockedSCVector::nblocks() const
{
  return d->blocks()->nblock();
}

RefSCVector
BlockedSCVector::block(int i)
{
  return vecs_[i];
}

RefSCMatrixSubblockIter
BlockedSCVector::local_blocks(SCMatrixSubblockIter::Access access)
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
BlockedSCVector::all_blocks(SCMatrixSubblockIter::Access access)
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
BlockedSCVector::save(StateOut&s)
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
BlockedSCVector::restore(StateIn&s)
{
  int ndimt, ndim = n();
  s.get(ndimt);
  if (ndimt != ndim) {
      cerr << indent << "BlockedSCVector::restore(): bad dimension" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          cerr << indent
               << "BlockedSCVector::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      cerr << indent
           << "BlockedSCVector::restore(): no subblocks--cannot restore"
           << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
