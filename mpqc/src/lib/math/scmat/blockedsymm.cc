//
// blockedsymm.cc
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
// BlockedSymmSCMatrix member functions

#define CLASSNAME BlockedSymmSCMatrix
#define PARENTS public SymmSCMatrix
#include <util/class/classi.h>
void *
BlockedSymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SymmSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
BlockedSymmSCMatrix::resize(SCDimension *a)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  
  d = a;
  
  mats_ = new RefSymmSCMatrix[d->blocks()->nblock()];
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (d->blocks()->size(i))
      mats_[i] = subkit->symmmatrix(d->blocks()->subdim(i));
}

BlockedSymmSCMatrix::BlockedSymmSCMatrix(const RefSCDimension&a,
                                         BlockedSCMatrixKit*k):
  SymmSCMatrix(a,k),
  subkit(k->subkit()),
  mats_(0)
{
  resize(a);
}

BlockedSymmSCMatrix::~BlockedSymmSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
}

double
BlockedSymmSCMatrix::get_element(int i,int j)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d->blocks()->elem_to_block(i,block_i,elem_i);
  d->blocks()->elem_to_block(j,block_j,elem_j);

  if (block_i != block_j)
    return 0;
  
  return mats_[block_i]->get_element(elem_i,elem_j);
}

void
BlockedSymmSCMatrix::set_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d->blocks()->elem_to_block(i,block_i,elem_i);
  d->blocks()->elem_to_block(j,block_j,elem_j);

  if (block_i != block_j)
    return;
  
  mats_[block_i]->set_element(elem_i,elem_j,a);
}

void
BlockedSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d->blocks()->elem_to_block(i,block_i,elem_i);
  d->blocks()->elem_to_block(j,block_j,elem_j);

  if (block_i != block_j)
    return;
  
  mats_[block_i]->accumulate_element(elem_i,elem_j,a);
}

SCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  cerr << indent << "BlockedSymmSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

SymmSCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er)
{
  cerr << indent << "BlockedSymmSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

void
BlockedSymmSCMatrix::assign_subblock(SCMatrix*sb,
                                     int br, int er, int bc, int ec)
{
  cerr << indent << "BlockedSymmSCMatrix::assign_subblock:"
       << " cannot assign subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  cerr << indent << "BlockedSymmSCMatrix::assign_subblock:"
                 << " cannot assign subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SCMatrix*sb,
                                         int br, int er, int bc, int ec)
{
  cerr << indent << "BlockedSymmSCMatrix::accumulate_subblock:"
                 << " cannot accumulate subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  cerr << indent << "BlockedSymmSCMatrix::accumulate_subblock:"
                 << " cannot accumulate subblock\n";
  abort();
}

SCVector *
BlockedSymmSCMatrix::get_row(int i)
{
  cerr << indent << "BlockedSymmSCMatrix::get_row: cannot get row\n";
  abort();

  return 0;
}

void
BlockedSymmSCMatrix::assign_row(SCVector *v, int i)
{
  cerr << indent << "BlockedSymmSCMatrix::assign_row: cannot assign row\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  cerr << indent << "BlockedSymmSCMatrix::accumulate_row:"
                 << " cannot accumulate row\n";
  abort();
}

double
BlockedSymmSCMatrix::invert_this()
{
  double res=1;
  
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      res *= mats_[i]->invert_this();

  return res;
}

double
BlockedSymmSCMatrix::determ_this()
{
  double res=1;
  
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSymmSCMatrix::trace()
{
  double res=0;
  
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      res += mats_[i]->determ_this();

  return res;
}

double
BlockedSymmSCMatrix::solve_this(SCVector*v)
{
  double res=1;
  
  BlockedSCVector* lv =
    BlockedSCVector::require_castdown(v,"BlockedSymmSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
    cerr << indent << "BlockedSymmSCMatrix::solve_this(SCVector*v): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      res *= mats_[i]->solve_this(lv->vecs_[i].pointer());

  return res;
}

void
BlockedSymmSCMatrix::scale(double s)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->scale(s);
}

void
BlockedSymmSCMatrix::gen_invert_this()
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->gen_invert_this();
}

double
BlockedSymmSCMatrix::scalar_product(SCVector*a)
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

  double result = 0.0;
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      result += mats_[i]->scalar_product(la->vecs_[i].pointer());

  return result;
}

void
BlockedSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "BlockedSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  BlockedDiagSCMatrix* la = BlockedDiagSCMatrix::require_castdown(a,name);
  BlockedSCMatrix* lb = BlockedSCMatrix::require_castdown(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
    cerr << indent << "BlockedSymmSCMatrix::"
         << "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->diagonalize(la->mats_[i].pointer(),lb->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::accumulate(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSymmSCMatrix* la = BlockedSymmSCMatrix::require_castdown(a,
                                   "BlockedSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedSymmSCMatrix::accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate(la->mats_[i].pointer());
}

// computes this += a * a.t
void
BlockedSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  int i, zero=0;
  
  // make sure that the argument is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
    cerr << indent << "BlockedSymmSCMatrix::"
         << "accumulate_symmetric_product(SCMatrix*a): bad dim\n";
    abort();
  }

  int mxnb = (d->blocks()->nblock() > la->nblocks_)
             ? d->blocks()->nblock()
             : la->nblocks_;
  int &mi = (d->blocks()->nblock()==1) ? zero : i;

  for (i=0; i < mxnb; i++)
    if (mats_[mi].nonnull() && la->mats_[i].nonnull())
      mats_[mi]->accumulate_symmetric_product(la->mats_[i].pointer());
}

// computes this += a + a.t
void
BlockedSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  BlockedSCMatrix* la
    = BlockedSCMatrix::require_castdown(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
    cerr << indent << "BlockedSymmSCMatrix::"
         << "accumulate_symmetric_sum(SCMatrix*a): bad dim\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate_symmetric_sum(la->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = BlockedSCVector::require_castdown(a,"BlockedSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
    cerr << indent << "BlockedSymmSCMatrix::"
         << "accumulate_symmetric_outer_product(SCMatrix*a): bad dim\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->accumulate_symmetric_outer_product(la->vecs_[i].pointer());
}

// this += a * b * transpose(a)
void
BlockedSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b,
                                          SCMatrix::Transform t)
{
  int i, zero=0;
  
  // do the necessary castdowns
  BlockedSCMatrix*la
    = BlockedSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  BlockedSymmSCMatrix*lb = require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (t == SCMatrix::NormalTransform) {
      if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
          cerr << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim\n";
          abort();
        }
    }
  else {
      if (!dim()->equiv(la->coldim()) || !lb->dim()->equiv(la->rowdim())) {
          cerr << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim\n";
          abort();
        }
    }

  int mxnb = (d->blocks()->nblock() > la->nblocks_)
             ? d->blocks()->nblock()
             : la->nblocks_;

  int &mi = (d->blocks()->nblock()==1) ? zero : i;
  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;
  
  for (i=0; i < mxnb; i++) {
    if (mats_[mi].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
                             
    mats_[mi]->accumulate_transform(la->mats_[i].pointer(),
                                    lb->mats_[bi].pointer(),t);
  }
}

// this += a * b * transpose(a)
void
BlockedSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b,
                                          SCMatrix::Transform t)
{
  int i, zero=0;
  
  // do the necessary castdowns
  BlockedSCMatrix*la
    = BlockedSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  BlockedDiagSCMatrix*lb
    = BlockedDiagSCMatrix::require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
    cerr << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim\n";
    abort();
  }

  int mxnb = (d->blocks()->nblock() > la->nblocks_)
             ? d->blocks()->nblock()
             : la->nblocks_;

  int &mi = (d->blocks()->nblock()==1) ? zero : i;
  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;

  for (i=0; i < mxnb; i++) {
    if (mats_[mi].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
                             
    mats_[mi]->accumulate_transform(la->mats_[i].pointer(),
                                    lb->mats_[bi].pointer());
  }
}

void
BlockedSymmSCMatrix::accumulate_transform(SymmSCMatrix*a,SymmSCMatrix*b)
{
  SymmSCMatrix::accumulate_transform(a,b);
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp& op)
{
  BlockedSCElementOp *bop = BlockedSCElementOp::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp2& op,
                                SymmSCMatrix* m)
{
  BlockedSymmSCMatrix *lm = BlockedSymmSCMatrix::require_castdown(m,
                                          "BlockedSymSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    cerr << indent << "BlockedSymmSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = BlockedSCElementOp2::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSymmSCMatrix::element_op(const RefSCElementOp3& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  BlockedSymmSCMatrix *lm = BlockedSymmSCMatrix::require_castdown(m,
                                        "BlockedSymSCMatrix::element_op");
  BlockedSymmSCMatrix *ln = BlockedSymmSCMatrix::require_castdown(n,
                                        "BlockedSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    cerr << indent << "BlockedSymmSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = BlockedSCElementOp3::castdown(op.pointer());

  int nb = d->blocks()->nblock();
  
  op->defer_collect(1);
  for (int i=0; i < nb; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i].nonnull())
      mats_[i]->element_op(op,lm->mats_[i].pointer(),
                            ln->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSymmSCMatrix::vprint(const char *title, ostream& os, int prec)
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
BlockedSymmSCMatrix::dim(int i)
{
  return d->blocks()->subdim(i);
}

int
BlockedSymmSCMatrix::nblocks() const
{
  return d->blocks()->nblock();
}

RefSymmSCMatrix
BlockedSymmSCMatrix::block(int i)
{
  return mats_[i];
}

RefSCMatrixSubblockIter
BlockedSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
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
BlockedSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
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
BlockedSymmSCMatrix::save(StateOut&s)
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
BlockedSymmSCMatrix::restore(StateIn& s)
{
  int ndimt, ndim = n();
  s.get(ndimt);
  if (ndimt != ndim) {
      cerr << indent
           << "BlockedSymmSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          cerr << indent
               << "BlockedSymmSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      cerr << indent
           << "BlockedSymmSCMatrix::restore(): no subblocks--cannot restore"
           << endl;
      abort();
    }
}

void
BlockedSymmSCMatrix::convert_accumulate(SymmSCMatrix*a)
{
  // ok, are we converting a non-blocked matrix to a blocked one, or are
  // we converting a blocked matrix from one specialization to another?
  
  if (BlockedSymmSCMatrix::castdown(a)) {
    BlockedSymmSCMatrix *ba = BlockedSymmSCMatrix::castdown(a);
    if (ba->nblocks() == this->nblocks()) {
      for (int i=0; i < nblocks(); i++)
        mats_[i]->convert_accumulate(ba->mats_[i]);
    } else {
      cerr << indent
           << "BlockedSymmSCMatrix::convert_accumulate: "
           << "I can't do that"
           << endl;
      abort();
    }
  }
  else {
    if (nblocks()==1) {
      mats_[0]->convert_accumulate(a);
    } else {
      cerr << indent
           << "BlockedSymmSCMatrix::convert_accumulate: "
           << "I can't do that either"
           << endl;
      abort();
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
