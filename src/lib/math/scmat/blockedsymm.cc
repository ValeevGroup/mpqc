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
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// BlockedSymmSCMatrix member functions

static ClassDesc BlockedSymmSCMatrix_cd(
  typeid(BlockedSymmSCMatrix),"BlockedSymmSCMatrix",1,"public SymmSCMatrix",
  0, 0, 0);

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
BlockedSymmSCMatrix::get_element(int i,int j) const
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

void
BlockedSymmSCMatrix::assign_p(const double*a)
{
  // if this is a single, complete block, then use the block's
  // assign member.  Otherwise, use the generic member.
  if (d->blocks()->nblock() == 1
      && (dim()->n() == mats_[0]->dim()->n())) {
      mats_[0]->assign_p(a);
    }
  else {
      SymmSCMatrix::assign_p(a);
    }
}

void
BlockedSymmSCMatrix::assign_pp(const double**a)
{
  // if this is a single, complete block, then use the block's
  // assign member.  Otherwise, use the generic member.
  if (d->blocks()->nblock() == 1
      && (dim()->n() == mats_[0]->dim()->n())) {
      mats_[0]->assign_pp(a);
    }
  else {
      SymmSCMatrix::assign_pp(a);
    }
}

void
BlockedSymmSCMatrix::convert_p(double*a) const
{
  // if this is a single, complete block, then use the block's
  // convert member.  Otherwise, use the generic member.
  if (d->blocks()->nblock() == 1
      && (dim()->n() == mats_[0]->dim()->n())) {
      mats_[0]->convert_p(a);
    }
  else {
      SymmSCMatrix::convert_p(a);
    }
}

void
BlockedSymmSCMatrix::convert_pp(double**a) const
{
  // if this is a single, complete block, then use the block's
  // convert member.  Otherwise, use the generic member.
  if (d->blocks()->nblock() == 1
      && (dim()->n() == mats_[0]->dim()->n())) {
      mats_[0]->convert_pp(a);
    }
  else {
      SymmSCMatrix::convert_pp(a);
    }
}

SCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

SymmSCMatrix *
BlockedSymmSCMatrix::get_subblock(int br, int er)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

void
BlockedSymmSCMatrix::assign_subblock(SCMatrix*sb,
                                     int br, int er, int bc, int ec)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::assign_subblock:"
       << " cannot assign subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::assign_subblock:"
                 << " cannot assign subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SCMatrix*sb,
                                         int br, int er, int bc, int ec)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate_subblock:"
                 << " cannot accumulate subblock\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate_subblock:"
                 << " cannot accumulate subblock\n";
  abort();
}

SCVector *
BlockedSymmSCMatrix::get_row(int i)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::get_row: cannot get row\n";
  abort();

  return 0;
}

void
BlockedSymmSCMatrix::assign_row(SCVector *v, int i)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::assign_row: cannot assign row\n";
  abort();
}

void
BlockedSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate_row:"
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
      res += mats_[i]->trace();

  return res;
}

double
BlockedSymmSCMatrix::solve_this(SCVector*v)
{
  double res=1;

  BlockedSCVector* lv =
    require_dynamic_cast<BlockedSCVector*>(v,"BlockedSymmSCMatrix::solve_this");

  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::solve_this(SCVector*v): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      res *= mats_[i]->solve_this(lv->vecs_[i].pointer());

  return res;
}

void
BlockedSymmSCMatrix::assign_val(double s)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->assign(s);
}

void
BlockedSymmSCMatrix::assign_s(SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  BlockedSymmSCMatrix* la = require_dynamic_cast<BlockedSymmSCMatrix*>(a,
                                   "BlockedSymmSCMatrix::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::assign_s(SymmSCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->assign(la->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::scale(double s)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->scale(s);
}

void
BlockedSymmSCMatrix::gen_invert_this(double condition_number_threshold)
{
  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->gen_invert_this(condition_number_threshold);
}

double
BlockedSymmSCMatrix::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  BlockedSCVector* la
    = require_dynamic_cast<BlockedSCVector*>(a,"BlockedSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSCVector::scale_product(SCVector*a): "
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
  BlockedDiagSCMatrix* la = require_dynamic_cast<BlockedDiagSCMatrix*>(a,name);
  BlockedSCMatrix* lb = require_dynamic_cast<BlockedSCMatrix*>(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::"
         << "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++)
    if (mats_[i].nonnull())
      mats_[i]->diagonalize(la->mats_[i].pointer(),lb->mats_[i].pointer());
}

void
BlockedSymmSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const BlockedSymmSCMatrix* la = require_dynamic_cast<const BlockedSymmSCMatrix*>(a,
                                   "BlockedSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate(SymmSCMatrix*a): "
         << "dimensions don't match\n";
    ExEnv::errn() << indent << "this->dim():" << std::endl;
    dim()->print(ExEnv::errn());
    ExEnv::errn() << indent << "la->dim():" << std::endl;
    la->dim()->print(ExEnv::errn());
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
    = require_dynamic_cast<BlockedSCMatrix*>(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::"
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
    = require_dynamic_cast<BlockedSCMatrix*>(a,"BlockedSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::"
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
    = require_dynamic_cast<BlockedSCVector*>(a,"BlockedSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::"
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
    = require_dynamic_cast<BlockedSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  BlockedSymmSCMatrix*lb = require_dynamic_cast<BlockedSymmSCMatrix*>(
      b,"%s::accumulate_transform", class_name());

  // check the dimensions
  if (t == SCMatrix::NormalTransform) {
      if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
          ExEnv::outn() << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim (not transposed)\n";
          ExEnv::outn() << "target dim:" << endl;
          dim()->print();
          ExEnv::outn() << "source dim" << endl;
          b->dim()->print();
          ExEnv::outn() << "transform dims" << endl;
          a->rowdim()->print();
          a->coldim()->print();
          abort();
        }
    }
  else {
      if (!dim()->equiv(la->coldim()) || !lb->dim()->equiv(la->rowdim())) {
          ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim\n";
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
    = require_dynamic_cast<BlockedSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  BlockedDiagSCMatrix*lb
    = require_dynamic_cast<BlockedDiagSCMatrix*>(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::accumulate_transform: bad dim\n";
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
BlockedSymmSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  BlockedSCElementOp *bop = dynamic_cast<BlockedSCElementOp*>(op.pointer());

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
BlockedSymmSCMatrix::element_op(const Ref<SCElementOp2>& op,
                                SymmSCMatrix* m)
{
  BlockedSymmSCMatrix *lm = require_dynamic_cast<BlockedSymmSCMatrix*>(m,
                                          "BlockedSymSCMatrix::element_op");
  if (!dim()->equiv(lm->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = dynamic_cast<BlockedSCElementOp2*>(op.pointer());

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
BlockedSymmSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  BlockedSymmSCMatrix *lm = require_dynamic_cast<BlockedSymmSCMatrix*>(m,
                                        "BlockedSymSCMatrix::element_op");
  BlockedSymmSCMatrix *ln = require_dynamic_cast<BlockedSymmSCMatrix*>(n,
                                        "BlockedSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = dynamic_cast<BlockedSCElementOp3*>(op.pointer());

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
BlockedSymmSCMatrix::vprint(const char *title, ostream& os, int prec) const
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
BlockedSymmSCMatrix::dim(int i) const
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

Ref<SCMatrixSubblockIter>
BlockedSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->local_blocks(access));
    }
  Ref<SCMatrixSubblockIter> ret = iter.pointer();
  return ret;
}

Ref<SCMatrixSubblockIter>
BlockedSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access,nblocks());
  for (int i=0; i<nblocks(); i++) {
      if (block(i).null())
          iter->set_iter(i, new SCMatrixNullSubblockIter(access));
      else
          iter->set_iter(i, block(i)->all_blocks(access));
    }
  Ref<SCMatrixSubblockIter> ret = iter.pointer();
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
      ExEnv::errn() << indent
           << "BlockedSymmSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          ExEnv::errn() << indent
               << "BlockedSymmSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      ExEnv::errn() << indent
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

  if (dynamic_cast<BlockedSymmSCMatrix*>(a)) {
    BlockedSymmSCMatrix *ba = dynamic_cast<BlockedSymmSCMatrix*>(a);
    if (ba->nblocks() == this->nblocks()) {
      for (int i=0; i < nblocks(); i++)
        mats_[i]->convert_accumulate(ba->mats_[i]);
    } else {
      ExEnv::errn() << indent
           << "BlockedSymmSCMatrix::convert_accumulate: "
           << "can't convert from BlockedSymmSCMatrix with different nblock"
           << endl;
      abort();
    }
  }
  else {
    if (nblocks()==1) {
      mats_[0]->convert_accumulate(a);
    } else {
      ExEnv::errn() << indent
           << "BlockedSymmSCMatrix::convert_accumulate: "
           << "can't convert from SymmSCMatrix when nblocks != 1"
           << endl;
      abort();
    }
  }
}

void
BlockedSymmSCMatrix::eigensystem(SymmSCMatrix* s, DiagSCMatrix*a, SCMatrix*b) {

  const char* name = "BlockedSymmSCMatrix::eigensystem";
  // make sure that the arguments is of the correct type
  BlockedSymmSCMatrix* ls = require_dynamic_cast<BlockedSymmSCMatrix*>(s,name);
  BlockedDiagSCMatrix* la = require_dynamic_cast<BlockedDiagSCMatrix*>(a,name);
  BlockedSCMatrix* lb = require_dynamic_cast<BlockedSCMatrix*>(b,name);

  if (!dim()->equiv(ls->dim()) ||
      !dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) ||
      !dim()->equiv(lb->rowdim())) {
    ExEnv::errn() << indent << "BlockedSymmSCMatrix::"
                  << "eigensystem(SymmSCMatrix*s,DiagSCMatrix*a,SCMatrix*b): bad dims\n";
    abort();
  }

  for (int i=0; i < d->blocks()->nblock(); i++) {
    if (mats_[i].nonnull())
      mats_[i]->eigensystem(ls->mats_[i].pointer(),
                            la->mats_[i].pointer(),
                            lb->mats_[i].pointer());
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
