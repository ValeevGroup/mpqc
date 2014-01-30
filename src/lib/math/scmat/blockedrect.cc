//
// blockedrect.cc
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

#include <cmath>
#include <cassert>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

#include <math/scmat/local.h>
#include <math/scmat/repl.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// BlockedSCMatrix member functions

static ClassDesc BlockedSCMatrix_cd(
  typeid(BlockedSCMatrix),"BlockedSCMatrix",1,"public SCMatrix",
  0, 0, 0);

void
BlockedSCMatrix::resize(SCDimension *a, SCDimension *b)
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }

  d1 = a;
  d2 = b;

  if (!a || !b || !a->blocks()->nblock() || !b->blocks()->nblock())
    return;

  if (a->blocks()->nblock() > 1 && b->blocks()->nblock() == 1) {
    nblocks_ = d1->blocks()->nblock();
    mats_ = new RefSCMatrix[d1->blocks()->nblock()];
    for (int i=0; i < d1->blocks()->nblock(); i++)
      if (d1->blocks()->size(i) && d2->blocks()->size(0))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(i),
                                  d2->blocks()->subdim(0));

  } else if (a->blocks()->nblock() == 1 && b->blocks()->nblock() > 1) {
    nblocks_ = d2->blocks()->nblock();
    mats_ = new RefSCMatrix[d2->blocks()->nblock()];
    for (int i=0; i < d2->blocks()->nblock(); i++)
      if (d2->blocks()->size(i) && d1->blocks()->size(0))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(0),
                                  d2->blocks()->subdim(i));

  } else if (a->blocks()->nblock() == b->blocks()->nblock()) {
    nblocks_ = d2->blocks()->nblock();
    mats_ = new RefSCMatrix[d1->blocks()->nblock()];
    for (int i=0; i < d1->blocks()->nblock(); i++)
      if (d2->blocks()->size(i) && d1->blocks()->size(i))
        mats_[i] = subkit->matrix(d1->blocks()->subdim(i),
                                  d2->blocks()->subdim(i));

  } else {
    ExEnv::errn() << indent << "BlockedSCMatrix::resize: wrong number of blocks\n";
    abort();
  }

}

BlockedSCMatrix::BlockedSCMatrix(const RefSCDimension&a,
                                 const RefSCDimension&b,
                                 BlockedSCMatrixKit*k):
  SCMatrix(a,b,k),
  subkit(k->subkit()),
  mats_(0), nblocks_(0)
{
  resize(a,b);
}

BlockedSCMatrix::~BlockedSCMatrix()
{
  if (mats_) {
    delete[] mats_;
    mats_=0;
  }
  nblocks_=0;
}

void
BlockedSCMatrix::assign_val(double v)
{
  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->assign(v);
}

double
BlockedSCMatrix::get_element(int i,int j) const
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);

  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    return mats_[block_j]->get_element(elem_i,elem_j);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1) {
    return mats_[block_i]->get_element(elem_i,elem_j);
  } else if (d1->blocks()->nblock() == d2->blocks()->nblock()
             && block_i == block_j) {
    return mats_[block_i]->get_element(elem_i,elem_j);
  } else {
    return 0;
  }
}

void
BlockedSCMatrix::set_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);

  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    mats_[block_j]->set_element(elem_i,elem_j,a);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1) {
    mats_[block_i]->set_element(elem_i,elem_j,a);
  } else if (d1->blocks()->nblock() == d2->blocks()->nblock()
             && block_i == block_j) {
    mats_[block_i]->set_element(elem_i,elem_j,a);
  }
}

void
BlockedSCMatrix::accumulate_element(int i,int j,double a)
{
  int block_i, block_j;
  int elem_i, elem_j;

  d1->blocks()->elem_to_block(i,block_i,elem_i);
  d2->blocks()->elem_to_block(j,block_j,elem_j);

  if (d1->blocks()->nblock() == 1 && d2->blocks()->nblock() > 1) {
    mats_[block_j]->accumulate_element(elem_i,elem_j,a);

  } else if (d1->blocks()->nblock() > 1 && d2->blocks()->nblock() == 1) {
    mats_[block_i]->accumulate_element(elem_i,elem_j,a);
  }
  else if (d1->blocks()->nblock() == d2->blocks()->nblock()
           && block_i == block_j) {
    mats_[block_i]->accumulate_element(elem_i,elem_j,a);
  }
}

void
BlockedSCMatrix::assign_p(const double*a)
{
  // if this is a single, complete block, then use the block's
  // assign member.  Otherwise, use the generic member.
  if (nblocks_ == 1
      && (rowdim()->n() == mats_[0]->rowdim()->n())
      && (coldim()->n() == mats_[0]->coldim()->n())) {
      mats_[0]->assign_p(a);
    }
  else {
      SCMatrix::assign_p(a);
    }
}

void
BlockedSCMatrix::assign_pp(const double**a)
{
  // if this is a single, complete block, then use the block's
  // assign member.  Otherwise, use the generic member.
  if (nblocks_ == 1
      && (rowdim()->n() == mats_[0]->rowdim()->n())
      && (coldim()->n() == mats_[0]->coldim()->n())) {
      mats_[0]->assign_pp(a);
    }
  else {
      SCMatrix::assign_pp(a);
    }
}

void
BlockedSCMatrix::convert_p(double*a) const
{
  // if this is a single, complete block, then use the block's
  // convert member.  Otherwise, use the generic member.
  if (nblocks_ == 1
      && (rowdim()->n() == mats_[0]->rowdim()->n())
      && (coldim()->n() == mats_[0]->coldim()->n())) {
      mats_[0]->convert_p(a);
    }
  else {
      SCMatrix::convert_p(a);
    }
}

void
BlockedSCMatrix::convert_pp(double**a) const
{
  // if this is a single, complete block, then use the block's
  // convert member.  Otherwise, use the generic member.
  if (nblocks_ == 1
      && (rowdim()->n() == mats_[0]->rowdim()->n())
      && (coldim()->n() == mats_[0]->coldim()->n())) {
      mats_[0]->convert_pp(a);
    }
  else {
      SCMatrix::convert_pp(a);
    }
}

void
BlockedSCMatrix::convert_accumulate(SCMatrix*a)
{
  // ok, are we converting a non-blocked matrix to a blocked one, or are
  // we converting a blocked matrix from one specialization to another?

  MPQC_ASSERT(a->nrow() == this->nrow());
  MPQC_ASSERT(a->ncol() == this->ncol());

  if (dynamic_cast<BlockedSCMatrix*>(a)) {
    BlockedSCMatrix *ba = dynamic_cast<BlockedSCMatrix*>(a);
    if (ba->nblocks() == this->nblocks()) {
      for (int i=0; i < nblocks(); i++)
        mats_[i]->convert_accumulate(ba->mats_[i]);
    } else {
      ExEnv::errn() << indent
           << "BlockedSCMatrix::convert_accumulate: "
           << "can't convert from BlockedSCMatrix with different nblock"
           << endl;
      abort();
    }
  }
  else {
    if (nblocks()==1) {
      mats_[0]->convert_accumulate(a);
    } else {
      int row_offset = 0, col_offset = 0;
      for (int i=0; i < nblocks(); i++) {
        mats_[i]->accumulate_subblock(a,
                                      0, mats_[i]->nrow()-1, 0, mats_[i]->ncol()-1,
                                      row_offset, col_offset);
        row_offset += mats_[i]->nrow();
        col_offset += mats_[i]->ncol();
      }
    }
  }
}

SCMatrix *
BlockedSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::get_subblock: cannot get subblock\n";
  abort();
  return 0;
}

void
BlockedSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                               int source_br, int source_bc)
{
  BlockedSCMatrix* b_sb
      = require_dynamic_cast<BlockedSCMatrix*>(
          sb,"BlockedSCMatrix::assign_subblock");
  if (! (nblocks() == 1 && b_sb->nblocks() == 1)) {
      ExEnv::errn() << indent << "BlockedSCMatrix::assign_subblock:"
                    << " cannot assign subblock\n";
      abort();
    }
  else {
      block(0).assign_subblock(b_sb->block(0),br,er,bc,ec,
                               source_br,source_bc);
    }
}

void
BlockedSCMatrix::accumulate_subblock(SCMatrix*sb,
                                     int br, int er, int bc, int ec,
                                     int source_br, int source_bc)
{
  BlockedSCMatrix* b_sb
      = require_dynamic_cast<BlockedSCMatrix*>(
          sb,"BlockedSCMatrix::accumulate_subblock");
  if (! (nblocks() == 1 && b_sb->nblocks() == 1)) {
      ExEnv::errn() << indent << "BlockedSCMatrix::accumulate_subblock:"
                    << " cannot accumulate subblock\n";
      abort();
    }
  else {
      block(0).accumulate_subblock(b_sb->block(0),br,er,bc,ec,
                                   source_br,source_bc);
    }
}

SCVector *
BlockedSCMatrix::get_row(int i)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::get_row: cannot get row\n";
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_row(SCVector *v, int i)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::assign_row: cannot assign row\n";
  abort();
}

void
BlockedSCMatrix::accumulate_row(SCVector *v, int i)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::accumulate_row: cannot accumulate row\n";
  abort();
}

SCVector *
BlockedSCMatrix::get_column(int i)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::get_column: cannot get column\n";
  abort();

  return 0;
}

void
BlockedSCMatrix::assign_column(SCVector *v, int i)
{
  ExEnv::errn() << indent << "BlockedSCMatrix::assign_column: cannot assign column\n";
  abort();
}

void
BlockedSCMatrix::accumulate_column(SCVector *v, int i)
{
  ExEnv::errn() << indent
       << "BlockedSCMatrix::accumulate_column: cannot accumulate column\n";
  abort();
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
BlockedSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "BlockedSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  BlockedSCVector* la = require_dynamic_cast<BlockedSCVector*>(a,name);
  BlockedSCVector* lb = require_dynamic_cast<BlockedSCVector*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::accumulate_outer_product(SCVector*,SCVector*): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < d1->blocks()->nblock(); i++)
    if (mats_[i])
      mats_[i]->accumulate_outer_product(la->vecs_[i], lb->vecs_[i]);
}

void
BlockedSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b)
{
  int i, zero = 0;

  const char* name = "BlockedSCMatrix::accumulate_product_rr";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = require_dynamic_cast<BlockedSCMatrix*>(a,name);
  BlockedSCMatrix* lb = require_dynamic_cast<BlockedSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b): "
         << "dimensions don't match\n";
    ExEnv::errn() << "this->rowdim: " << std::endl; rowdim().print();
    ExEnv::errn() << "this->coldim: " << std::endl; coldim().print();
    ExEnv::errn() << "   a->rowdim: " << std::endl; la->rowdim().print();
    ExEnv::errn() << "   a->coldim: " << std::endl; la->coldim().print();
    ExEnv::errn() << "   b->rowdim: " << std::endl; lb->rowdim().print();
    ExEnv::errn() << "   b->coldim: " << std::endl; lb->coldim().print();
    abort();
  }

  // find out the number of blocks we need to process.
  int mxnb = (nblocks_ > la->nblocks_) ? nblocks_ : la->nblocks_;

  int nrba = la->d1->blocks()->nblock();
  int ncba = la->d2->blocks()->nblock();
  int nrbb = lb->d1->blocks()->nblock();
  int ncbb = lb->d2->blocks()->nblock();

  int &mi = (nrba==1 && ncba > 1 && nrbb > 1 && ncbb==1) ? zero : i;
  int &ai = (nrba==1 && ncba==1) ? zero : i;
  int &bi = (nrbb==1 && ncbb==1) ? zero : i;

  for (i=0; i < mxnb; i++) {
    if (mats_[mi].null() || la->mats_[ai].null() || lb->mats_[bi].null())
      continue;
    mats_[mi]->accumulate_product(la->mats_[ai], lb->mats_[bi]);
  }
}

void
BlockedSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b)
{
  int i, zero=0;

  const char* name = "BlockedSCMatrix::accumulate_product_rs";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = require_dynamic_cast<BlockedSCMatrix*>(a,name);
  BlockedSymmSCMatrix* lb = require_dynamic_cast<BlockedSymmSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b): "
         << "dimensions don't match\n";
    rowdim()->print();
    coldim()->print();
    la->rowdim()->print();
    la->coldim()->print();
    lb->dim()->print();
    abort();
  }

  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;

  for (i=0; i < nblocks_; i++) {
    if (mats_[i].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[bi]);
  }
}


void
BlockedSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b)
{
  int i, zero=0;

  const char* name = "BlockedSCMatrix::accumulate_product_rd";
  // make sure that the arguments are of the correct type
  BlockedSCMatrix* la = require_dynamic_cast<BlockedSCMatrix*>(a,name);
  BlockedDiagSCMatrix* lb = require_dynamic_cast<BlockedDiagSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b): "
         << "dimensions don't match\n";
    abort();
  }

  int &bi = (lb->d->blocks()->nblock()==1) ? zero : i;

  for (i=0; i < nblocks_; i++) {
    if (mats_[i].null() || la->mats_[i].null() || lb->mats_[bi].null())
      continue;
    mats_[i]->accumulate_product(la->mats_[i], lb->mats_[bi]);
  }
}

void
BlockedSCMatrix::accumulate(const SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const BlockedSCMatrix* la
    = require_dynamic_cast<const BlockedSCMatrix*>(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix::accumulate(SCMatrix*a): "
         << "dimensions don't match\n";
    rowdim().print();
    coldim().print();
    la->rowdim().print();
    la->coldim().print();
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const BlockedSymmSCMatrix* la
    = require_dynamic_cast<const BlockedSymmSCMatrix*>(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix::accumulate(SymmSCMatrix*a): "
         << "dimensions don't match\n";
    rowdim().print();
    coldim().print();
    la->dim().print();
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const BlockedDiagSCMatrix* la
    = require_dynamic_cast<const BlockedDiagSCMatrix*>(a,"BlockedSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix::accumulate(DiagSCMatrix*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->accumulate(la->mats_[i]);
}

void
BlockedSCMatrix::accumulate(const SCVector*a)
{
  // make sure that the arguments is of the correct type
  const BlockedSCVector* la
    = require_dynamic_cast<const BlockedSCVector*>(a,"BlockedSCVector::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
    ExEnv::errn() << indent << "BlockedSCMatrix::accumulate(SCVector*a): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->accumulate(la->vecs_[i]);
}

void
BlockedSCMatrix::transpose_this()
{
  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->transpose_this();

  RefSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
}

// hack, hack, hack!  One day we'll get svd working everywhere.
double
BlockedSCMatrix::invert_this()
{
  int i;
  double res=1;

  // if this matrix is block diagonal, then give a normal inversion a shot
  if (d1->blocks()->nblock() == d2->blocks()->nblock()) {
    for (i=0; i < nblocks_; i++)
      if (mats_[i]) res *= mats_[i]->invert_this();
    return res;
  }

  // ok, let's make sure that the matrix is at least square
  if (d1->n() != d2->n()) {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::invert_this: SVD not implemented yet\n";
    abort();
  }

  if (d1->blocks()->nblock() == 1) {
    RefSCMatrix tdim = subkit->matrix(d1->blocks()->subdim(0),
                                      d1->blocks()->subdim(0));
    tdim->convert(this);
    res = tdim->invert_this();
    transpose_this();

    // d1 and d2 were swapped by now
    for (i=0; i < d1->blocks()->nblock(); i++)
      if (mats_[i])
        mats_[i]->convert(tdim.get_subblock(d1->blocks()->start(i),
                                            d1->blocks()->fence(i)-1,
                                            0, d2->n()-1));

    return res;

  } else if (d2->blocks()->nblock() == 1) {
    RefSCMatrix tdim = subkit->matrix(d2->blocks()->subdim(0),
                                      d2->blocks()->subdim(0));

    tdim->convert(this);
    res = tdim->invert_this();
    transpose_this();

    // d1 and d2 were swapped by now
    for (i=0; i < d2->blocks()->nblock(); i++)
      if (mats_[i])
        mats_[i]->convert(tdim.get_subblock(0, d1->n()-1,
                                            d2->blocks()->start(i),
                                            d2->blocks()->fence(i)-1));

    return res;

  } else {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::invert_this: SVD not implemented yet\n";
    abort();
  }

  return 0.0;
}

void
BlockedSCMatrix::gen_invert_this(double condition_number_threshold)
{
  // if this matrix is block-diagonal use SVD
  if (d1->blocks()->nblock() == d2->blocks()->nblock()) {
    for (int i=0; i < nblocks_; i++)
      if (mats_[i])
        mats_[i]->gen_invert_this(condition_number_threshold);
  }
  else {
    ExEnv::errn() << indent
         << "BlockedSCMatrix::gen_invert_this: SVD not implemented yet for non-block-diagonal matrices\n";
    abort();
  }
}

double
BlockedSCMatrix::determ_this()
{
  double res=1;

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      res *= mats_[i]->determ_this();

  return res;
}

double
BlockedSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      ret += mats_[i]->trace();

  return ret;
}

void
BlockedSCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  BlockedSCMatrix* lU =
    require_dynamic_cast<BlockedSCMatrix*>(U,"BlockedSCMatrix::svd_this");
  BlockedSCMatrix* lV =
    require_dynamic_cast<BlockedSCMatrix*>(V,"BlockedSCMatrix::svd_this");
  BlockedDiagSCMatrix* lsigma =
    require_dynamic_cast<BlockedDiagSCMatrix*>(sigma,"BlockedSCMatrix::svd_this");

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->svd_this(lU->mats_[i], lsigma->mats_[i], lV->mats_[i]);
}

double
BlockedSCMatrix::solve_this(SCVector*v)
{
  double res=1;

  BlockedSCVector* lv =
    require_dynamic_cast<BlockedSCVector*>(v,"BlockedSCMatrix::solve_this");

  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix::solve_this(SCVector*v): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      res *= mats_[i]->solve_this(lv->vecs_[i]);

  return res;
}

void
BlockedSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  BlockedSymmSCMatrix* lS =
    require_dynamic_cast<BlockedSymmSCMatrix*>(S,"BlockedSCMatrix::schmidt_orthog");

  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix::schmidt_orthog(): "
         << "dimensions don't match\n";
    abort();
  }

  for (int i=0; i < nblocks_; i++)
    if (mats_[i])
      mats_[i]->schmidt_orthog(lS->mats_[i].pointer(),
                               lS->dim()->blocks()->subdim(i).n());
}

int
BlockedSCMatrix::schmidt_orthog_tol(SymmSCMatrix *S, double tol, double *res)
{
  ExEnv::err0()
               << "ERROR: BlockedSCMatrix::schmidt_orthog_tol doesn't exist"
               << endl;
  abort();
  return 0;
}

void
BlockedSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  BlockedSCElementOp *bop = dynamic_cast<BlockedSCElementOp*>(op.pointer());

  op->defer_collect(1);
  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op);
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCMatrix::element_op(const Ref<SCElementOp2>& op,
                          SCMatrix* m)
{
  BlockedSCMatrix *lm
    = require_dynamic_cast<BlockedSCMatrix*>(m,"BlockedSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp2 *bop = dynamic_cast<BlockedSCElementOp2*>(op.pointer());

  op->defer_collect(1);
  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op,lm->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCMatrix::element_op(const Ref<SCElementOp3>& op,
                          SCMatrix* m,SCMatrix* n)
{
  BlockedSCMatrix *lm
    = require_dynamic_cast<BlockedSCMatrix*>(m,"BlockedSCMatrix::element_op");
  BlockedSCMatrix *ln
    = require_dynamic_cast<BlockedSCMatrix*>(n,"BlockedSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
    ExEnv::errn() << indent << "BlockedSCMatrix: bad element_op\n";
    abort();
  }

  BlockedSCElementOp3 *bop = dynamic_cast<BlockedSCElementOp3*>(op.pointer());

  op->defer_collect(1);
  for (int i=0; i < nblocks_; i++) {
    if (bop)
      bop->working_on(i);
    if (mats_[i])
      mats_[i]->element_op(op,lm->mats_[i].pointer(),
                              ln->mats_[i].pointer());
  }
  op->defer_collect(0);
  if (op->has_collect()) op->collect(messagegrp());
}

void
BlockedSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int len = (title) ? strlen(title) : 0;
  char *newtitle = new char[len + 80];

  for (int i=0; i < nblocks_; i++) {
    if (mats_[i].null())
      continue;

    sprintf(newtitle,"%s:  block %d",title,i+1);
    mats_[i]->print(newtitle, os, prec);
  }

  delete[] newtitle;
}

RefSCDimension
BlockedSCMatrix::rowdim(int i) const
{
  return d1->blocks()->subdim(i);
}

RefSCDimension
BlockedSCMatrix::coldim(int i) const
{
  return d2->blocks()->subdim(i);
}

int
BlockedSCMatrix::nblocks() const
{
  return nblocks_;
}

RefSCMatrix
BlockedSCMatrix::block(int i)
{
  if (mats_)
      return mats_[i];
  else
      return (SCMatrix*)0;
}

Ref<SCMatrixSubblockIter>
BlockedSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access, nblocks());
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
BlockedSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  Ref<SCMatrixCompositeSubblockIter> iter
      = new SCMatrixCompositeSubblockIter(access, nblocks());
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
BlockedSCMatrix::save(StateOut&s)
{
  int nr = nrow();
  int nc = ncol();
  s.put(nr);
  s.put(nc);
  int has_subblocks = 1;
  s.put(has_subblocks);
  s.put(nblocks());
  for (int i=0; i<nblocks(); i++) {
      block(i).save(s);
    }
}

void
BlockedSCMatrix::restore(StateIn& s)
{
  int nrt, nr = nrow();
  int nct, nc = ncol();
  s.get(nrt);
  s.get(nct);
  if (nrt != nr || nct != nc) {
      ExEnv::errn() << indent << "BlockedSCMatrix::restore(): bad dimensions" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (has_subblocks) {
      int nblock;
      s.get(nblock);
      if (nblock != nblocks()) {
          ExEnv::errn() << indent
               << "BlockedSCMatrix::restore(): nblock differs\n" << endl;
          abort();
        }
      for (int i=0; i<nblocks(); i++) {
          block(i).restore(s);
        }
    }
  else {
      ExEnv::errn() << indent
           << "BlockedSCMatrix::restore(): no subblocks--cannot restore"
           << endl;
      abort();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
