//
// distsymm.cc
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
#include <math/scmat/disthql.h>

using namespace std;
using namespace sc;

extern "C" { int DBmalloc_chain_check(const char *, int, int); }

/////////////////////////////////////////////////////////////////////////////
// DistSymmSCMatrix member functions

static ClassDesc DistSymmSCMatrix_cd(
  typeid(DistSymmSCMatrix),"DistSymmSCMatrix",1,"public SymmSCMatrix",
  0, 0, 0);

DistSymmSCMatrix::DistSymmSCMatrix(const RefSCDimension&a,DistSCMatrixKit*k):
  SymmSCMatrix(a,k)
{
  init_blocklist();
}

int
DistSymmSCMatrix::block_to_node(int i, int j) const
{
  if (j>i) {
      ExEnv::errn() << indent << "DistSymmSCMatrix::block_to_node: j>i" << endl;
      abort();
    }

  return ((i*(i+1))/2 + j)%messagegrp()->n();
}

Ref<SCMatrixBlock>
DistSymmSCMatrix::block_to_block(int i, int j) const
{
  if (j>i) {
      ExEnv::errn() << indent << "DistSymmSCMatrix::block_to_block: j>i" << endl;
      abort();
    }

  int offset = (i*(i+1))/2 + j;
  int nproc = messagegrp()->n();

  if ((offset%nproc) != messagegrp()->me()) return 0;

  SCMatrixBlockListIter I;
  for (I=blocklist->begin(); I!=blocklist->end(); I++) {
      if (I.block()->blocki == i && I.block()->blockj == j)
          return I.block();
    }

  ExEnv::errn() << indent << "DistSymmSCMatrix::block_to_block: internal error" << endl;
  abort();
  return 0;
}

double *
DistSymmSCMatrix::find_element(int i, int j) const
{
  if (j>i) {
      int tmp = i; i=j; j=tmp;
    }

  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  int bj, oj;
  d->blocks()->elem_to_block(j, bj, oj);

  Ref<SCMatrixBlock> ablk = block_to_block(bi, bj);
  if (ablk.nonnull()) {
      if (bi != bj) {
          Ref<SCMatrixRectBlock> blk
              = dynamic_cast<SCMatrixRectBlock*>(ablk.pointer());
          if (blk.null()) return 0;
          return &blk->data[oi*(blk->jend-blk->jstart)+oj];
        }
      else {
          Ref<SCMatrixLTriBlock> blk
              = dynamic_cast<SCMatrixLTriBlock*>(ablk.pointer());
          if (blk.null()) return 0;
          return &blk->data[(oi*(oi+1))/2+oj];
        }
    }
  return 0;
}

int
DistSymmSCMatrix::element_to_node(int i, int j) const
{
  if (j>i) {
      int tmp = i; i=j; j=tmp;
    }

  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  int bj, oj;
  d->blocks()->elem_to_block(j, bj, oj);

  return block_to_node(bi,bj);
}

void
DistSymmSCMatrix::init_blocklist()
{
  int i, j, index;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  SCMatrixBlock *b;
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d->blocks()->nblock(); i++) {
      for (j=0; j<i; j++, index++) {
          if (index%nproc != me) continue;
          b = new SCMatrixRectBlock(d->blocks()->start(i),
                                    d->blocks()->fence(i),
                                    d->blocks()->start(j),
                                    d->blocks()->fence(j));
          b->blocki = i;
          b->blockj = j;
          blocklist->insert(b);
        }
      if (index%nproc == me) {
          b = new SCMatrixLTriBlock(d->blocks()->start(i),
                                    d->blocks()->fence(i));
          b->blocki = i;
          b->blockj = i;
          blocklist->insert(b);
        }
      index++;
    }
}

DistSymmSCMatrix::~DistSymmSCMatrix()
{
}

double
DistSymmSCMatrix::get_element(int i,int j) const
{
  double res;
  double *e = find_element(i,j);
  if (e) {
      res = *e;
      messagegrp()->bcast(res, messagegrp()->me());
    }
  else {
      messagegrp()->bcast(res, element_to_node(i, j));
    }
  return res;
}

void
DistSymmSCMatrix::set_element(int i,int j,double a)
{
  double *e = find_element(i,j);
  if (e) {
      *e = a;
    }
}

void
DistSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  double *e = find_element(i,j);
  if (e) {
      *e += a;
    }
}

SymmSCMatrix *
DistSymmSCMatrix::get_subblock(int br, int er)
{
  error("get_subblock");
  return 0;
}

SCMatrix *
DistSymmSCMatrix::get_subblock(int, int, int, int)
{
  error("get_subblock");
  return 0;
}

void
DistSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  error("assign_subblock");
}

void
DistSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  error("accumulate_subblock");
}

void
DistSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  error("accumulate_subblock");
}

void
DistSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  error("accumulate_subblock");
}

SCVector *
DistSymmSCMatrix::get_row(int i)
{
  error("get_row");
  return 0;
}

void
DistSymmSCMatrix::assign_row(SCVector *v, int i)
{
  error("assign_row");
}

void
DistSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  error("accumulate_row");
}

void
DistSymmSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const DistSymmSCMatrix* la
    = require_dynamic_cast<const DistSymmSCMatrix*>(a,"DistSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "DistSymmSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter i1, i2;
  for (i1=la->blocklist->begin(),i2=blocklist->begin();
       i1!=la->blocklist->end() && i2!=blocklist->end();
       i1++,i2++) {
      int n = i1.block()->ndat();
      if (n != i2.block()->ndat()) {
          ExEnv::errn() << indent
               << "DistSymmSCMatrixListSubblockIter::accumulate block "
               << "mismatch: internal error" << endl;
          abort();
        }
      double *dat1 = i1.block()->dat();
      double *dat2 = i2.block()->dat();
      for (int i=0; i<n; i++) {
          dat2[i] += dat1[i];
        }
    }
}

double
DistSymmSCMatrix::invert_this()
{
  RefDiagSCMatrix refa = kit()->diagmatrix(d);
  RefSCMatrix refb = kit()->matrix(d,d);
  diagonalize(refa.pointer(),refb.pointer());
  double determ = 1.0;
  for (int i=0; i<dim()->n(); i++) {
      double val = refa->get_element(i);
      determ *= val;
    }
  Ref<SCElementOp> op = new SCElementInvert(1.0e-12);
  refa->element_op(op.pointer());
  assign(0.0);
  accumulate_transform(refb.pointer(), refa.pointer());
  return determ;
}

double
DistSymmSCMatrix::determ_this()
{
  return invert_this();
}

double
DistSymmSCMatrix::trace()
{
  double ret=0.0;
  Ref<SCMatrixSubblockIter> I = local_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixLTriBlock> b = dynamic_cast<SCMatrixLTriBlock*>(I->block());
      if (b.nonnull() && b->blocki == b->blockj) {
          int ni = b->end-b->start;
          double *data = b->data;
          for (int i=0; i<ni; i++) {
              data += i;
              ret += *data;
            }
        }
    }
  messagegrp()->sum(ret);

  return ret;
}

double
DistSymmSCMatrix::solve_this(SCVector*v)
{
  DistSCVector* lv =
    require_dynamic_cast<DistSCVector*>(v,"DistSymmSCMatrix::solve_this");

  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
      ExEnv::errn() << indent << "DistSymmSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match\n";
      abort();
    }

  error("no solve this");

  return 0.0;
}

void
DistSymmSCMatrix::gen_invert_this(double condition_number_threshold)
{
  invert_this();
}

void
DistSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  const char* name = "DistSymmSCMatrix::diagonalize";
  // make sure that the argument are of the correct type
  require_dynamic_cast<DistDiagSCMatrix*>(a,name);
  DistSCMatrix* lb = require_dynamic_cast<DistSCMatrix*>(b,name);

  int n = dim()->n();
  int me = messagegrp()->me();
  int nproc = messagegrp()->n();

  RefSCMatrix arect = kit()->matrix(dim(),dim());
  DistSCMatrix *rect = dynamic_cast<DistSCMatrix*>(arect.pointer());
  rect->assign(0.0);
  rect->accumulate(this);

  // This sets up the index list of columns to be stored on this node
  int nvec = n/nproc + (me<(n%nproc)?1:0);
  int *ivec = new int[nvec];
  for (int i=0; i<nvec; i++) {
      ivec[i] = i*nproc + me;
    }

  rect->create_vecform(DistSCMatrix::Col,nvec);
  rect->vecform_op(DistSCMatrix::CopyToVec,ivec);
  lb->create_vecform(DistSCMatrix::Col,nvec);

  double *d = new double[n];
  dist_diagonalize(n, rect->nvec, rect->vec[0], d, lb->vec[0],
                   messagegrp());

  // put d into the diagonal matrix
  a->assign(d);

  lb->vecform_op(DistSCMatrix::CopyFromVec, ivec);
  lb->delete_vecform();
  rect->delete_vecform();
  arect = 0;
  delete[] ivec;
}

void
DistSymmSCMatrix::eigensystem(SymmSCMatrix*s, DiagSCMatrix*a, SCMatrix*b)
{
  throw FeatureNotImplemented("DistSymmSCMatrix::eigensystem", __FILE__, __LINE__, this->class_desc());
}

// computes this += a + a.t
void
DistSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  DistSCMatrix* la
    = require_dynamic_cast<DistSCMatrix*>(a,"DistSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "DistSymmSCMatrix::"
           << "accumulate_symmetric_sum(SCMatrix*a): bad dim\n";
      abort();
    }

  Ref<SCMatrixSubblockIter> I = all_blocks(SCMatrixSubblockIter::Accum);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixBlock> block = I->block();
      // see if i've got this block
      Ref<SCMatrixBlock> localblock
          = la->block_to_block(block->blocki,block->blockj);
      if (localblock.nonnull()) {
          // the diagonal blocks require special handling
          if (block->blocki == block->blockj) {
              int n = la->rowblocks()->size(block->blocki);
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<n; i++) {
                  for (int j=0; j<=i; j++) {
                      double tmp = 0.0;
                      tmp += dat2[i*n+j];
                      tmp += dat2[j*n+i];
                      *dat1 += tmp;
                      dat1++;
                    }
                }
            }
          else {
              int n = block->ndat();
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<n; i++) {
                  dat1[i] += dat2[i];
                }
            }
        }
      // now for the transpose
      if (block->blocki != block->blockj) {
          localblock = la->block_to_block(block->blockj,block->blocki);
          if (localblock.nonnull()) {
              int nr = la->rowblocks()->size(block->blocki);
              int nc = la->rowblocks()->size(block->blockj);
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<nr; i++) {
                  for (int j=0; j<nc; j++) {
                      *dat1++ += dat2[j*nr+i];
                    }
                }
            }
        }
    }
}

void
DistSymmSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

void
DistSymmSCMatrix::element_op(const Ref<SCElementOp2>& op,
                              SymmSCMatrix* m)
{
  DistSymmSCMatrix *lm
      = require_dynamic_cast<DistSymmSCMatrix*>(m,"DistSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::errn() << indent << "DistSymmSCMatrix: bad element_op\n";
      abort();
    }
  SCMatrixBlockListIter i, j;
  for (i = blocklist->begin(), j = lm->blocklist->begin();
       i != blocklist->end();
       i++, j++) {
      op->process_base(i.block(), j.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

void
DistSymmSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  DistSymmSCMatrix *lm
      = require_dynamic_cast<DistSymmSCMatrix*>(m,"DistSymSCMatrix::element_op");
  DistSymmSCMatrix *ln
      = require_dynamic_cast<DistSymmSCMatrix*>(n,"DistSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::errn() << indent << "DistSymmSCMatrix: bad element_op\n";
      abort();
    }
  SCMatrixBlockListIter i, j, k;
  for (i = blocklist->begin(),
           j = lm->blocklist->begin(),
           k = ln->blocklist->begin();
       i != blocklist->end();
       i++, j++, k++) {
      op->process_base(i.block(), j.block(), k.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

Ref<SCMatrixSubblockIter>
DistSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new SCMatrixListSubblockIter(access, blocklist);
}

Ref<SCMatrixSubblockIter>
DistSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  return new DistSCMatrixListSubblockIter(access, blocklist, messagegrp());
}

void
DistSymmSCMatrix::error(const char *msg)
{
  ExEnv::errn() << "DistSymmSCMatrix: error: " << msg << endl;
}

Ref<DistSCMatrixKit>
DistSymmSCMatrix::skit()
{
  return dynamic_cast<DistSCMatrixKit*>(kit().pointer());
}

void
DistSymmSCMatrix::convert_accumulate(SymmSCMatrix*a)
{
  SymmSCMatrix::convert_accumulate(a);

#if 0
  DistSymmSCMatrix *d = require_dynamic_cast<DistSymmSCMatrix*>(a,
                                 "DistSymmSCMatrix::convert_accumulate");

  SCMatrixBlockListIter i, j;
  for (i = blocklist->begin(), j = d->blocklist->begin();
       i != blocklist->end();
       i++, j++) {
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
