//
// distdiag.cc
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

#include <iostream.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

#define DEBUG 0

/////////////////////////////////////////////////////////////////////////////
// DistDiagSCMatrix member functions

#define CLASSNAME DistDiagSCMatrix
#define PARENTS public DiagSCMatrix
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DistDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DiagSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

DistDiagSCMatrix::DistDiagSCMatrix(const RefSCDimension&a,DistSCMatrixKit*k):
  DiagSCMatrix(a,k)
{
  init_blocklist();
}

int
DistDiagSCMatrix::block_to_node(int i) const
{
  return (i)%messagegrp()->n();
}

RefSCMatrixBlock
DistDiagSCMatrix::block_to_block(int i) const
{
  int offset = i;
  int nproc = messagegrp()->n();

  if ((offset%nproc) != messagegrp()->me()) return 0;

  SCMatrixBlockListIter I;
  for (I=blocklist->begin(); I!=blocklist->end(); I++) {
      if (I.block()->blocki == i)
          return I.block();
    }

  cerr << indent
       << "DistDiagSCMatrix::block_to_block: internal error" << endl;
  abort();
  return 0;
}

double *
DistDiagSCMatrix::find_element(int i) const
{
  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  if (DEBUG)
      ExEnv::out() << messagegrp()->me() << ": "
                   << "find_element(" << i << "): "
                   << "block = " << bi << ", "
                   << "offset = " << oi
                   << endl;

  RefSCMatrixDiagBlock blk = block_to_block(bi);
  if (blk.nonnull()) {
      if (DEBUG)
          ExEnv::out() << messagegrp()->me() << ": ndat = " << blk->ndat() << endl;
      if (oi >= blk->ndat()) {
          cerr << messagegrp()->me() << ": DistDiagSCMatrix::find_element"
               << ": internal error" << endl;
          abort();
        }
      return &blk->dat()[oi];
    }
  else {
      if (DEBUG)
          ExEnv::out() << messagegrp()->me() << ": can't find" << endl;
      return 0;
    }
}

int
DistDiagSCMatrix::element_to_node(int i) const
{
  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  return block_to_node(bi);
}

void
DistDiagSCMatrix::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  SCMatrixBlock *b;
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc != me) continue;
      b = new SCMatrixDiagBlock(d->blocks()->start(i),
                                d->blocks()->fence(i),
                                d->blocks()->start(i));
      b->blocki = i;
      b->blockj = i;
      blocklist->insert(b);
    }
}

DistDiagSCMatrix::~DistDiagSCMatrix()
{
}

double
DistDiagSCMatrix::get_element(int i) const
{
  double res;
  double *e = find_element(i);
  if (e) {
      res = *e;
      messagegrp()->bcast(res, messagegrp()->me());
    }
  else {
      messagegrp()->bcast(res, element_to_node(i));
    }
  return res;
}

void
DistDiagSCMatrix::set_element(int i,double a)
{
  double *e = find_element(i);
  if (e) {
      *e = a;
    }
}

void
DistDiagSCMatrix::accumulate_element(int i,double a)
{
  double *e = find_element(i);
  if (e) {
      *e += a;
    }
}

void
DistDiagSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const DistDiagSCMatrix* la
    = DistDiagSCMatrix::require_const_castdown(a,"DistDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "DistDiagSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter i1, i2;
  for (i1=la->blocklist->begin(),i2=blocklist->begin();
       i1!=la->blocklist->end() && i2!=blocklist->end();
       i1++,i2++) {
      int n = i1.block()->ndat();
      if (n != i2.block()->ndat()) {
          cerr << indent << "DistDiagSCMatrix::accumulate "
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
DistDiagSCMatrix::invert_this()
{
  RefSCMatrixSubblockIter I = local_blocks(SCMatrixSubblockIter::Read);
  double det = 1.0;
  for (I->begin(); I->ready(); I->next()) {
      int n = I->block()->ndat();
      double *data = I->block()->dat();
      for (int i=0; i<n; i++) {
          det *= data[i];
          data[i] = 1.0/data[i];
        }
    }
  GrpProductReduce<double> gred;
  messagegrp()->reduce(&det, 1, gred);
  return det;
}

double
DistDiagSCMatrix::determ_this()
{
  RefSCMatrixSubblockIter I = local_blocks(SCMatrixSubblockIter::Read);
  double det = 1.0;
  for (I->begin(); I->ready(); I->next()) {
      int n = I->block()->ndat();
      double *data = I->block()->dat();
      for (int i=0; i<n; i++) {
          det *= data[i];
        }
    }
  GrpProductReduce<double> gred;
  messagegrp()->reduce(det, gred);
  return det;
}

double
DistDiagSCMatrix::trace()
{
  double ret=0.0;
  RefSCMatrixSubblockIter I = local_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      int n = I->block()->ndat();
      double *data = I->block()->dat();
      for (int i=0; i<n; i++) {
          ret += data[i];
        }
    }
  messagegrp()->sum(ret);
  return ret;
}

void
DistDiagSCMatrix::gen_invert_this()
{
  RefSCMatrixSubblockIter I = local_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      int n = I->block()->ndat();
      double *data = I->block()->dat();
      for (int i=0; i<n; i++) {
          if (fabs(data[i]) > 1.0e-8)
              data[i] = 1.0/data[i];
          else
              data[i] = 0.0;
        }
    }
}

void
DistDiagSCMatrix::element_op(const RefSCElementOp& op)
{
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

void
DistDiagSCMatrix::element_op(const RefSCElementOp2& op,
                              DiagSCMatrix* m)
{
  DistDiagSCMatrix *lm
      = DistDiagSCMatrix::require_castdown(m,"DistDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      cerr << indent << "DistDiagSCMatrix: bad element_op\n";
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
DistDiagSCMatrix::element_op(const RefSCElementOp3& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  DistDiagSCMatrix *lm
      = DistDiagSCMatrix::require_castdown(m,"DistDiagSCMatrix::element_op");
  DistDiagSCMatrix *ln
      = DistDiagSCMatrix::require_castdown(n,"DistDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      cerr << indent << "DistDiagSCMatrix: bad element_op\n";
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

RefSCMatrixSubblockIter
DistDiagSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new SCMatrixListSubblockIter(access, blocklist);
}

RefSCMatrixSubblockIter
DistDiagSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  return new DistSCMatrixListSubblockIter(access, blocklist, messagegrp());
}

void
DistDiagSCMatrix::error(const char *msg)
{
  cerr << indent << "DistDiagSCMatrix: error: " << msg << endl;
}

RefDistSCMatrixKit
DistDiagSCMatrix::skit()
{
  return DistSCMatrixKit::castdown(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
