//
// distvect.cc
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
#include <iomanip.h>

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// DistSCVector member functions

#define CLASSNAME DistSCVector
#define PARENTS public SCVector
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DistSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

DistSCVector::DistSCVector(const RefSCDimension&a, DistSCMatrixKit*k):
  SCVector(a,k)
{
  init_blocklist();
}

int
DistSCVector::block_to_node(int i) const
{
  return (i)%messagegrp()->n();
}

RefSCMatrixBlock
DistSCVector::block_to_block(int i) const
{
  int offset = i;
  int nproc = messagegrp()->n();

  if ((offset%nproc) != messagegrp()->me()) return 0;

  SCMatrixBlockListIter I;
  for (I=blocklist->begin(); I!=blocklist->end(); I++) {
      if (I.block()->blocki == i)
          return I.block();
    }

  ExEnv::err() << indent << "DistSCVector::block_to_block: internal error" << endl;
  abort();
  return 0;
}

double *
DistSCVector::find_element(int i) const
{
  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  RefSCVectorSimpleBlock blk = block_to_block(bi);
  if (blk.nonnull()) {
      return &blk->dat()[oi];
    }
  else {
      return 0;
    }
}

int
DistSCVector::element_to_node(int i) const
{
  int bi, oi;
  d->blocks()->elem_to_block(i, bi, oi);

  return block_to_node(bi);
}

void
DistSCVector::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  SCMatrixBlock *b;
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc != me) continue;
      b = new SCVectorSimpleBlock(d->blocks()->start(i),
                                  d->blocks()->fence(i));
      b->blocki = i;
      blocklist->insert(b);
    }
}

DistSCVector::~DistSCVector()
{
}

double
DistSCVector::get_element(int i) const
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
DistSCVector::set_element(int i,double a)
{
  double *e = find_element(i);
  if (e) {
      *e = a;
    }
}

void
DistSCVector::accumulate_element(int i,double a)
{
  double *e = find_element(i);
  if (e) {
      *e += a;
    }
}

void
DistSCVector::accumulate(const SCVector*a)
{
  // make sure that the argument is of the correct type
  const DistSCVector* la
    = DistSCVector::require_const_castdown(a,"DistSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent << "DistSCVector::accumulate(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter i1, i2;
  for (i1=la->blocklist->begin(),i2=blocklist->begin();
       i1!=la->blocklist->end() && i2!=blocklist->end();
       i1++,i2++) {
      int n = i1.block()->ndat();
      if (n != i2.block()->ndat()) {
          ExEnv::err() << indent << "DistSCVector::accumulate "
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

void
DistSCVector::accumulate(const SCMatrix*a)
{
  // make sure that the argument is of the correct type
  const DistSCMatrix* la
    = DistSCMatrix::require_const_castdown(a,"DistSCVector::accumulate");

  // make sure that the dimensions match
  if (!((la->rowdim()->equiv(dim()) && la->coldim()->n() == 1)
        || (la->coldim()->equiv(dim()) && la->rowdim()->n() == 1))) {
      ExEnv::err() << indent << "DistSCVector::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter I, J;
  for (I = la->blocklist->begin(), J = blocklist->begin();
       I != la->blocklist->end() && J != blocklist->end();
       I++,J++) {
      int n = I.block()->ndat();
      if (n != J.block()->ndat()) {
          ExEnv::err() << indent << "DistSCVector::accumulate(SCMatrix*a): "
               << "block lists do not match" << endl;
          abort();
        }
      double *dati = I.block()->dat();
      double *datj = J.block()->dat();
      for (int i=0; i<n; i++) datj[i] += dati[i];
    }
}

void
DistSCVector::assign_v(SCVector*a)
{
  // make sure that the argument is of the correct type
  DistSCVector* la
    = DistSCVector::require_castdown(a,"DistSCVector::assign_v");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent << "DistSCVector::assign_v(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter i1, i2;
  for (i1=la->blocklist->begin(),i2=blocklist->begin();
       i1!=la->blocklist->end() && i2!=blocklist->end();
       i1++,i2++) {
      int n = i1.block()->ndat();
      if (n != i2.block()->ndat()) {
          ExEnv::err() << indent << "DistSCVector::assign "
               << "mismatch: internal error" << endl;
          abort();
        }
      double *dat1 = i1.block()->dat();
      double *dat2 = i2.block()->dat();
      for (int i=0; i<n; i++) {
          dat2[i] = dat1[i];
        }
    }
}

void
DistSCVector::assign_p(const double*a)
{
  SCMatrixBlockListIter I;
  for (I=blocklist->begin();
       I!=blocklist->end();
       I++) {
      RefSCVectorSimpleBlock b = SCVectorSimpleBlock::castdown(I.block());
      if (b.null()) {
          ExEnv::err() << indent << "DistSCVector::assign "
               << "mismatch: internal error" << endl;
          abort();
        }
      int n = b->ndat();
      const double *dat1 = &a[b->istart];
      double *dat2 = b->dat();
      for (int i=0; i<n; i++) {
          dat2[i] = dat1[i];
        }
    }
}

double
DistSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  DistSCVector* la
    = DistSCVector::require_castdown(a,"DistSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent << "DistSCVector::scalar_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  double result = 0.0;
  SCMatrixBlockListIter i1, i2;
  for (i1=la->blocklist->begin(),i2=blocklist->begin();
       i1!=la->blocklist->end() && i2!=blocklist->end();
       i1++,i2++) {
      int n = i1.block()->ndat();
      if (n != i2.block()->ndat()) {
          ExEnv::err() << indent << "DistSCVector::scalar_product: block mismatch: "
               << "internal error" << endl;
          abort();
        }
      double *dat1 = i1.block()->dat();
      double *dat2 = i2.block()->dat();
      for (int i=0; i<n; i++) {
          result += dat2[i] * dat1[i];
        }
    }
  messagegrp()->sum(result);
  return result;
}

void
DistSCVector::element_op(const RefSCElementOp& op)
{
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

void
DistSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  DistSCVector *lm
      = DistSCVector::require_castdown(m, "DistSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::err() << indent << "DistSCVector: bad element_op\n";
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
DistSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  DistSCVector *lm
      = DistSCVector::require_castdown(m, "DistSCVector::element_op");
  DistSCVector *ln
      = DistSCVector::require_castdown(n, "DistSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::err() << indent << "DistSCVector: bad element_op\n";
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

void
DistSCVector::accumulate_product_rv(SCMatrix *pa, SCVector *pb)
{
  const char* name = "DistSCMatrix::accumulate_product_rv";
  // make sure that the arguments are of the correct type
  DistSCMatrix* a = DistSCMatrix::require_castdown(pa,name);
  DistSCVector* b = DistSCVector::require_castdown(pb,name);

  // make sure that the dimensions match
  if (!dim()->equiv(a->rowdim()) || !a->coldim()->equiv(b->dim())) {
      ExEnv::err() << indent
           << "DistSCVector::accumulate_product_rv(SCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  a->create_vecform(DistSCMatrix::Row);
  a->vecform_op(DistSCMatrix::CopyToVec);

  int n = dim()->n();
  double *res = new double[n];

  for (int i=0; i<n; i++) res[i] = 0.0;

  RefSCMatrixSubblockIter I = b->all_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      RefSCVectorSimpleBlock blk
          = SCVectorSimpleBlock::castdown(I->block());
      int n = blk->iend - blk->istart;
      int joff = blk->istart;
      double *data = blk->data;
      for (int i=0; i<a->nvec; i++) {
          double *aveci = a->vec[i];
          for (int j=0; j<n; j++) {
              res[i+a->vecoff] += aveci[j+joff]*data[j];
            }
        }
    }

  a->delete_vecform();

  messagegrp()->sum(res, n);
  I = local_blocks(SCMatrixSubblockIter::Accum);
  for (I->begin(); I->ready(); I->next()) {
      RefSCVectorSimpleBlock blk
          = SCVectorSimpleBlock::castdown(I->block());
      int n = blk->iend - blk->istart;
      int ioff = blk->istart;
      double *data = blk->data;
      for (int i=0; i<n; i++) {
          data[i] += res[i+ioff];
        }
    }
  I = 0;

  delete[] res;
}

void
DistSCVector::convert(double *res) const
{
  int n = dim()->n();
  for (int i=0; i<n; i++) res[i] = 0.0;

  RefSCMatrixSubblockIter I = ((DistSCVector*)this)->local_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      RefSCVectorSimpleBlock blk
          = SCVectorSimpleBlock::castdown(I->block());
      int ni = blk->iend - blk->istart;
      int ioff = blk->istart;
      double *data = blk->data;
      for (int i=0; i<ni; i++) {
          res[i+ioff] = data[i];
        }
    }

  messagegrp()->sum(res, n);
}

void
DistSCVector::convert(SCVector *v)
{
  SCVector::convert(v);
}

RefSCMatrixSubblockIter
DistSCVector::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new SCMatrixListSubblockIter(access, blocklist);
}

RefSCMatrixSubblockIter
DistSCVector::all_blocks(SCMatrixSubblockIter::Access access)
{
  return new DistSCMatrixListSubblockIter(access, blocklist, messagegrp());
}

void
DistSCVector::vprint(const char *title, ostream& os, int prec) const
{
  double *data = new double[dim()->n()];
  convert(data);

  int i;
  int lwidth;
  double max=this->maxabs();

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec+5+(int) max;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if (messagegrp()->me() == 0) {
    if (title)
      os << endl << indent << title << endl;
    else
      os << endl;

    if (n()==0) {
      os << indent << "empty vector\n";
      return;
    }

    for (i=0; i<n(); i++)
      os << indent << setw(5) << i+1 << setw(lwidth) << data[i] << endl;
    os << endl;

    os.flush();
  }

  delete[] data;
}

void
DistSCVector::error(const char *msg)
{
  ExEnv::err() << indent << "DistSCVector: error: " << msg << endl;
}

RefDistSCMatrixKit
DistSCVector::skit()
{
  return DistSCMatrixKit::castdown(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
