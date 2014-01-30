//
// distrect.cc
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
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/misc/consumableresources.h>
#include <math/scmat/dist.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

#define DEBUG 0

/////////////////////////////////////////////////////////////////////////////

static void
fail(const char *m)
{
  ExEnv::errn() << indent << "distrect.cc: error: " << m << endl;
  abort();
}

/////////////////////////////////////////////////////////////////////////////
// DistSCMatrix member functions

static ClassDesc DistSCMatrix_cd(
  typeid(DistSCMatrix),"DistSCMatrix",1,"public SCMatrix",
  0, 0, 0);

DistSCMatrix::DistSCMatrix(const RefSCDimension&a,const RefSCDimension&b,
                           DistSCMatrixKit*k):
  SCMatrix(a,b,k)
{
  init_blocklist();
}

int
DistSCMatrix::block_to_node(int i, int j) const
{
  return (i*d2->blocks()->nblock() + j)%messagegrp()->n();
}

Ref<SCMatrixBlock>
DistSCMatrix::block_to_block(int i, int j) const
{
  int offset = i*d2->blocks()->nblock() + j;
  int nproc = messagegrp()->n();

  if ((offset%nproc) != messagegrp()->me()) return 0;

  SCMatrixBlockListIter I;
  for (I=blocklist->begin(); I!=blocklist->end(); I++) {
      if (I.block()->blocki == i && I.block()->blockj == j)
          return I.block();
    }

  ExEnv::errn() << indent << "DistSCMatrix::block_to_block: internal error" << endl;
  abort();
  return 0;
}

double *
DistSCMatrix::find_element(int i, int j) const
{
  int bi, oi;
  d1->blocks()->elem_to_block(i, bi, oi);

  int bj, oj;
  d2->blocks()->elem_to_block(j, bj, oj);

  Ref<SCMatrixRectBlock> blk
      = dynamic_cast<SCMatrixRectBlock*>(block_to_block(bi, bj).pointer());
  if (blk) {
      return &blk->data[oi*(blk->jend-blk->jstart)+oj];
    }
  else {
      return 0;
    }
}

int
DistSCMatrix::element_to_node(int i, int j) const
{
  int bi, oi;
  d1->blocks()->elem_to_block(i, bi, oi);

  int bj, oj;
  d2->blocks()->elem_to_block(j, bj, oj);

  return block_to_node(bi,bj);
}

void
DistSCMatrix::init_blocklist()
{
  int i, j, index;
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d1->blocks()->nblock(); i++) {
      for (j=0; j<d2->blocks()->nblock(); j++, index++) {
          if (block_to_node(i,j) == me) {
              Ref<SCMatrixBlock> b
                  = new SCMatrixRectBlock(d1->blocks()->start(i),
                                          d1->blocks()->fence(i),
                                          d2->blocks()->start(j),
                                          d2->blocks()->fence(j));
              b->blocki = i;
              b->blockj = j;
              blocklist->append(b);
            }
        }
    }
}

DistSCMatrix::~DistSCMatrix()
{
}

double
DistSCMatrix::get_element(int i,int j) const
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
DistSCMatrix::set_element(int i,int j,double a)
{
  double *e = find_element(i,j);
  if (e) {
      *e = a;
    }
}

void
DistSCMatrix::accumulate_element(int i,int j,double a)
{
  double *e = find_element(i,j);
  if (e) {
      *e += a;
    }
}

SCMatrix *
DistSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  error("get_subblock");
  return 0;
}

void
DistSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                              int source_br, int source_bc)
{
  error("assign_subblock");
}

void
DistSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                  int source_br, int source_bc)
{
  error("accumulate_subblock");
}

SCVector *
DistSCMatrix::get_row(int i)
{
  error("get_row");
  return 0;
}

void
DistSCMatrix::assign_row(SCVector *v, int i)
{
  error("assign_row");
}

void
DistSCMatrix::accumulate_row(SCVector *v, int i)
{
  error("accumulate_row");
}

SCVector *
DistSCMatrix::get_column(int i)
{
  double *col = new double[nrow()];

  Ref<SCMatrixSubblockIter> iter = local_blocks(SCMatrixSubblockIter::Read);
  for (iter->begin(); iter->ready(); iter->next()) {
    SCMatrixRectBlock *b = dynamic_cast<SCMatrixRectBlock*>(iter->block());
    if (b->jstart > i || b->jend <= i)
      continue;

    int joff = i-b->jstart;
    int jlen = b->jend-b->jstart;
    int ist = 0;

    for (int ii=b->istart; ii < b->iend; ii++,ist++)
      col[ii] = b->data[ist*jlen+joff];
  }

  SCVector * rcol = kit_->vector(rowdim());
  rcol->assign(col);

  delete[] col;

  return rcol;
}

void
DistSCMatrix::assign_column(SCVector *v, int i)
{
  error("assign_column");
}

void
DistSCMatrix::accumulate_column(SCVector *v, int i)
{
  error("accumulate_column");
}

void
DistSCMatrix::accumulate(const SCMatrix*a)
{
  // make sure that the argument is of the correct type
  const DistSCMatrix* la
    = require_dynamic_cast<const DistSCMatrix*>(a,"DistSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "DistSCMatrix::accumulate(SCMatrix*a): "
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
               << "DistSCMatrixListSubblockIter::accumulate block mismatch: "
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

void
DistSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const DistSymmSCMatrix* la
    = require_dynamic_cast<const DistSymmSCMatrix*>(a,"DistSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "DistSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  Ref<SCMatrixSubblockIter> I = ((SymmSCMatrix*)a)->all_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixBlock> block = I->block();
      if (DEBUG)
          ExEnv::outn() << messagegrp()->me() << ": "
               << block->class_name()
               << "(" << block->blocki << ", " << block->blockj << ")"
               << endl;
      // see if i've got this block
      Ref<SCMatrixBlock> localblock
          = block_to_block(block->blocki,block->blockj);
      if (localblock) {
          // the diagonal blocks require special handling
          if (block->blocki == block->blockj) {
              int n = rowblocks()->size(block->blocki);
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<n; i++) {
                  for (int j=0; j<i; j++) {
                      double tmp = *dat1;
                      dat2[i*n+j] += tmp;
                      dat2[j*n+i] += tmp;
                      dat1++;
                    }
                  dat2[i*n+i] = *dat1++;
                }
            }
          else {
              int n = block->ndat();
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<n; i++) {
                  dat2[i] += dat1[i];
                }
            }
        }
      // now for the transpose
      if (block->blocki != block->blockj) {
          localblock = block_to_block(block->blockj,block->blocki);
          if (localblock) {
              int nr = rowblocks()->size(block->blocki);
              int nc = rowblocks()->size(block->blockj);
              double *dat1 = block->dat();
              double *dat2 = localblock->dat();
              for (int i=0; i<nr; i++) {
                  for (int j=0; j<nc; j++) {
                      dat2[j*nr+i] += *dat1++;
                    }
                }
            }
        }
    }
}

void
DistSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const DistDiagSCMatrix* la
    = require_dynamic_cast<const DistDiagSCMatrix*>(a,"DistSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "DistSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  Ref<SCMatrixSubblockIter> I = ((DiagSCMatrix*)a)->all_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixBlock> block = I->block();
      // see if i've got this block
      Ref<SCMatrixBlock> localblock
          = block_to_block(block->blocki,block->blockj);
      if (localblock) {
          int n = rowblocks()->size(block->blocki);
          double *dat1 = block->dat();
          double *dat2 = localblock->dat();
          for (int i=0; i<n; i++) {
              dat2[i*n+i] += *dat1++;
            }
        }
    }
}

void
DistSCMatrix::accumulate(const SCVector*a)
{
  // make sure that the argument is of the correct type
  const DistSCVector* la
    = require_dynamic_cast<const DistSCVector*>(a,"DistSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
      ExEnv::errn() << indent << "DistSCMatrix::accumulate(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  SCMatrixBlockListIter I, J;
  for (I = blocklist->begin(), J = la->blocklist->begin();
       I != blocklist->end() && J != la->blocklist->end();
       I++,J++) {
      int n = I.block()->ndat();
      if (n != J.block()->ndat()) {
          ExEnv::errn() << indent << "DistSCMatrix::accumulate(SCVector*a): "
               << "block lists do not match" << endl;
          abort();
        }
      double *dati = I.block()->dat();
      double *datj = J.block()->dat();
      for (int i=0; i<n; i++) dati[i] += datj[i];
    }
}

void
DistSCMatrix::accumulate_product_rr(SCMatrix*pa,SCMatrix*pb)
{
  const char* name = "DistSCMatrix::accumulate_product_rr";
  // make sure that the arguments are of the correct type
  DistSCMatrix* a = require_dynamic_cast<DistSCMatrix*>(pa,name);
  DistSCMatrix* b = require_dynamic_cast<DistSCMatrix*>(pb,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(a->rowdim()) || !coldim()->equiv(b->coldim()) ||
      !a->coldim()->equiv(b->rowdim())) {
      ExEnv::errn() << indent
           << "DistSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b): "
           << "dimensions don't match\n";
      ExEnv::err0() << indent << "rowdim():" << endl;
      rowdim().print();
      ExEnv::err0() << indent << "coldim():" << endl;
      coldim().print();
      ExEnv::err0() << indent << "a->rowdim():" << endl;
      a->rowdim().print();
      ExEnv::err0() << indent << "a->coldim():" << endl;
      a->coldim().print();
      ExEnv::err0() << indent << "b->rowdim():" << endl;
      b->rowdim().print();
      ExEnv::err0() << indent << "b->coldim():" << endl;
      b->coldim().print();
      abort();
    }

  // i need this in row form and a in row form
  create_vecform(Row);
  vecform_zero();
  a->create_vecform(Row);
  a->vecform_op(CopyToVec);

  Ref<SCMatrixSubblockIter> I = b->all_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixRectBlock> blk
          = dynamic_cast<SCMatrixRectBlock*>(I->block());
      int kk,k,jj,j,i,nj;
      nj = blk->jend - blk->jstart;
      double *data = blk->data;
      for (i=0; i<nvec; i++) {
          double *veci = vec[i];
          double *aveci = a->vec[i];
          for (j=blk->jstart,jj=0; j<blk->jend; j++,jj++) {
              double tmp = 0.0;
              for (k=blk->istart,kk=0; k<blk->iend; k++,kk++) {
                  tmp += aveci[k] * data[kk*nj+jj];
                }
              veci[j] += tmp;
            }
        }
    }

  vecform_op(AccumFromVec);
  delete_vecform();
  a->delete_vecform();
}

void
DistSCMatrix::create_vecform(Form f, int nvectors)
{
  // determine with rows/cols go on this node
  form = f;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  int n1=0, n2=0;
  if (form == Row) { n1 = nrow(); n2 = ncol(); }
  if (form == Col) { n1 = ncol(); n2 = nrow(); }
  if (nvectors == -1) {
      nvec = n1/nproc;
      vecoff = nvec*me;
      int nremain = n1%nproc;
      if (me < nremain) {
          vecoff += me;
          nvec++;
        }
      else {
          vecoff += nremain;
        }
    }
  else {
      nvec = nvectors;
      vecoff = 0;
    }

  // allocate storage
  vec = new double*[nvec];
  vec[0] = allocate<double>(nvec*n2);
  int i;
  for (i=1; i<nvec; i++) {
      vec[i] = &vec[0][i*n2];
    }
}

void
DistSCMatrix::vecform_zero()
{
  int n=0;
  if (form == Row) { n = nvec * ncol(); }
  if (form == Col) { n = nvec * nrow(); }

  double *dat = vec[0];
  for (int i=0; i<n; i++) dat[i] = 0.0;
}

void
DistSCMatrix::delete_vecform()
{
  deallocate(vec[0]);
  delete[] vec;
  vec = 0;
  nvec = 0;
}

void
DistSCMatrix::vecform_op(VecOp op, int *ivec)
{
  Ref<SCMatrixSubblockIter> i;
  if (op == CopyToVec || op == AccumToVec) {
      i = all_blocks(SCMatrixSubblockIter::Read);
    }
  else {
      if (op == CopyFromVec) assign(0.0);
      i = all_blocks(SCMatrixSubblockIter::Accum);
    }
  for (i->begin(); i->ready(); i->next()) {
      Ref<SCMatrixRectBlock> b = dynamic_cast<SCMatrixRectBlock*>(i->block());
      if (DEBUG)
          ExEnv::outn() << messagegrp()->me() << ": "
               << "got block " << b->blocki << ' ' << b->blockj << endl;
      int b1start, b2start, b1end, b2end;
      if (form == Row) {
          b1start = b->istart;
          b2start = b->jstart;
          b1end = b->iend;
          b2end = b->jend;
        }
      else {
          b1start = b->jstart;
          b2start = b->istart;
          b1end = b->jend;
          b2end = b->iend;
        }
      int nbj = b->jend - b->jstart;
      int start, end;
      if (ivec) {
          start = b1start;
          end = b1end;
        }
      else {
          start = b1start > vecoff ? b1start : vecoff;
          end = b1end > vecoff+nvec ? vecoff+nvec : b1end;
        }
      double *dat = b->data;
      int off = -b1start;
      for (int j=start; j<end; j++) {
          double *vecj;
          if (ivec) {
              vecj = 0;
              for (int ii=0; ii<nvec; ii++) {
                  if (ivec[ii] == j) { vecj = vec[ii]; break; }
                }
              if (!vecj) continue;
              if (DEBUG)
                  ExEnv::outn() << messagegrp()->me() << ": getting ["
                               << j << ","
                               << b2start << "-" << b2end << ")" << endl;
            }
          else {
              vecj = vec[j-vecoff];
            }
          for (int k=b2start; k<b2end; k++) {
              int blockoffset;
              if (DEBUG)
                  ExEnv::outn() << messagegrp()->me() << ": "
                       << "using vec[" << j-vecoff << "]"
                       << "[" << k << "]" << endl;
              if (form == Row) {
                  blockoffset = (j+off)*nbj+k - b2start;
                  if (DEBUG)
                      ExEnv::outn() << messagegrp()->me() << ": "
                           << "Row datum offset is "
                           << "(" << j << "+" << off << ")*" << nbj << "+" << k
                           << "-" << b2start
                           << " = " << blockoffset << "(" << b->ndat() << ") "
                           << " -> " << dat[blockoffset] << endl;
                }
              else {
                  blockoffset = (k-b2start)*nbj+j+off;
                }
              if (blockoffset >= b->ndat()) {
                  fail("bad offset");
                }
              double *datum = &dat[blockoffset];
              if (op == CopyToVec) {
                  if (DEBUG)
                      ExEnv::outn() << messagegrp()->me() << ": "
                           << "copying " << *datum << " "
                           << "to " << j << " " << k << endl;
                  vecj[k] = *datum;
                }
              else if (op == CopyFromVec) {
                  *datum = vecj[k];
                }
              else if (op == AccumToVec) {
                  vecj[k] += *datum;
                }
              else if (op == AccumFromVec) {
                  *datum += vecj[k];
                }
            }
        }
    }
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
DistSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "DistSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  DistSCVector* la = require_dynamic_cast<DistSCVector*>(a,name);
  DistSCVector* lb = require_dynamic_cast<DistSCVector*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent
           << "DistSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  Ref<SCMatrixSubblockIter> I = a->all_blocks(SCMatrixSubblockIter::Read);
  Ref<SCMatrixSubblockIter> J = b->all_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCVectorSimpleBlock> vi
          = dynamic_cast<SCVectorSimpleBlock*>(I->block());
      int ni = vi->iend - vi->istart;
      for (J->begin(); J->ready(); J->next()) {
          Ref<SCVectorSimpleBlock> vj
              = dynamic_cast<SCVectorSimpleBlock*>(J->block());
          Ref<SCMatrixRectBlock> rij
              = dynamic_cast<SCMatrixRectBlock*>(block_to_block(vi->blocki,
                                                       vj->blocki).pointer());
          // if the block is held locally sum in the outer prod contrib
          if (rij) {
              int nj = vj->iend - vj->istart;
              double *dat = rij->data;
              for (int i=0; i<ni; i++) {
                  for (int j=0; j<nj; j++) {
                      *dat += vi->data[i]*vj->data[j];
                    }
                }
            }
        }
    }
}

void
DistSCMatrix::transpose_this()
{
  RefSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;

  Ref<SCMatrixBlockList> oldlist = blocklist;
  init_blocklist();

  assign(0.0);

  Ref<SCMatrixSubblockIter> I
      = new DistSCMatrixListSubblockIter(SCMatrixSubblockIter::Read,
                                         oldlist,
                                         messagegrp());
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixRectBlock> remote
          = dynamic_cast<SCMatrixRectBlock*>(I->block());
      Ref<SCMatrixRectBlock> local
          = dynamic_cast<SCMatrixRectBlock*>(block_to_block(remote->blockj,
                                                    remote->blocki).pointer());
      if (local) {
          int ni = local->iend - local->istart;
          int nj = local->jend - local->jstart;
          for (int i=0; i<ni; i++) {
              for (int j=0; j<nj; j++) {
                  local->data[i*nj+j] = remote->data[j*ni+i];
                }
            }
        }
    }

}

double
DistSCMatrix::invert_this()
{
  if (nrow() != ncol()) {
      ExEnv::errn() << indent << "DistSCMatrix::invert_this: matrix is not square\n";
      abort();
    }
  RefSymmSCMatrix refs = kit()->symmmatrix(d1);
  refs->assign(0.0);
  refs->accumulate_symmetric_product(this);
  double determ2 = refs->invert_this();
  transpose_this();
  RefSCMatrix reft = copy();
  assign(0.0);
  ((SCMatrix*)this)->accumulate_product(reft.pointer(), refs.pointer());
  return sqrt(fabs(determ2));
}

void
DistSCMatrix::gen_invert_this(double condition_number_threshold)
{
  invert_this();
}

double
DistSCMatrix::determ_this()
{
  if (nrow() != ncol()) {
    ExEnv::errn() << indent << "DistSCMatrix::determ_this: matrix is not square\n";
    abort();
  }
  return invert_this();
}

double
DistSCMatrix::trace()
{
  if (nrow() != ncol()) {
    ExEnv::errn() << indent << "DistSCMatrix::trace: matrix is not square\n";
    abort();
  }

  double ret=0.0;
  Ref<SCMatrixSubblockIter> I = local_blocks(SCMatrixSubblockIter::Read);
  for (I->begin(); I->ready(); I->next()) {
      Ref<SCMatrixRectBlock> b = dynamic_cast<SCMatrixRectBlock*>(I->block());
      if (b->blocki == b->blockj) {
          int ni = b->iend-b->istart;
          for (int i=0; i<ni; i++) {
              ret += b->data[i*ni+i];
            }
        }
    }
  messagegrp()->sum(ret);

  return ret;
}

double
DistSCMatrix::solve_this(SCVector*v)
{
  error("no solve_this");

  // make sure that the dimensions match
  if (!rowdim()->equiv(v->dim())) {
      ExEnv::errn() << indent << "DistSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match\n";
      abort();
    }

  return 0.0;
}

void
DistSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  error("no schmidt_orthog");
}

int
DistSCMatrix::schmidt_orthog_tol(SymmSCMatrix *S, double tol, double *res)
{
  error("no schmidt_orthog_tol");
  return 0;
}

void
DistSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
//       ExEnv::outn() << "rect elemop processing a block of type "
//            << i.block()->class_name() << endl;
      op->process_base(i.block());
    }
  if (op->has_collect()) op->collect(messagegrp());
}

void
DistSCMatrix::element_op(const Ref<SCElementOp2>& op,
                          SCMatrix* m)
{
  DistSCMatrix *lm
      = require_dynamic_cast<DistSCMatrix*>(m,"DistSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
      ExEnv::errn() << indent << "DistSCMatrix: bad element_op\n";
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
DistSCMatrix::element_op(const Ref<SCElementOp3>& op,
                          SCMatrix* m,SCMatrix* n)
{
  DistSCMatrix *lm
      = require_dynamic_cast<DistSCMatrix*>(m,"DistSCMatrix::element_op");
  DistSCMatrix *ln
      = require_dynamic_cast<DistSCMatrix*>(n,"DistSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
      ExEnv::errn() << indent << "DistSCMatrix: bad element_op\n";
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
DistSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  // cast so the non const vprint member can be called
  ((DistSCMatrix*)this)->vprint(title,os,prec);
}

void
DistSCMatrix::vprint(const char *title, ostream& os, int prec)
{
  int i,j;
  int lwidth;
  double max=this->maxabs();

  int me = messagegrp()->me();

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec+5+(int) max;

  os.setf(ios::fixed,ios::floatfield); os.precision(prec);
  os.setf(ios::right,ios::adjustfield);

  if (messagegrp()->me() == 0) {
      if (title) os << endl << indent << title << endl;
      else os << endl;
    }

  if (nrow()==0 || ncol()==0) {
      if (me == 0) os << indent << "empty matrix\n";
      return;
    }

  create_vecform(Row);
  vecform_op(CopyToVec);

  int nc = ncol();

  int tmp = 0;
  if (me != 0) {
      messagegrp()->recv(me-1, tmp);
    }
  else {
      os << indent;
      for (i=0; i<nc; i++) os << setw(lwidth) << i;
      os << endl;
    }
  for (i=0; i<nvec; i++) {
      os << indent << setw(5) << i+vecoff;
      for (j=0; j<nc; j++) os << setw(lwidth) << vec[i][j];
      os << endl;
    }
  if (messagegrp()->n() > 1) {
      // send the go ahead to the next node
      int dest = me+1;
      if (dest == messagegrp()->n()) dest = 0;
      messagegrp()->send(dest, tmp);
      // make node zero wait on the last node
      if (me == 0) messagegrp()->recv(messagegrp()->n()-1, tmp);
    }

  delete_vecform();
}

Ref<SCMatrixSubblockIter>
DistSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new SCMatrixListSubblockIter(access, blocklist);
}

Ref<SCMatrixSubblockIter>
DistSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  return new DistSCMatrixListSubblockIter(access, blocklist, messagegrp());
}

void
DistSCMatrix::error(const char *msg)
{
  ExEnv::errn() << "DistSCMatrix: error: " << msg << endl;
}

Ref<DistSCMatrixKit>
DistSCMatrix::skit()
{
  return dynamic_cast<DistSCMatrixKit*>(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
