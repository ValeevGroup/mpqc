//
// replvect.cc
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

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/misc/consumableresources.h>
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// ReplSCVector member functions

static ClassDesc ReplSCVector_cd(
  typeid(ReplSCVector),"ReplSCVector",1,"public SCVector",
  0, 0, 0);

ReplSCVector::ReplSCVector(const RefSCDimension&a,ReplSCMatrixKit*k):
  SCVector(a,k)
{
  vector = allocate<double>(a->n());
  init_blocklist();
}

void
ReplSCVector::before_elemop()
{
  // zero out the blocks not in my block list
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc == me) continue;
      memset(&vector[d->blocks()->start(i)], 0,
             sizeof(double)*(d->blocks()->fence(i)
                             - d->blocks()->start(i)));
    }
}

void
ReplSCVector::after_elemop()
{
  messagegrp()->sum(vector, d->n());
}

void
ReplSCVector::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc != me) continue;
      blocklist->insert(
          new SCVectorSimpleSubBlock(d->blocks()->start(i),
                                     d->blocks()->fence(i),
                                     d->blocks()->start(i),
                                     vector));
    }
}

ReplSCVector::~ReplSCVector()
{
  if (vector)
      deallocate(vector);
  vector=0;
}

double
ReplSCVector::get_element(int i) const
{
  return vector[i];
}

void
ReplSCVector::set_element(int i,double a)
{
  vector[i] = a;
}

void
ReplSCVector::accumulate_element(int i,double a)
{
  vector[i] += a;
}

void
ReplSCVector::accumulate_product_rv(SCMatrix*a,SCVector*b)
{
  const char* name = "ReplSCVector::accumulate_product_rv";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = require_dynamic_cast<ReplSCMatrix*>(a,name);
  ReplSCVector* lb = require_dynamic_cast<ReplSCVector*>(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
      ExEnv::out0() << indent << "dim():" << endl << incindent;
      dim().print();
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "la->rowdim():" << endl << incindent;
      la->rowdim().print();
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "la->coldim():" << endl << incindent;
      la->coldim().print();
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "lb->dim():" << endl << incindent;
      lb->dim().print();
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent
           << "ReplSCVector::accumulate_product_rv(SCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  cmat_mxm(la->rows, 0,
           &lb->vector, 1,
           &vector, 1,
           n(), la->ncol(), 1,
           1);
}

void
ReplSCVector::accumulate_product_sv(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "ReplSCVector::accumulate_product_sv";
  // make sure that the arguments are of the correct type
  ReplSymmSCMatrix* la = require_dynamic_cast<ReplSymmSCMatrix*>(a,name);
  ReplSCVector* lb = require_dynamic_cast<ReplSCVector*>(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
      ExEnv::errn() << indent
           << "ReplSCVector::accumulate_product_sv(SymmSCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  double** adat = la->rows;
  double* bdat = lb->vector;
  double tmp;
  int n = dim()->n();
  int i, j;
  for (i=0; i<n; i++) {
      tmp = 0.0;
      for (j=0; j<=i; j++) {
          tmp += adat[i][j] * bdat[j];
        }
      for (; j<n; j++) {
          tmp += adat[j][i] * bdat[j];
        }
      vector[i] += tmp;
    }
}

void
ReplSCVector::accumulate(const SCVector*a)
{
  // make sure that the argument is of the correct type
  const ReplSCVector* la
    = require_dynamic_cast<const ReplSCVector*>(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCVector::accumulate(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] += la->vector[i];
}

void
ReplSCVector::accumulate(const SCMatrix*a)
{
  // make sure that the argument is of the correct type
  const ReplSCMatrix *la
    = require_dynamic_cast<const ReplSCMatrix*>(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!((la->rowdim()->equiv(dim()) && la->coldim()->n() == 1)
        || (la->coldim()->equiv(dim()) && la->rowdim()->n() == 1))) {
      ExEnv::errn() << indent << "ReplSCVector::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] += la->matrix[i];
}

void
ReplSCVector::assign_val(double a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = a;
}

void
ReplSCVector::assign_v(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = require_dynamic_cast<ReplSCVector*>(a,"ReplSCVector::assign_v");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCVector::assign_v(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = la->vector[i];
}

void
ReplSCVector::assign_p(const double*a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) vector[i] = a[i];
}

double
ReplSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = require_dynamic_cast<ReplSCVector*>(a,"ReplSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCVector::scalar_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  double result = 0.0;
  for (i=0; i<nelem; i++) result += vector[i] * la->vector[i];
  return result;
}

void
ReplSCVector::element_op(const Ref<SCElementOp>& op)
{
  if (op->has_side_effects()) before_elemop();
  SCMatrixBlockListIter i;
  for (i = blocklist->begin(); i != blocklist->end(); i++) {
      op->process_base(i.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_collect()) op->collect(messagegrp());
}

void
ReplSCVector::element_op(const Ref<SCElementOp2>& op,
                          SCVector* m)
{
  ReplSCVector *lm
      = require_dynamic_cast<ReplSCVector*>(m, "ReplSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::errn() << indent << "ReplSCVector: bad element_op\n";
      abort();
    }

  if (op->has_side_effects()) before_elemop();
  if (op->has_side_effects_in_arg()) lm->before_elemop();
  SCMatrixBlockListIter i, j;
  for (i = blocklist->begin(), j = lm->blocklist->begin();
       i != blocklist->end();
       i++, j++) {
      op->process_base(i.block(), j.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg()) lm->after_elemop();
  if (op->has_collect()) op->collect(messagegrp());
}

void
ReplSCVector::element_op(const Ref<SCElementOp3>& op,
                          SCVector* m,SCVector* n)
{
  ReplSCVector *lm
      = require_dynamic_cast<ReplSCVector*>(m, "ReplSCVector::element_op");
  ReplSCVector *ln
      = require_dynamic_cast<ReplSCVector*>(n, "ReplSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::errn() << indent << "ReplSCVector: bad element_op\n";
      abort();
    }
  if (op->has_side_effects()) before_elemop();
  if (op->has_side_effects_in_arg1()) lm->before_elemop();
  if (op->has_side_effects_in_arg2()) ln->before_elemop();
  SCMatrixBlockListIter i, j, k;
  for (i = blocklist->begin(),
           j = lm->blocklist->begin(),
           k = ln->blocklist->begin();
       i != blocklist->end();
       i++, j++, k++) {
      op->process_base(i.block(), j.block(), k.block());
    }
  if (op->has_side_effects()) after_elemop();
  if (op->has_side_effects_in_arg1()) lm->after_elemop();
  if (op->has_side_effects_in_arg2()) ln->after_elemop();
  if (op->has_collect()) op->collect(messagegrp());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
ReplSCVector::vprint(const char *title, ostream& os, int prec) const
{
  int i;
  int lwidth;
  double max=this->maxabs();

  if (messagegrp()->me() != 0) return;

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec + 5 + (int) max;

  if (title)
    os << endl << indent << title << endl;
  else
    os << endl;

  if (n()==0) {
    os << indent << "empty vector\n";
    return;
  }

  for (i=0; i < n(); i++)
    os << indent
       << scprintf("%5d %*.*f\n",i+1,lwidth,prec,vector[i]);
  os << endl;

  os.flush();
}

Ref<SCMatrixSubblockIter>
ReplSCVector::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          vector, d->n());
}

Ref<SCMatrixSubblockIter>
ReplSCVector::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << indent << "ReplSCVector::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  Ref<SCMatrixBlockList> allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCVectorSimpleSubBlock(0, d->n(), 0, vector));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          vector, d->n());
}

Ref<ReplSCMatrixKit>
ReplSCVector::skit()
{
  return dynamic_cast<ReplSCMatrixKit*>(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
