//
// repldiag.cc
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
#include <math/scmat/repl.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// ReplDiagSCMatrix member functions

static ClassDesc ReplDiagSCMatrix_cd(
  typeid(ReplDiagSCMatrix),"ReplDiagSCMatrix",1,"public DiagSCMatrix",
  0, 0, 0);

ReplDiagSCMatrix::ReplDiagSCMatrix(const RefSCDimension&a,ReplSCMatrixKit*k):
  DiagSCMatrix(a,k)
{
  matrix = new double[a->n()];
  init_blocklist();
}

void
ReplDiagSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc == me) continue;
      memset(&matrix[d->blocks()->start(i)], 0,
             sizeof(double)*(d->blocks()->fence(i)
                             - d->blocks()->start(i)));
    }
}

void
ReplDiagSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d->n());
}

void
ReplDiagSCMatrix::init_blocklist()
{
  int i;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0; i<d->blocks()->nblock(); i++) {
      if (i%nproc != me) continue;
      blocklist->insert(
          new SCMatrixDiagSubBlock(d->blocks()->start(i),
                                   d->blocks()->fence(i),
                                   d->blocks()->start(i),
                                   matrix));
    }
}

ReplDiagSCMatrix::~ReplDiagSCMatrix()
{
  if (matrix)
      delete[] matrix;
  matrix=0;
}

double
ReplDiagSCMatrix::get_element(int i) const
{
  return matrix[i];
}

void
ReplDiagSCMatrix::set_element(int i,double a)
{
  matrix[i] = a;
}

void
ReplDiagSCMatrix::accumulate_element(int i,double a)
{
  matrix[i] += a;
}

void
ReplDiagSCMatrix::assign_val(double val)
{
  int n = d->n();
  for (int i=0; i<n; i++) matrix[i] = val;
}

void
ReplDiagSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const ReplDiagSCMatrix* la
    = require_dynamic_cast<const ReplDiagSCMatrix*>(a,"ReplDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplDiagSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = n();
  for (int i=0; i<nelem; i++) matrix[i] += la->matrix[i];
}

double
ReplDiagSCMatrix::invert_this()
{
  double det = 1.0;
  int nelem = n();
  for (int i=0; i<nelem; i++) {
      det *= matrix[i];
      matrix[i] = 1.0/matrix[i];
    }
  return det;
}

double
ReplDiagSCMatrix::determ_this()
{
  double det = 1.0;
  int nelem = n();
  for (int i=0; i < nelem; i++) {
    det *= matrix[i];
  }
  return det;
}

double
ReplDiagSCMatrix::trace()
{
  double tr = 0;
  int nelem = n();
  for (int i=0; i < nelem; i++) {
    tr += matrix[i];
  }
  return tr;
}

void
ReplDiagSCMatrix::gen_invert_this(double condition_number_threshold)
{
  int nelem = n();
  const double sigma_max = this->maxabs();
  const double sigma_min_threshold = sigma_max / condition_number_threshold;
  for (int i=0; i < nelem; i++) {
    if (fabs(matrix[i]) > sigma_min_threshold)
      matrix[i] = 1.0/matrix[i];
    else
      matrix[i] = 0;
  }
}

void
ReplDiagSCMatrix::element_op(const Ref<SCElementOp>& op)
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
ReplDiagSCMatrix::element_op(const Ref<SCElementOp2>& op,
                              DiagSCMatrix* m)
{
  ReplDiagSCMatrix *lm
      = require_dynamic_cast<ReplDiagSCMatrix*>(m,"ReplDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::errn() << indent << "ReplDiagSCMatrix: bad element_op\n";
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
ReplDiagSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  ReplDiagSCMatrix *lm
      = require_dynamic_cast<ReplDiagSCMatrix*>(m,"ReplDiagSCMatrix::element_op");
  ReplDiagSCMatrix *ln
      = require_dynamic_cast<ReplDiagSCMatrix*>(n,"ReplDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::errn() << indent << "ReplDiagSCMatrix: bad element_op\n";
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
ReplDiagSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int i;
  int lwidth;
  double max=this->maxabs();

  if (messagegrp()->me() != 0) return;

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec+5+(int) max;

  if (title)
    os << endl << indent << title << endl;
  else
    os << endl;

  if (n()==0) {
    os << indent << "empty matrix\n";
    return;
  }

  for (i=0; i<n(); i++)
    os << indent
       << scprintf("%5d %*.*f\n",i+1,lwidth,prec,matrix[i]);
  os << endl;

  os.flush();
}

Ref<SCMatrixSubblockIter>
ReplDiagSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          matrix, d->n());
}

Ref<SCMatrixSubblockIter>
ReplDiagSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << "ReplDiagSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  Ref<SCMatrixBlockList> allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCMatrixDiagSubBlock(0, d->n(), 0, matrix));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          matrix, d->n());
}

Ref<ReplSCMatrixKit>
ReplDiagSCMatrix::skit()
{
  return dynamic_cast<ReplSCMatrixKit*>(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
