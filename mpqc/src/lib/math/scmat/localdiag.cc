//
// localdiag.cc
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
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// LocalDiagSCMatrix member functions

#define CLASSNAME LocalDiagSCMatrix
#define PARENTS public DiagSCMatrix
#include <util/class/classi.h>
void *
LocalDiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DiagSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalDiagSCMatrix::LocalDiagSCMatrix(const RefSCDimension&a,
                                     LocalSCMatrixKit *kit):
  DiagSCMatrix(a,kit)
{
  resize(a->n());
}

LocalDiagSCMatrix::~LocalDiagSCMatrix()
{
}

void
LocalDiagSCMatrix::resize(int n)
{
  block = new SCMatrixDiagBlock(0,n);
}

double *
LocalDiagSCMatrix::get_data()
{
  return block->data;
}

double
LocalDiagSCMatrix::get_element(int i) const
{
  return block->data[i];
}

void
LocalDiagSCMatrix::set_element(int i,double a)
{
  block->data[i] = a;
}

void
LocalDiagSCMatrix::accumulate_element(int i,double a)
{
  block->data[i] += a;
}

void
LocalDiagSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the argument is of the correct type
  const LocalDiagSCMatrix* la
    = LocalDiagSCMatrix::require_castdown(a,"LocalDiagSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "LocalDiagSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = n();
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

double
LocalDiagSCMatrix::invert_this()
{
  double det = 1.0;
  int nelem = n();
  double*data = block->data;
  for (int i=0; i<nelem; i++) {
      det *= data[i];
      data[i] = 1.0/data[i];
    }
  return det;
}

double
LocalDiagSCMatrix::determ_this()
{
  double det = 1.0;
  int nelem = n();
  double *data = block->data;
  for (int i=0; i < nelem; i++) {
    det *= data[i];
  }
  return det;
}

double
LocalDiagSCMatrix::trace()
{
  double det = 0;
  int nelem = n();
  double *data = block->data;
  for (int i=0; i < nelem; i++) {
    det += data[i];
  }
  return det;
}

void
LocalDiagSCMatrix::gen_invert_this()
{
  int nelem = n();
  double *data = block->data;
  for (int i=0; i < nelem; i++) {
    if (fabs(data[i]) > 1.0e-8)
      data[i] = 1.0/data[i];
    else
      data[i] = 0;
  }
}

void
LocalDiagSCMatrix::element_op(const RefSCElementOp& op)
{
  op->process_spec_diag(block.pointer());
}

void
LocalDiagSCMatrix::element_op(const RefSCElementOp2& op,
                              DiagSCMatrix* m)
{
  LocalDiagSCMatrix *lm
      = LocalDiagSCMatrix::require_castdown(m,"LocalDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      cerr << indent << "LocalDiagSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_diag(block.pointer(), lm->block.pointer());
}

void
LocalDiagSCMatrix::element_op(const RefSCElementOp3& op,
                              DiagSCMatrix* m,DiagSCMatrix* n)
{
  LocalDiagSCMatrix *lm
      = LocalDiagSCMatrix::require_castdown(m,"LocalDiagSCMatrix::element_op");
  LocalDiagSCMatrix *ln
      = LocalDiagSCMatrix::require_castdown(n,"LocalDiagSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      cerr << indent << "LocalDiagSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_diag(block.pointer(),
                        lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalDiagSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int i;
  int lwidth;
  double max=this->maxabs();

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec + 5 + (int) max;

  if (title)
    os << node0 << endl << indent << title << endl;
  else
    os << node0 << endl;

  if (n()==0) {
    os << node0 << indent << "empty matrix\n";
    return;
  }

  for (i=0; i<n(); i++)
    os << node0 << indent
       << scprintf("%5d %*.*f\n",i+1,lwidth,prec,block->data[i]);
  os << node0 << endl;

  os.flush();
}

RefSCMatrixSubblockIter
LocalDiagSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  if (messagegrp()->n() > 1) {
      cerr << indent
           << "LocalDiagSCMatrix::local_blocks: not valid for local matrices"
           << endl;
      abort();
    }
  RefSCMatrixSubblockIter iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

RefSCMatrixSubblockIter
LocalDiagSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << indent << "LocalDiagSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  return local_blocks(access);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
