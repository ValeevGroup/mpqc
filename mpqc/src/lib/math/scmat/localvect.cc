//
// localvect.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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
// LocalSCVector member functions

#define CLASSNAME LocalSCVector
#define PARENTS public SCVector
#include <util/class/classi.h>
void *
LocalSCVector::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCVector::_castdown(cd);
  return do_castdowns(casts,cd);
}

LocalSCVector::LocalSCVector(const RefSCDimension&a,
                             LocalSCMatrixKit *kit):
  SCVector(a,kit)
{
  resize(a->n());
}

void
LocalSCVector::resize(int n)
{
  block = new SCVectorSimpleBlock(0,n);
}

LocalSCVector::~LocalSCVector()
{
}

double *
LocalSCVector::get_data()
{
  return block->data;
}

double
LocalSCVector::get_element(int i)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      cerr << indent << "LocalSCVector::get_element: bad index\n";
      abort();
    }
  return block->data[i];
}

void
LocalSCVector::set_element(int i,double a)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      cerr << indent << "LocalSCVector::set_element: bad index\n";
      abort();
    }
  block->data[i] = a;
}

void
LocalSCVector::accumulate_element(int i,double a)
{
  int size = block->iend - block->istart;
  if (i < 0 || i >= size) {
      cerr << indent << "LocalSCVector::accumulate_element: bad index\n";
      abort();
    }
  block->data[i] += a;
}

void
LocalSCVector::accumulate_product(SCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = LocalSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->rowdim()) || !la->coldim()->equiv(lb->dim())) {
      cerr << indent
           << "LocalSCVector:: accumulate_product(SCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  cmat_mxm(la->rows, 0,
           &(lb->block->data), 1,
           &(block->data), 1,
           n(), la->ncol(), 1,
           1);
}

void
LocalSCVector::accumulate_product(SymmSCMatrix*a,SCVector*b)
{
  const char* name = "LocalSCVector::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSymmSCMatrix* la = LocalSymmSCMatrix::require_castdown(a,name);
  LocalSCVector* lb = LocalSCVector::require_castdown(b,name);

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim()) || !la->dim()->equiv(lb->dim())) {
      cerr << indent
           << "LocalSCVector:: accumulate_product(SymmSCMatrix*a,SCVector*b): "
           << "dimensions don't match\n";
      abort();
    }

  double* thisdat = block->data;
  double** adat = la->rows;
  double* bdat = lb->block->data;
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
      thisdat[i] += tmp;
    }
}

void
LocalSCVector::accumulate(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "LocalSCVector::accumulate(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCVector::accumulate(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (!((la->rowdim()->equiv(dim()) && la->coldim()->n() == 1)
        || (la->coldim()->equiv(dim()) && la->rowdim()->n() == 1))) {
      cerr << indent << "LocalSCVector::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCVector::assign(double a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = a;
}

void
LocalSCVector::assign(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::assign");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "LocalSCVector::assign(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = la->block->data[i];
}

void
LocalSCVector::assign(const double*a)
{
  int nelem = d->n();
  int i;
  for (i=0; i<nelem; i++) block->data[i] = a[i];
}

double
LocalSCVector::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      cerr << indent << "LocalSCVector::scalar_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = d->n();
  int i;
  double result = 0.0;
  for (i=0; i<nelem; i++) result += block->data[i] * la->block->data[i];
  return result;
}

void
LocalSCVector::element_op(const RefSCElementOp& op)
{
  op->process_spec_vsimp(block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp2& op,
                          SCVector* m)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");

  if (!dim()->equiv(lm->dim())) {
      cerr << indent << "LocalSCVector: bad element_op\n";
      abort();
    }
  op->process_spec_vsimp(block.pointer(), lm->block.pointer());
}

void
LocalSCVector::element_op(const RefSCElementOp3& op,
                          SCVector* m,SCVector* n)
{
  LocalSCVector *lm
      = LocalSCVector::require_castdown(m, "LocalSCVector::element_op");
  LocalSCVector *ln
      = LocalSCVector::require_castdown(n, "LocalSCVector::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      cerr << indent << "LocalSCVector: bad element_op\n";
      abort();
    }
  op->process_spec_vsimp(block.pointer(),
                         lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSCVector::vprint(const char *title, ostream& os, int prec)
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
    os << node0 << indent << "empty vector\n";
    return;
  }

  for (i=0; i<n(); i++)
    os << node0 << indent
       << scprintf("%5d %*.*f\n",i+1,lwidth,prec,block->data[i]);
  os << node0 << endl;

  os.flush();
}

RefSCMatrixSubblockIter
LocalSCVector::local_blocks(SCMatrixSubblockIter::Access access)
{
  if (messagegrp()->n() > 1) {
      cerr << indent
           << "LocalSCVector::local_blocks: not valid for local matrices"
           << endl;
      abort();
    }
  RefSCMatrixSubblockIter iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

RefSCMatrixSubblockIter
LocalSCVector::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      cerr << indent << "LocalVectorSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  return local_blocks(access);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
