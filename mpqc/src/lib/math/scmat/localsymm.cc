//
// localsymm.cc
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
#include <math/scmat/offset.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////
// LocalSymmSCMatrix member functions

#define CLASSNAME LocalSymmSCMatrix
#define PARENTS public SymmSCMatrix
#include <util/class/classi.h>
void *
LocalSymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SymmSCMatrix::_castdown(cd);
  return do_castdowns(casts,cd);
}

static double **
init_symm_rows(double *data, int n)
{
  double** r = new double*[n];
  for (int i=0; i<n; i++) r[i] = &data[(i*(i+1))/2];
  return r;
}

LocalSymmSCMatrix::LocalSymmSCMatrix(const RefSCDimension&a,
                                     LocalSCMatrixKit *kit):
  SymmSCMatrix(a,kit),
  rows(0)
{
  resize(a->n());
}

LocalSymmSCMatrix::~LocalSymmSCMatrix()
{
  if (rows) delete[] rows;
}

int
LocalSymmSCMatrix::compute_offset(int i,int j) const
{
  if (i<0 || j<0 || i>=d->n() || j>=d->n()) {
      ExEnv::err() << indent << "LocalSymmSCMatrix: index out of bounds\n";
      abort();
    }
  return ij_offset(i,j);
}

void
LocalSymmSCMatrix::resize(int n)
{
  block = new SCMatrixLTriBlock(0,n);
  rows = init_symm_rows(block->data,n);
}

double *
LocalSymmSCMatrix::get_data()
{
  return block->data;
}

double **
LocalSymmSCMatrix::get_rows()
{
  return rows;
}

double
LocalSymmSCMatrix::get_element(int i,int j) const
{
  return block->data[compute_offset(i,j)];
}

void
LocalSymmSCMatrix::set_element(int i,int j,double a)
{
  block->data[compute_offset(i,j)] = a;
}

void
LocalSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  block->data[compute_offset(i,j)] += a;
}

SCMatrix *
LocalSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::get_subblock: trying to get too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  RefSCDimension dnrow = (nsrow==n()) ? dim().pointer():new SCDimension(nsrow);
  RefSCDimension dncol = (nscol==n()) ? dim().pointer():new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  LocalSCMatrix *lsb =
    LocalSCMatrix::require_castdown(sb, "LocalSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = get_element(i+br,j+bc);
      
  return sb;
}

SymmSCMatrix *
LocalSymmSCMatrix::get_subblock(int br, int er)
{
  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::get_subblock: trying to get too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  RefSCDimension dnrow = new SCDimension(nsrow);

  SymmSCMatrix * sb = kit()->symmmatrix(dnrow);
  sb->assign(0.0);

  LocalSymmSCMatrix *lsb =
    LocalSymmSCMatrix::require_castdown(sb, "LocalSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      lsb->rows[i][j] = get_element(i+br,j+br);
      
  return sb;
}

void
LocalSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  LocalSCMatrix *lsb =
    LocalSCMatrix::require_castdown(sb, "LocalSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::assign_subblock: trying to assign too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+bc,lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  LocalSymmSCMatrix *lsb = LocalSymmSCMatrix::require_castdown(sb,
                                        "LocalSymmSCMatrix::assign_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::assign_subblock: trying to assign too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::accumulate_subblock: trying to "
         << "accumulate too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

void
LocalSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  LocalSCMatrix *lsb = LocalSCMatrix::require_castdown(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::accumulate_subblock: trying to "
         << "accumulate too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

SCVector *
LocalSymmSCMatrix::get_row(int i)
{
  if (i >= n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::get_row: trying to get invalid row "
         << i << " max " << n() << endl;
    abort();
  }
  
  SCVector * v = kit()->vector(dim());

  LocalSCVector *lv =
    LocalSCVector::require_castdown(v, "LocalSymmSCMatrix::get_row");

  for (int j=0; j < n(); j++)
    lv->set_element(j,get_element(i,j));
      
  return v;
}

void
LocalSymmSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::assign_row: trying to assign invalid row "
         << i << " max " << n() << endl;
    abort();
  }
  
  if (v->n() != n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::assign_row: vector is wrong size "
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }
  
  LocalSCVector *lv =
    LocalSCVector::require_castdown(v, "LocalSymmSCMatrix::assign_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::accumulate_row: trying to "
         << "accumulate invalid row "
         << i << " max " << n() << endl;
    abort();
  }
  
  if (v->n() != n()) {
    ExEnv::err() << indent
         << "LocalSymmSCMatrix::accumulate_row: vector is wrong size"
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }
  
  LocalSCVector *lv =
    LocalSCVector::require_castdown(v, "LocalSymmSCMatrix::accumulate_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,get_element(i,j)+lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const LocalSymmSCMatrix* la
    = LocalSymmSCMatrix::require_const_castdown(a,"LocalSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = (this->n() * (this->n() + 1))/2;
  for (int i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

double
LocalSymmSCMatrix::invert_this()
{
  return cmat_invert(rows,1,n());
}

double
LocalSymmSCMatrix::determ_this()
{
  return cmat_determ(rows,1,n());
}

double
LocalSymmSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < n(); i++) ret += rows[i][i];
  return ret;
}

double
LocalSymmSCMatrix::solve_this(SCVector*v)
{
  LocalSCVector* lv =
    LocalSCVector::require_castdown(v,"LocalSymmSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match\n";
      abort();
    }

  return cmat_solve_lin(rows,1,lv->block->data,n());
}

void
LocalSymmSCMatrix::gen_invert_this()
{
  if (n() == 0) return;

  double *evals = new double[n()];
  double **evecs = cmat_new_square_matrix(n());
  
  cmat_diag(rows,evals,evecs,n(),1,1.0e-15);

  for (int i=0; i < n(); i++) {
    if (fabs(evals[i]) > 1.0e-8)
      evals[i] = 1.0/evals[i];
    else
      evals[i] = 0;
  }

  cmat_transform_diagonal_matrix(rows, n(), evals, n(), evecs, 0);
  
  delete[] evals;
  cmat_delete_matrix(evecs);  
}

void
LocalSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  if (n() == 0) return;

  const char* name = "LocalSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  LocalDiagSCMatrix* la = LocalDiagSCMatrix::require_castdown(a,name);
  LocalSCMatrix* lb = LocalSCMatrix::require_castdown(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::"
           << "diagonalize(DiagSCMatrix*a,SCMatrix*b): bad dims";
      abort();
    }

  double *eigvals;
  double **eigvecs;
  if (!la) {
      eigvals = new double[n()];
    }
  else {
      eigvals = la->block->data;
    }

  if (!lb) {
      eigvecs = cmat_new_square_matrix(n());
    }
  else {
      eigvecs = lb->rows;
    }

  cmat_diag(rows,eigvals,eigvecs,n(),1,1.0e-15);

  if (!la) delete[] eigvals;
  if (!lb) cmat_delete_matrix(eigvecs);
}

// computes this += a * a.t
void
LocalSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::"
           << "accumulate_symmetric_product(SCMatrix*a): bad dim";
      abort();
    }

  cmat_symmetric_mxm(rows,n(),la->rows,la->ncol(),1);
}

// computes this += a + a.t
void
LocalSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = LocalSCMatrix::require_castdown(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::"
           << "accumulate_symmetric_sum(SCMatrix*a): bad dim";
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double** adat = la->rows;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i][j] + adat[j][i];
        }
    }
}

void
LocalSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::"
           << "accumulate_symmetric_outer_product(SCMatrix*a): bad dim";
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double* adat = la->block->data;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i]*adat[j];
        }
    }
}

// this += a * b * transpose(a)
void
LocalSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b,
                                        SCMatrix::Transform t)
{
  // do the necessary castdowns
  LocalSCMatrix*la
    = LocalSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  LocalSymmSCMatrix*lb = require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::accumulate_transform: bad dim" << endl;
      ExEnv::err() << indent << "this dimension:" << endl << incindent;
      dim()->print(ExEnv::err());
      ExEnv::err() << decindent << indent
           << "a row and col dimension:" << endl << incindent;
      a->rowdim()->print(ExEnv::err());
      a->coldim()->print(ExEnv::err());
      ExEnv::err() << decindent << indent << "b dimension:" << endl << incindent;
      b->dim()->print(ExEnv::err());
      ExEnv::err() << decindent;
      abort();
    }

  cmat_transform_symmetric_matrix(rows,n(),lb->rows,lb->n(),la->rows,1);
}

// this += a * b * transpose(a)
void
LocalSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b,
                                        SCMatrix::Transform t)
{
  // do the necessary castdowns
  LocalSCMatrix*la
    = LocalSCMatrix::require_castdown(a,"%s::accumulate_transform",
                                      class_name());
  LocalDiagSCMatrix*lb
    = LocalDiagSCMatrix::require_castdown(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::accumulate_transform: bad dim\n";
      abort();
    }

  cmat_transform_diagonal_matrix(rows,n(),lb->block->data,lb->n(),la->rows,1);
}

void
LocalSymmSCMatrix::accumulate_transform(SymmSCMatrix*a,SymmSCMatrix*b)
{
  SymmSCMatrix::accumulate_transform(a,b);
}

double
LocalSymmSCMatrix::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  LocalSCVector* la
    = LocalSCVector::require_castdown(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::err() << indent
           << "LocalSCVector::scalar_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = n();
  double* adat = la->block->data;
  double result = 0.0;
  for (int i=0; i<nelem; i++) {
      for (int j=0; j<i; j++) {
          result += 2.0 * rows[i][j] * adat[i] * adat[j];
        }
      result += rows[i][i] * adat[i] * adat[i];
    }
  return result;
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp& op)
{
  op->process_spec_ltri(block.pointer());
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp2& op,
                              SymmSCMatrix* m)
{
  LocalSymmSCMatrix *lm
      = LocalSymmSCMatrix::require_castdown(m,"LocalSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::err() << indent << "LocalSymmSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_ltri(block.pointer(), lm->block.pointer());
}

void
LocalSymmSCMatrix::element_op(const RefSCElementOp3& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  LocalSymmSCMatrix *lm
      = LocalSymmSCMatrix::require_castdown(m,"LocalSymSCMatrix::element_op");
  LocalSymmSCMatrix *ln
      = LocalSymmSCMatrix::require_castdown(n,"LocalSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::err() << indent << "LocalSymmSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_ltri(block.pointer(),
                        lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH (with a bit of hacking)
void
LocalSymmSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int ii,jj,kk,nn;
  int i,j;
  int lwidth,width;
  double max=this->maxabs();

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec + 5 + (int) max;
  width = 75/(lwidth+SCFormIO::getindent(os));

  if (title)
    os << node0 << endl << indent << title << endl;
  else
    os << node0 << endl;

  if (n()==0) {
    os << node0 << indent << "empty matrix\n";
    return;
  }

  for (ii=jj=0;;) {
    ii++; jj++;
    kk=width*jj;
    nn = (n() > kk) ? kk : n();

    // print column indices
    os << node0 << indent;
    for (i=ii; i <= nn; i++)
      os << node0 << scprintf("%*d",lwidth,i);
    os << node0 << endl;

    // print the rows
    for (i=ii-1; i < n() ; i++) {
      os << node0 << indent << scprintf("%5d",i+1);
      for (j=ii-1; j<nn && j<=i; j++)
        os << node0 << scprintf("%*.*f",lwidth,prec,rows[i][j]);
      os << node0 << endl;
    }
    os << node0 << endl;

    if (n() <= kk) {
      os.flush();
      return;
    }
    ii=kk;
  }
}

RefSCMatrixSubblockIter
LocalSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  if (messagegrp()->n() > 1) {
      ExEnv::err() << indent
           << "LocalSymmSCMatrix::local_blocks: not valid for local matrices"
           << endl;
      abort();
    }
  RefSCMatrixSubblockIter iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

RefSCMatrixSubblockIter
LocalSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::err() << indent << "LocalSymmSCMatrix::all_blocks: "
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
