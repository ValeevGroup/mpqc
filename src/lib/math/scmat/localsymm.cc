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
#include <algorithm>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>
#include <math/scmat/offset.h>
#include <math/scmat/predicate.h>

#include <math/scmat/mops.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// LocalSymmSCMatrix member functions

static ClassDesc LocalSymmSCMatrix_cd(
  typeid(LocalSymmSCMatrix),"LocalSymmSCMatrix",1,"public SymmSCMatrix",
  0, 0, 0);

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
      ExEnv::errn() << indent << "LocalSymmSCMatrix: index out of bounds\n";
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
    ExEnv::errn() << indent
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
    require_dynamic_cast<LocalSCMatrix*>(sb, "LocalSymmSCMatrix::get_subblock");

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
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::get_subblock: trying to get too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }

  RefSCDimension dnrow = new SCDimension(nsrow);

  SymmSCMatrix * sb = kit()->symmmatrix(dnrow);
  sb->assign(0.0);

  LocalSymmSCMatrix *lsb =
    require_dynamic_cast<LocalSymmSCMatrix*>(sb, "LocalSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      lsb->rows[i][j] = get_element(i+br,j+br);

  return sb;
}

void
LocalSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  LocalSCMatrix *lsb =
    require_dynamic_cast<LocalSCMatrix*>(sb, "LocalSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::errn() << indent
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
  LocalSymmSCMatrix *lsb = require_dynamic_cast<LocalSymmSCMatrix*>(sb,
                                        "LocalSymmSCMatrix::assign_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::errn() << indent
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
  LocalSCMatrix *lsb = require_dynamic_cast<LocalSCMatrix*>(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::errn() << indent
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
  LocalSCMatrix *lsb = require_dynamic_cast<LocalSCMatrix*>(sb,
                                  "LocalSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::errn() << indent
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
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::get_row: trying to get invalid row "
         << i << " max " << n() << endl;
    abort();
  }

  SCVector * v = kit()->vector(dim());

  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSymmSCMatrix::get_row");

  for (int j=0; j < n(); j++)
    lv->set_element(j,get_element(i,j));

  return v;
}

void
LocalSymmSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::assign_row: trying to assign invalid row "
         << i << " max " << n() << endl;
    abort();
  }

  if (v->n() != n()) {
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::assign_row: vector is wrong size "
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }

  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSymmSCMatrix::assign_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::accumulate_row: trying to "
         << "accumulate invalid row "
         << i << " max " << n() << endl;
    abort();
  }

  if (v->n() != n()) {
    ExEnv::errn() << indent
         << "LocalSymmSCMatrix::accumulate_row: vector is wrong size"
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }

  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSymmSCMatrix::accumulate_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,get_element(i,j)+lv->get_element(j));
}

void
LocalSymmSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const LocalSymmSCMatrix* la
    = require_dynamic_cast<const LocalSymmSCMatrix*>(a,"LocalSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent
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
    require_dynamic_cast<LocalSCVector*>(v,"LocalSymmSCMatrix::solve_this");

  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
      ExEnv::errn() << indent
           << "LocalSymmSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match\n";
      abort();
    }

  return cmat_solve_lin(rows,1,lv->block->data,n());
}

void
LocalSymmSCMatrix::gen_invert_this(double condition_number_threshold)
{
  if (n() == 0) return;

  double *evals = new double[n()];
  double **evecs = cmat_new_square_matrix(n());

  cmat_diag(rows,evals,evecs,n(),1,1.0e-15);
  const double sigma_max = * std::max_element(evals, evals+n(), fabs_less<double>());
  const double sigma_min_threshold = sigma_max / condition_number_threshold;
  for (int i=0; i < n(); i++) {
    if (fabs(evals[i]) > sigma_min_threshold)
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
  LocalDiagSCMatrix* la = require_dynamic_cast<LocalDiagSCMatrix*>(a,name);
  LocalSCMatrix* lb = require_dynamic_cast<LocalSCMatrix*>(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent
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

void
LocalSymmSCMatrix::eigensystem(SymmSCMatrix*s, DiagSCMatrix*a, SCMatrix*b)
{
  if (n() == 0) return;

  const char* name = "LocalSymmSCMatrix::eigensystem";
  // make sure that the arguments is of the correct type
  LocalSymmSCMatrix* ls = require_dynamic_cast<LocalSymmSCMatrix*>(s,name);
  LocalDiagSCMatrix* la = require_dynamic_cast<LocalDiagSCMatrix*>(a,name);
  LocalSCMatrix* lb = require_dynamic_cast<LocalSCMatrix*>(b,name);

  if (!dim()->equiv(ls->dim()) ||
      !dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) ||
      !dim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent
           << "LocalSymmSCMatrix::eigensystem: bad dims";
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

  cmat_eigensystem(this->rows, ls->rows, eigvals, eigvecs,n(), 1);

  if (!la) delete[] eigvals;
  if (!lb) cmat_delete_matrix(eigvecs);
}

// computes this += a * a.t
void
LocalSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  LocalSCMatrix* la
    = require_dynamic_cast<LocalSCMatrix*>(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
      ExEnv::errn() << indent
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
    = require_dynamic_cast<LocalSCMatrix*>(a,"LocalSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent
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
    = require_dynamic_cast<LocalSCVector*>(a,"LocalSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent
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
  int i,j,k;
  int ii,jj;
  int nc, nr;

  // do the necessary castdowns
  LocalSCMatrix*la
    = require_dynamic_cast<LocalSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  LocalSymmSCMatrix*lb = require_dynamic_cast<LocalSymmSCMatrix*>(
      b,"%s::accumulate_transform", class_name());

  // check the dimensions
  if (t == SCMatrix::NormalTransform) {
    if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "LocalSymmSCMatrix::accumulate_transform: bad dim\n";
      abort();
    }

    nc = lb->n();
    nr = la->nrow();
  } else {
    if (!dim()->equiv(la->coldim()) || !lb->dim()->equiv(la->rowdim())) {
      ExEnv::errn() << indent << "LocalSymmSCMatrix::accumulate_transform: bad dim\n";
      abort();
    }

    nc = lb->n();
    nr = la->ncol();
  }

  if (nr==0 || nc==0)
    return;

  int nproc = messagegrp()->n();

  double **ablock = cmat_new_square_matrix(D1);
  double **bblock = cmat_new_square_matrix(D1);
  double **cblock = cmat_new_square_matrix(D1);

  double **temp = cmat_new_rect_matrix(D1,nc);

  for (i=0; i < nr; i += D1) {
      int ni = nr-i;
      if (ni > D1) ni = D1;

      memset(temp[0], 0, sizeof(double)*D1*nc);

      for (j=0; j < nc; j+= D1) {
          int nj = nc-j;
          if (nj > D1) nj = D1;

          for (k=0; k < nc; k += D1) {

              int nk = nc-k;
              if (nk > D1) nk = D1;

              if (t == SCMatrix::NormalTransform)
                  copy_block(ablock, la->rows, i, ni, k, nk);
              else
                  copy_trans_block(ablock, la->rows, i, ni, k, nk);

              copy_sym_block(bblock, lb->rows, j, nj, k, nk);
              copy_block(cblock, temp, 0, ni, j, nj);
              mult_block(ablock, bblock, cblock, ni, nj, nk);
              return_block(temp, cblock, 0, ni, j, nj);
            }
        }

      // now do ab * a~
      for (j=0; j <= i; j+= D1) {
          int nj = nr-j;
          if (nj > D1) nj = D1;

          memset(cblock[0], 0, sizeof(double)*D1*D1);

          for (k=0; k < nc; k += D1) {

              int nk = nc-k;
              if (nk > D1) nk = D1;

              copy_block(ablock, temp, 0, ni, k, nk);
              if (t == SCMatrix::NormalTransform)
                  copy_block(bblock, la->rows, j, nj, k, nk);
              else
                  copy_trans_block(bblock, la->rows, j, nj, k, nk);

              mult_block(ablock, bblock, cblock, ni, nj, nk);
            }

          // copy cblock(i,j) into result
          if (j==i) {
              for (ii=0; ii < ni; ii++)
                  for (jj=0; jj <= ii; jj++)
                      rows[i+ii][j+jj] += cblock[ii][jj];
            } else {
                for (ii=0; ii < ni; ii++)
                    for (jj=0; jj < nj; jj++)
                        rows[i+ii][j+jj] += cblock[ii][jj];
              }
        }
    }

  cmat_delete_matrix(temp);

  cmat_delete_matrix(ablock);
  cmat_delete_matrix(bblock);
  cmat_delete_matrix(cblock);
}

// this += a * b * transpose(a)
void
LocalSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b,
                                        SCMatrix::Transform t)
{
  // do the necessary castdowns
  LocalSCMatrix*la
    = require_dynamic_cast<LocalSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  LocalDiagSCMatrix*lb
    = require_dynamic_cast<LocalDiagSCMatrix*>(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent
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
    = require_dynamic_cast<LocalSCVector*>(a,"LocalSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent
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
LocalSymmSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  op->process_spec_ltri(block.pointer());
}

void
LocalSymmSCMatrix::element_op(const Ref<SCElementOp2>& op,
                              SymmSCMatrix* m)
{
  LocalSymmSCMatrix *lm
      = require_dynamic_cast<LocalSymmSCMatrix*>(m,"LocalSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::errn() << indent << "LocalSymmSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_ltri(block.pointer(), lm->block.pointer());
}

void
LocalSymmSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  LocalSymmSCMatrix *lm
      = require_dynamic_cast<LocalSymmSCMatrix*>(m,"LocalSymSCMatrix::element_op");
  LocalSymmSCMatrix *ln
      = require_dynamic_cast<LocalSymmSCMatrix*>(n,"LocalSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::errn() << indent << "LocalSymmSCMatrix: bad element_op\n";
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
    os << endl << indent << title << endl;
  else
    os << endl;

  if (n()==0) {
    os << indent << "empty matrix\n";
    return;
  }

  for (ii=jj=0;;) {
    ii++; jj++;
    kk=width*jj;
    nn = (n() > kk) ? kk : n();

    // print column indices
    os << indent;
    for (i=ii; i <= nn; i++)
      os << scprintf("%*d",lwidth,i);
    os << endl;

    // print the rows
    for (i=ii-1; i < n() ; i++) {
      os << indent << scprintf("%5d",i+1);
      for (j=ii-1; j<nn && j<=i; j++)
        os << scprintf("%*.*f",lwidth,prec,rows[i][j]);
      os << endl;
    }
    os << endl;

    if (n() <= kk) {
      os.flush();
      return;
    }
    ii=kk;
  }
}

Ref<SCMatrixSubblockIter>
LocalSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  if (messagegrp()->n() > 1) {
      ExEnv::errn() << indent
           << "LocalSymmSCMatrix::local_blocks: not valid for local matrices"
           << endl;
      abort();
    }
  Ref<SCMatrixSubblockIter> iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

Ref<SCMatrixSubblockIter>
LocalSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << indent << "LocalSymmSCMatrix::all_blocks: "
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
