//
// localrect.cc
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

using namespace std;
using namespace sc;

extern "C" {
    int sing_(double *q, int *lq, int *iq, double *s, double *p,
              int *lp, int *ip, double *a, int *la, int *m, int *n,
              double *w);
};

/////////////////////////////////////////////////////////////////////////////
// LocalSCMatrix member functions

static ClassDesc LocalSCMatrix_cd(
  typeid(LocalSCMatrix),"LocalSCMatrix",1,"public SCMatrix",
  0, 0, 0);

static double **
init_rect_rows(double *data, int ni,int nj)
{
  double** r = new double*[ni];
  int i;
  for (i=0; i<ni; i++) r[i] = &data[i*nj];
  return r;
}

LocalSCMatrix::LocalSCMatrix(const RefSCDimension&a,const RefSCDimension&b,
                             LocalSCMatrixKit*kit):
  SCMatrix(a,b,kit),
  rows(0)
{
  resize(a->n(),b->n());
}

LocalSCMatrix::~LocalSCMatrix()
{
  if (rows) delete[] rows;
}

int
LocalSCMatrix::compute_offset(int i,int j) const
{
  if (i<0 || j<0 || i>=d1->n() || j>=d2->n()) {
      ExEnv::errn() << indent << "LocalSCMatrix: index out of bounds\n";
      abort();
    }
  return i*(d2->n()) + j;
}

void
LocalSCMatrix::resize(int nr, int nc)
{
  block = new SCMatrixRectBlock(0,nr,0,nc);
  if (rows) delete[] rows;
  rows = init_rect_rows(block->data,nr,nc);
}

double *
LocalSCMatrix::get_data()
{
  return block->data;
}

double **
LocalSCMatrix::get_rows()
{
  return rows;
}

double
LocalSCMatrix::get_element(int i,int j) const
{
  int off = compute_offset(i,j);
  return block->data[off];
}

void
LocalSCMatrix::set_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  block->data[off] = a;
}

void
LocalSCMatrix::accumulate_element(int i,int j,double a)
{
  int off = compute_offset(i,j);
  block->data[off] += a;
}

SCMatrix *
LocalSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::get_subblock: trying to get too big a subblock (" <<
          nsrow << "," << nscol << ") from (" <<
          nrow() << "," << ncol() << ")\n";
      abort();
    }
  
  RefSCDimension dnrow;
  if (nsrow==nrow()) dnrow = rowdim();
  else dnrow = new SCDimension(nsrow);

  RefSCDimension dncol;
  if (nscol==ncol()) dncol = coldim();
  else dncol = new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  LocalSCMatrix *lsb =
    require_dynamic_cast<LocalSCMatrix*>(sb, "LocalSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = rows[i+br][j+bc];
      
  return sb;
}

void
LocalSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                               int source_br, int source_bc)
{
  LocalSCMatrix *lsb =
    require_dynamic_cast<LocalSCMatrix*>(sb, "LocalSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::assign_subblock: trying to assign too big a " <<
          "subblock (" << nsrow << "," << nscol << " to (" <<
          nrow() << "," << ncol() << ")\n";
      abort();
    }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] = lsb->rows[source_br + i][source_bc + j];
}

void
LocalSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                   int source_br, int source_bc)
{
  LocalSCMatrix *lsb =
    require_dynamic_cast<LocalSCMatrix*>(sb, "LocalSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate_subblock: " <<
          "trying to accumulate too big a subblock (" <<
          nsrow << "," << nscol << " to (" <<
          nrow() << "," << ncol() << ")\n";
      abort();
    }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] += lsb->rows[source_br + i][source_bc + j]; 
}

SCVector *
LocalSCMatrix::get_row(int i)
{
  if (i >= nrow()) {
      ExEnv::errn() << indent << "LocalSCMatrix::get_row: trying to get invalid row " <<
          i << " max " << nrow() << endl;
      abort();
    }
  
  SCVector * v = kit()->vector(coldim());

  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::get_row");

  for (int j=0; j < ncol(); j++)
    lv->set_element(j,rows[i][j]);
      
  return v;
}

void
LocalSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= nrow()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::assign_row: trying to assign invalid row " <<
          i << " max " << nrow() << endl;
      abort();
    }
  
  if (v->n() != ncol()) {
      ExEnv::errn() << indent << "LocalSCMatrix::assign_row: vector is wrong size " <<
          " is " << v->n() << ", should be " << ncol() << endl;
      abort();
    }
  
  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::assign_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] = lv->get_element(j);
}

void
LocalSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= nrow()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_row: trying to assign invalid row " <<
          i << " max " << nrow() << endl;
      abort();
    }
  
  if (v->n() != ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_row: vector is wrong size " <<
          "is " << v->n() << ", should be " << ncol() << endl;
      abort();
    }
  
  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::accumulate_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] += lv->get_element(j);
}

SCVector *
LocalSCMatrix::get_column(int i)
{
  if (i >= ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::get_column: trying to get invalid column " <<
          i << " max " << ncol() << endl;
      abort();
    }
  
  SCVector * v = kit()->vector(rowdim());

  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::get_column");

  for (int j=0; j < nrow(); j++)
    lv->set_element(j,rows[j][i]);
      
  return v;
}

void
LocalSCMatrix::assign_column(SCVector *v, int i)
{
  if (i >= ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::assign_column: trying to assign invalid column " <<
          i << " max " << ncol() << endl;
      abort();
    }
  
  if (v->n() != nrow()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::assign_column: vector is wrong size " <<
          "is " << v->n() << ", should be " << nrow() << endl;
      abort();
    }
  
  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::assign_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] = lv->get_element(j);
}

void
LocalSCMatrix::accumulate_column(SCVector *v, int i)
{
  if (i >= ncol()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_column: trying to assign invalid column "
           << i << " max " << ncol() << endl;
      abort();
    }
  
  if (v->n() != nrow()) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_column: vector is wrong size " <<
          "is " << v->n() << ", should be " << nrow() << endl;
      abort();
    }
  
  LocalSCVector *lv =
    require_dynamic_cast<LocalSCVector*>(v, "LocalSCMatrix::accumulate_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] += lv->get_element(j);
}

void
LocalSCMatrix::assign_val(double a)
{
  int n = d1->n() * d2->n();
  double *data = block->data;
  for (int i=0; i<n; i++) data[i] = a;
}

void
LocalSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = require_dynamic_cast<LocalSCMatrix*>(a,name);
  LocalSCMatrix* lb = require_dynamic_cast<LocalSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate_product: bad dim" << endl;
      ExEnv::errn() << indent << "this row and col dimension:" << endl;
      rowdim()->print(ExEnv::errn());
      coldim()->print(ExEnv::errn());
      ExEnv::errn() << indent << "a row and col dimension:" << endl;
      a->rowdim()->print(ExEnv::errn());
      a->coldim()->print(ExEnv::errn());
      ExEnv::errn() << indent << "b row and col dimension:" << endl;
      b->rowdim()->print(ExEnv::errn());
      b->coldim()->print(ExEnv::errn());
      abort();
    }

  cmat_mxm(la->rows, 0,
           lb->rows, 0,
           rows, 0,
           nrow(), la->ncol(), this->ncol(),
           1);
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
LocalSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "LocalSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  LocalSCVector* la = require_dynamic_cast<LocalSCVector*>(a,name);
  LocalSCVector* lb = require_dynamic_cast<LocalSCVector*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b): " <<
          "dimensions don't match" << endl;
      abort();
    }

  int nr = a->n();
  int nc = b->n();
  int i, j;
  double* adat = la->block->data;
  double* bdat = lb->block->data;
  double** thisdat = rows;
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          thisdat[i][j] += adat[i] * bdat[j];
        }
    }
}

void
LocalSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = require_dynamic_cast<LocalSCMatrix*>(a,name);
  LocalSymmSCMatrix* lb = require_dynamic_cast<LocalSymmSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b): " <<
          "dimensions don't match" << endl;
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double **bd = lb->rows;
  int ni = a->rowdim().n();
  int njk = b->dim().n();
  int i, j, k;
  for (i=0; i<ni; i++) {
      for (j=0; j<njk; j++) {
          for (k=0; k<=j; k++) {
              cd[i][k] += ad[i][j]*bd[j][k];
            }
          for (; k<njk; k++) {
              cd[i][k] += ad[i][j]*bd[k][j];
            }
        }
    }
}

void
LocalSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "LocalSCMatrix::accumulate_product_rd";
  // make sure that the arguments are of the correct type
  LocalSCMatrix* la = require_dynamic_cast<LocalSCMatrix*>(a,name);
  LocalDiagSCMatrix* lb = require_dynamic_cast<LocalDiagSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent <<
          "LocalSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b): " <<
          "dimensions don't match" << endl;
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double *bd = lb->block->data;
  int ni = a->rowdim().n();
  int nj = b->dim().n();
  int i, j;
  for (i=0; i<ni; i++) {
      for (j=0; j<nj; j++) {
          cd[i][j] += ad[i][j]*bd[j];
        }
    }
}

void
LocalSCMatrix::accumulate(const SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const LocalSCMatrix* la
    = require_dynamic_cast<const LocalSCMatrix*>(a,"LocalSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate(SCMatrix*a): " <<
          "dimensions don't match" << endl;
      abort();
    }

  int nelem = this->ncol() * this->nrow();
  int i;
  for (i=0; i<nelem; i++) block->data[i] += la->block->data[i];
}

void
LocalSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const LocalSymmSCMatrix* la
    = require_dynamic_cast<const LocalSymmSCMatrix*>(a,"LocalSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate(SymmSCMatrix*a): " <<
          "dimensions don't match" << endl;
      abort();
    }

  int n = this->ncol();
  double *dat = la->block->data;
  int i, j;
  for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
          double tmp = *dat;
          block->data[i*n+j] += tmp;
          block->data[j*n+i] += tmp;
          dat++;
        }
      block->data[i*n+i] += *dat++;
    }
}

void
LocalSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const LocalDiagSCMatrix* la
    = require_dynamic_cast<const LocalDiagSCMatrix*>(a,"LocalSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate(DiagSCMatrix*a): " <<
          "dimensions don't match\n";
      abort();
    }

  int n = this->ncol();
  double *dat = la->block->data;
  int i;
  for (i=0; i<n; i++) {
      block->data[i*n+i] += *dat++;
    }
}

void
LocalSCMatrix::accumulate(const SCVector*a)
{
  // make sure that the arguments is of the correct type
  const LocalSCVector* la
    = require_dynamic_cast<const LocalSCVector*>(a,"LocalSCVector::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
      ExEnv::errn() << indent << "LocalSCMatrix::accumulate(SCVector*a): " <<
          "dimensions don't match" << endl;
      abort();
    }

  int n = this->ncol();
  double *dat = la->block->data;
  int i;
  for (i=0; i<n; i++) {
      block->data[i*n+i] += *dat++;
    }
}

void
LocalSCMatrix::transpose_this()
{
  cmat_transpose_matrix(rows,nrow(),ncol());
  delete[] rows;
  rows = new double*[ncol()];
  cmat_matrix_pointers(rows,block->data,ncol(),nrow());
  RefSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;

  int itmp = block->istart;
  block->istart = block->jstart;
  block->jstart = itmp;

  itmp = block->iend;
  block->iend = block->jend;
  block->jend = itmp;
}

double
LocalSCMatrix::invert_this()
{
  if (nrow() != ncol()) {
      ExEnv::errn() << indent << "LocalSCMatrix::invert_this: matrix is not square\n";
      abort();
    }
  return cmat_invert(rows,0,nrow());
}

double
LocalSCMatrix::determ_this()
{
  if (nrow() != ncol()) {
      ExEnv::errn() << indent << "LocalSCMatrix::determ_this: matrix is not square\n";
      abort();
    }
  return cmat_determ(rows,0,nrow());
}

double
LocalSCMatrix::trace()
{
  if (nrow() != ncol()) {
      ExEnv::errn() << indent << "LocalSCMatrix::trace: matrix is not square\n";
      abort();
    }
  double ret=0;
  int i;
  for (i=0; i < nrow(); i++)
    ret += rows[i][i];
  return ret;
}

void
LocalSCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  LocalSCMatrix* lU =
    require_dynamic_cast<LocalSCMatrix*>(U,"LocalSCMatrix::svd_this");
  LocalSCMatrix* lV =
    require_dynamic_cast<LocalSCMatrix*>(V,"LocalSCMatrix::svd_this");
  LocalDiagSCMatrix* lsigma =
    require_dynamic_cast<LocalDiagSCMatrix*>(sigma,"LocalSCMatrix::svd_this");

  RefSCDimension mdim = rowdim();
  RefSCDimension ndim = coldim();
  int m = mdim.n();
  int n = ndim.n();

  RefSCDimension pdim;
  if (m == n && m == sigma->dim().n())
    pdim = sigma->dim();
  else if (m<n)
    pdim = mdim;
  else
    pdim = ndim;

  int p = pdim.n();

  if (!mdim->equiv(lU->rowdim()) ||
      !mdim->equiv(lU->coldim()) ||
      !ndim->equiv(lV->rowdim()) ||
      !ndim->equiv(lV->coldim()) ||
      !pdim->equiv(sigma->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix: svd_this: dimension mismatch\n";
      abort();
    }

  // form a fortran style matrix for the SVD routines
  double *dA = new double[m*n];
  double *dU = new double[m*m];
  double *dV = new double[n*n];
  double *dsigma = new double[n];
  double *w = new double[(3*p-1>m)?(3*p-1):m];

  int i,j;
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          dA[i + j*m] = this->block->data[i*n + j];
        }
    }

  int three = 3;

  sing_(dU, &m, &three, dsigma, dV, &n, &three, dA, &m, &m, &n, w);

  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          lU->block->data[i*m + j] = dU[i + j*m];
        }
    }

  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          lV->block->data[i*n + j] = dV[i + j*n];
        }
    }

  for (i=0; i<p; i++) {
      lsigma->block->data[i] = dsigma[i];
    }

  delete[] dA;
  delete[] dU;
  delete[] dV;
  delete[] dsigma;
  delete[] w;
}

double
LocalSCMatrix::solve_this(SCVector*v)
{
  LocalSCVector* lv =
    require_dynamic_cast<LocalSCVector*>(v,"LocalSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::solve_this(SCVector*v): " <<
          "dimensions don't match" << endl;
      abort();
    }

  return cmat_solve_lin(rows,0,lv->block->data,nrow());
}

void
LocalSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  LocalSymmSCMatrix* lS =
    require_dynamic_cast<LocalSymmSCMatrix*>(S,"LocalSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::schmidt_orthog(): " <<
          "dimensions don't match\n";
      abort();
    }

  cmat_schmidt(rows,lS->block->data,nrow(),nc);
}

int
LocalSCMatrix::schmidt_orthog_tol(SymmSCMatrix *S, double tol, double *res)
{
  LocalSymmSCMatrix* lS =
    require_dynamic_cast<LocalSymmSCMatrix*>(S,"LocalSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      ExEnv::errn() << indent << "LocalSCMatrix::schmidt_orthog(): " <<
          "dimensions don't match\n";
      abort();
    }

  return cmat_schmidt_tol(rows,lS->block->data,nrow(),ncol(),tol,res);
}

void
LocalSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  op->process_spec_rect(block.pointer());
}

void
LocalSCMatrix::element_op(const Ref<SCElementOp2>& op,
                          SCMatrix* m)
{
  LocalSCMatrix *lm
      = require_dynamic_cast<LocalSCMatrix*>(m,"LocalSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
      ExEnv::errn() << indent << "LocalSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_rect(block.pointer(), lm->block.pointer());
}

void
LocalSCMatrix::element_op(const Ref<SCElementOp3>& op,
                          SCMatrix* m,SCMatrix* n)
{
  LocalSCMatrix *lm
      = require_dynamic_cast<LocalSCMatrix*>(m,"LocalSCMatrix::element_op");
  LocalSCMatrix *ln
      = require_dynamic_cast<LocalSCMatrix*>(n,"LocalSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
      ExEnv::errn() << indent << "LocalSCMatrix: bad element_op\n";
      abort();
    }
  op->process_spec_rect(block.pointer(),
                        lm->block.pointer(), ln->block.pointer());
}

// from Ed Seidl at the NIH
void
LocalSCMatrix::vprint(const char *title, ostream& os, int prec) const
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

  if (nrow()==0 || ncol()==0) {
    os << indent << "empty matrix\n";
    return;
  }

  for (ii=jj=0;;) {
    ii++; jj++;
    kk=width*jj;
    nn = (ncol()>kk) ? kk : ncol();

    // print column indices
    os << indent;
    for (i=ii; i <= nn; i++)
      os << scprintf("%*d",lwidth,i);
    os << endl;

    // print the rows
    for (i=0; i < nrow() ; i++) {
      os << indent << scprintf("%5d",i+1);
      for (j=ii-1; j < nn; j++)
        os << scprintf("%*.*f",lwidth,prec,rows[i][j]);
      os << endl;
    }
    os << endl;

    if (ncol() <= kk) {
      os.flush();
      return;
    }

    ii=kk;
  }
}

Ref<SCMatrixSubblockIter>
LocalSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  if (messagegrp()->n() > 1) {
      ExEnv::errn() << indent
           << "LocalSCMatrix::local_blocks: not valid for local matrices"
           << endl;
      abort();
    }
  Ref<SCMatrixSubblockIter> iter
      = new SCMatrixSimpleSubblockIter(access, block.pointer());
  return iter;
}

Ref<SCMatrixSubblockIter>
LocalSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << indent << "LocalSCMatrix::all_blocks: "
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
