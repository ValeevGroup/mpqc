//
// replrect.cc
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
#include <math/scmat/mops.h>

using namespace std;
using namespace sc;

extern "C" {
    int sing_(double *q, int *lq, int *iq, double *s, double *p,
              int *lp, int *ip, double *a, int *la, int *m, int *n,
              double *w);
};

/////////////////////////////////////////////////////////////////////////////
// ReplSCMatrix member functions

static ClassDesc ReplSCMatrix_cd(
  typeid(ReplSCMatrix),"ReplSCMatrix",1,"public SCMatrix",
  0, 0, 0);

static double **
init_rect_rows(double *data, int ni,int nj)
{
  double** r = new double*[ni];
  int i;
  size_t row_begin_index = 0;
  for (i=0; i<ni; i++, row_begin_index+=nj) r[i] = &data[row_begin_index];
  return r;
}

ReplSCMatrix::ReplSCMatrix(const RefSCDimension&a,const RefSCDimension&b,
                           ReplSCMatrixKit*k):
  SCMatrix(a,b,k)
{
  size_t nr = a->n();
  size_t nc = b->n();

  matrix = allocate<double>(nr*nc);

  rows = init_rect_rows(matrix,nr,nc);

  init_blocklist();
}

void
ReplSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i, j, index;
  size_t nc = d2->n();
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0, index=0; i<d1->blocks()->nblock(); i++) {
      for (j=0; j<d2->blocks()->nblock(); j++, index++) {
          if (index%nproc == me) continue;
          for (int ii=d1->blocks()->start(i);
               ii<d1->blocks()->fence(i); ii++) {
              for (int jj=d2->blocks()->start(j);
                   jj<d2->blocks()->fence(j); jj++) {
                  matrix[ii*nc + jj] = 0.0;
                }
            }
        }
    }
}

void
ReplSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d1->n()*d2->n());
}

void
ReplSCMatrix::init_blocklist()
{
  int i, j, index;
  int nc = d2->n();
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d1->blocks()->nblock(); i++) {
      for (j=0; j<d2->blocks()->nblock(); j++, index++) {
          if (index%nproc == me) {
              blocklist->insert(
                  new SCMatrixRectSubBlock(d1->blocks()->start(i),
                                           d1->blocks()->fence(i),
                                           nc,
                                           d2->blocks()->start(j),
                                           d2->blocks()->fence(j),
                                           matrix));
            }
        }
    }
}

ReplSCMatrix::~ReplSCMatrix()
{
  if (matrix) deallocate(matrix);
  if (rows) delete[] rows;
}

size_t
ReplSCMatrix::compute_offset(int i,int j) const
{
  if (i<0 || j<0 || i>=d1->n() || j>=d2->n()) {
      ExEnv::errn() << indent << "ReplSCMatrix: index out of bounds" << endl;
      abort();
    }
  return i*static_cast<size_t>(d2->n()) + j;
}

double
ReplSCMatrix::get_element(int i,int j) const
{
  size_t off = compute_offset(i,j);
  return matrix[off];
}

void
ReplSCMatrix::set_element(int i,int j,double a)
{
  size_t off = compute_offset(i,j);
  matrix[off] = a;
}

void
ReplSCMatrix::accumulate_element(int i,int j,double a)
{
  size_t off = compute_offset(i,j);
  matrix[off] += a;
}

SCMatrix *
ReplSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::get_subblock: trying to get too big a"
         << "subblock (" << nsrow << "," << nscol
         << ") from (" << nrow() << "," << ncol() << ")" << endl;
    abort();
  }
  
  RefSCDimension dnrow = (nsrow==nrow()) ? rowdim().pointer()
                                         : new SCDimension(nsrow);
  RefSCDimension dncol = (nscol==ncol()) ? coldim().pointer()
                                         : new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  ReplSCMatrix *lsb =
    require_dynamic_cast<ReplSCMatrix*>(sb, "ReplSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = rows[i+br][j+bc];
      
  return sb;
}

void
ReplSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                              int source_br, int source_bc)
{
  ReplSCMatrix *lsb = require_dynamic_cast<ReplSCMatrix*>(sb,
                                      "ReplSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    ExEnv::errn() << indent
         << "ReplSCMatrix::assign_subblock: trying to assign too big a"
         << "subblock (" << nsrow << "," << nscol
         << ") to (" << nrow() << "," << ncol() << ")" << endl;;
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] = lsb->rows[source_br + i][source_bc + j];
}

void
ReplSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec,
                                  int source_br, int source_bc)
{
  ReplSCMatrix *lsb = require_dynamic_cast<ReplSCMatrix*>(sb,
                                      "ReplSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > nrow() || nscol > ncol()) {
    ExEnv::errn() << indent
         << "ReplSCMatrix::accumulate_subblock: trying to accumulate to big a"
         << "subblock (" << nsrow << "," << nscol
         << ") to (" << nrow() << "," << ncol() << ")" << endl;
    abort();
  }
  
  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      rows[i+br][j+bc] += lsb->rows[source_br + i][source_bc + j];
}

SCVector *
ReplSCMatrix::get_row(int i)
{
  if (i >= nrow()) {
    ExEnv::errn() << indent << "ReplSCMatrix::get_row: trying to get invalid row "
         << i << " max " << nrow() << endl;
    abort();
  }
  
  SCVector * v = kit()->vector(coldim());

  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::get_row");

  for (int j=0; j < ncol(); j++)
    lv->set_element(j,rows[i][j]);
      
  return v;
}

void
ReplSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    ExEnv::errn() << indent << "ReplSCMatrix::assign_row: trying to assign invalid row "
         << i << " max " << nrow() << endl;
    abort();
  }
  
  if (v->n() != ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::assign_row: vector is wrong size, "
         << "is " << v->n() << ", should be " << ncol() << endl;
    abort();
  }
  
  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::assign_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] = lv->get_element(j);
}

void
ReplSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= nrow()) {
    ExEnv::errn() << indent
         << "ReplSCMatrix::accumulate_row: trying to accumulate invalid row "
         << i << " max " << nrow() << endl;
    abort();
  }
  
  if (v->n() != ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::accumulate_row: vector is wrong size, "
         << "is " << v->n() << ", should be " << ncol() << endl;
    abort();
  }
  
  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::accumulate_row");

  for (int j=0; j < ncol(); j++)
    rows[i][j] += lv->get_element(j);
}

SCVector *
ReplSCMatrix::get_column(int i)
{
  if (i >= ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::get_column: trying to get invalid column "
         << i << " max " << ncol() << endl;
    abort();
  }
  
  SCVector * v = kit()->vector(rowdim());

  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::get_column");

  for (int j=0; j < nrow(); j++)
    lv->set_element(j,rows[j][i]);
      
  return v;
}

void
ReplSCMatrix::assign_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    ExEnv::errn() << indent
         << "ReplSCMatrix::assign_column: trying to assign invalid column "
         << i << " max " << ncol() << endl;
    abort();
  }
  
  if (v->n() != nrow()) {
    ExEnv::errn() << indent << "ReplSCMatrix::assign_column: vector is wrong size, "
         << "is " << v->n() << ", should be " << nrow() << endl;
    abort();
  }
  
  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::assign_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] = lv->get_element(j);
}

void
ReplSCMatrix::accumulate_column(SCVector *v, int i)
{
  if (i >= ncol()) {
    ExEnv::errn() << indent
         << "ReplSCMatrix::accumulate_column: trying to accumulate invalid"
         << " column" << i << " max " << ncol() << endl;
    abort();
  }
  
  if (v->n() != nrow()) {
    ExEnv::errn() << indent << "ReplSCMatrix::accumulate_column: vector is wrong size, "
         << "is " << v->n() << ", should be " << nrow() << endl;
    abort();
  }
  
  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSCMatrix::accumulate_column");

  for (int j=0; j < nrow(); j++)
    rows[j][i] += lv->get_element(j);
}

void
ReplSCMatrix::assign_val(double a)
{
  size_t n = static_cast<size_t>(d1->n()) * d2->n();
  for (size_t i=0; i<n; i++) matrix[i] = a;
}

void
ReplSCMatrix::assign_p(const double*m)
{
  size_t n = static_cast<size_t>(d1->n()) * d2->n();
  memcpy(matrix, m, sizeof(double)*n);
}

void
ReplSCMatrix::assign_pp(const double**m)
{
  int n1 = d1->n();
  int n2 = d2->n();
  for (int i=0; i < n1; i++)
      for (int j=0; j < n2; j++)
          rows[i][j] = m[i][j];
}

void
ReplSCMatrix::convert_p(double*m) const
{
  size_t n = static_cast<size_t>(d1->n()) * d2->n();
  memcpy(m, matrix, sizeof(double)*n);
}

void
ReplSCMatrix::convert_pp(double**m) const
{
  int n1 = d1->n();
  int n2 = d2->n();
  for (int i=0; i < n1; i++)
      for (int j=0; j < n2; j++)
          m[i][j] = rows[i][j];
}

void
ReplSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product_rr";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = require_dynamic_cast<ReplSCMatrix*>(a,name);
  ReplSCMatrix* lb = require_dynamic_cast<ReplSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->coldim()) ||
      !la->coldim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent
           << "ReplSCMatrix::accumulate_product_rr(SCMatrix*a,SCMatrix*b): "
           << "dimensions don't match" << endl;
      abort();
    }

#if 0
  cmat_transpose_matrix(lb->rows, la->ncol(), this->ncol());
  double** btrans;
  btrans = new double*[this->ncol()];
  btrans[0] = lb->rows[0];
  cmat_matrix_pointers(btrans,btrans[0],this->ncol(),la->ncol());

  Ref<SCElementOp> op = new SCElementDot(la->rows, btrans, la->ncol());
  element_op(op);

  cmat_transpose_matrix(btrans,this->ncol(),la->ncol());
  delete[] btrans;
#else
  int i,j,k;
  int ii,jj;

  int nr = la->nrow();
  int nc = lb->ncol();
  int ncc = la->ncol();

  if (nr==0 || nc==0 || ncc==0)
    return;
  
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  int mod = nr%nproc;
  int nirow = nr/nproc + ((mod <= me) ? 0 : 1);
  int istart = (nr/nproc)*me + ((mod <= me) ? mod : me);
  int iend = istart+nirow;

  double ** ablock = cmat_new_square_matrix(D1);
  double ** bblock = cmat_new_square_matrix(D1);
  double ** cblock = cmat_new_square_matrix(D1);

  for (i=istart; i < iend; i += D1) {
    int ni = iend-i;
    if (ni > D1) ni = D1;
    
    for (j=0; j < nc; j += D1) {
      int nj = nc-j;
      if (nj > D1) nj = D1;

      memset(cblock[0], 0, sizeof(double)*D1*D1);

      for (k=0; k < ncc; k += D1) {
        int nk = ncc-k;
        if (nk > D1) nk = D1;

        copy_block(ablock, la->rows, i, ni, k, nk);
        copy_trans_block(bblock, lb->rows, j, nj, k, nk);
        mult_block(ablock, bblock, cblock, ni, nj, nk);
      }

      for (ii=0; ii < ni; ii++)
        for (jj=0; jj < nj; jj++)
          rows[i+ii][j+jj] += cblock[ii][jj];
    }
  }

  for (i=0; i < nproc; i++) {
    nirow = nr/nproc + ((mod <= i) ? 0 : 1);
    istart = (nr/nproc)*i + ((mod <= i) ? mod : i);
    if (!nirow)
      break;
    messagegrp()->bcast(rows[istart], nirow*nc, i);
  }

  cmat_delete_matrix(ablock);
  cmat_delete_matrix(bblock);
  cmat_delete_matrix(cblock);
        
#endif  
}

// does the outer product a x b.  this must have rowdim() == a->dim() and
// coldim() == b->dim()
void
ReplSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b)
{
  const char* name = "ReplSCMatrix::accumulate_outer_product";
  // make sure that the arguments are of the correct type
  ReplSCVector* la = require_dynamic_cast<ReplSCVector*>(a,name);
  ReplSCVector* lb = require_dynamic_cast<ReplSCVector*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent
           << "ReplSCMatrix::accumulate_outer_product(SCVector*a,SCVector*b): "
           << "dimensions don't match" << endl;
      abort();
    }

  int nr = a->n();
  int nc = b->n();
  int i, j;
  double* adat = la->vector;
  double* bdat = lb->vector;
  double** thisdat = rows;
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          thisdat[i][j] += adat[i] * bdat[j];
        }
    }
}

void
ReplSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product_rs";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = require_dynamic_cast<ReplSCMatrix*>(a,name);
  ReplSymmSCMatrix* lb = require_dynamic_cast<ReplSymmSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent
           << "ReplSCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b): "
           << "dimensions don't match" << endl;
      ExEnv::err0() << indent << "rowdim():" << endl;
      rowdim().print();
      ExEnv::err0() << indent << "coldim():" << endl;
      coldim().print();
      ExEnv::err0() << indent << "la->rowdim():" << endl;
      la->rowdim().print();
      ExEnv::err0() << indent << "la->coldim():" << endl;
      la->coldim().print();
      ExEnv::err0() << indent << "lb->dim():" << endl;
      lb->dim().print();
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
ReplSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b)
{
  const char* name = "ReplSCMatrix::accumulate_product_rd";
  // make sure that the arguments are of the correct type
  ReplSCMatrix* la = require_dynamic_cast<ReplSCMatrix*>(a,name);
  ReplDiagSCMatrix* lb = require_dynamic_cast<ReplDiagSCMatrix*>(b,name);

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(lb->dim()) ||
      !la->coldim()->equiv(lb->dim())) {
      ExEnv::errn() << indent
           << "ReplSCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b): "
           << "dimensions don't match" << endl;
      abort();
    }

  double **cd = rows;
  double **ad = la->rows;
  double *bd = lb->matrix;
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
ReplSCMatrix::accumulate(const SCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const ReplSCMatrix* la
    = require_dynamic_cast<const ReplSCMatrix*>(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->rowdim()) || !coldim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match" << endl;
      abort();
    }

  size_t nelem = static_cast<size_t>(this->ncol()) * this->nrow();
  size_t i;
  for (i=0; i<nelem; i++) matrix[i] += la->matrix[i];
}

void
ReplSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const ReplSymmSCMatrix* la
    = require_dynamic_cast<const ReplSymmSCMatrix*>(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::accumulate(SymmSCMatrix*a): "
           << "dimensions don't match" << endl;
      abort();
    }

  size_t n = this->ncol();
  double *dat = la->matrix;
  int i, j;
  for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
          double tmp = *dat;
          matrix[i*n+j] += tmp;
          matrix[j*n+i] += tmp;
          dat++;
        }
      matrix[i*n+i] += *dat++;
    }
}

void
ReplSCMatrix::accumulate(const DiagSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const ReplDiagSCMatrix* la
    = require_dynamic_cast<const ReplDiagSCMatrix*>(a,"ReplSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!rowdim()->equiv(la->dim()) || !coldim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::accumulate(DiagSCMatrix*a): "
           << "dimensions don't match" << endl;
      abort();
    }

  size_t n = this->ncol();
  double *dat = la->matrix;
  int i;
  for (i=0; i<n; i++) {
      matrix[i*n+i] += *dat++;
    }
}

void
ReplSCMatrix::accumulate(const SCVector*a)
{
  // make sure that the arguments is of the correct type
  const ReplSCVector* la
    = require_dynamic_cast<const ReplSCVector*>(a,"ReplSCVector::accumulate");

  // make sure that the dimensions match
  if (!((rowdim()->equiv(la->dim()) && coldim()->n() == 1)
        || (coldim()->equiv(la->dim()) && rowdim()->n() == 1))) {
      ExEnv::errn() << indent << "ReplSCMatrix::accumulate(SCVector*a): "
           << "dimensions don't match" << endl;
      abort();
    }

  size_t n = this->ncol();
  int i;
  double *dat = la->vector;
  for (i=0; i<n; i++) {
      matrix[i*n+i] += dat[i];
    }
}

void
ReplSCMatrix::transpose_this()
{
  cmat_transpose_matrix(rows,nrow(),ncol());
  delete[] rows;
  rows = new double*[ncol()];
  cmat_matrix_pointers(rows,matrix,ncol(),nrow());
  RefSCDimension tmp = d1;
  d1 = d2;
  d2 = tmp;
  init_blocklist();
}

double
ReplSCMatrix::invert_this()
{
  if (nrow() != ncol()) {
      ExEnv::errn() << indent << "ReplSCMatrix::invert_this: matrix is not square" << endl;
      abort();
    }
  return cmat_invert(rows,0,nrow());
}

double
ReplSCMatrix::determ_this()
{
  if (nrow() != ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::determ_this: matrix is not square" << endl;
    abort();
  }
  return cmat_determ(rows,0,nrow());
}

double
ReplSCMatrix::trace()
{
  if (nrow() != ncol()) {
    ExEnv::errn() << indent << "ReplSCMatrix::trace: matrix is not square" << endl;
    abort();
  }
  double ret=0;
  int i;
  for (i=0; i < nrow(); i++)
    ret += rows[i][i];
  return ret;
}

void
ReplSCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  ReplSCMatrix* lU =
    require_dynamic_cast<ReplSCMatrix*>(U,"ReplSCMatrix::svd_this");
  ReplSCMatrix* lV =
    require_dynamic_cast<ReplSCMatrix*>(V,"ReplSCMatrix::svd_this");
  ReplDiagSCMatrix* lsigma =
    require_dynamic_cast<ReplDiagSCMatrix*>(sigma,"ReplSCMatrix::svd_this");

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
      ExEnv::errn() << indent << "ReplSCMatrix: svd_this: dimension mismatch" << endl;
      abort();
    }

  // form a fortran style matrix for the SVD routines
  double *dA = allocate<double>(static_cast<size_t>(m)*n);
  double *dU = allocate<double>(static_cast<size_t>(m)*m);
  double *dV = allocate<double>(static_cast<size_t>(n)*n);
  double *dsigma = new double[n];
  double *w = new double[(3*p-1>m)?(3*p-1):m];

  int i,j;
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          dA[i + j*static_cast<size_t>(m)] = this->rows[i][j];
        }
    }

  int three = 3;

  sing_(dU, &m, &three, dsigma, dV, &n, &three, dA, &m, &m, &n, w);

  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          lU->rows[i][j] = dU[i + j*static_cast<size_t>(m)];
        }
    }

  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          lV->rows[i][j] = dV[i + j*static_cast<size_t>(n)];
        }
    }

  for (i=0; i<p; i++) {
      lsigma->matrix[i] = dsigma[i];
    }

  deallocate(dA);
  deallocate(dU);
  deallocate(dV);
  delete[] dsigma;
  delete[] w;
}

double
ReplSCMatrix::solve_this(SCVector*v)
{
  ReplSCVector* lv =
    require_dynamic_cast<ReplSCVector*>(v,"ReplSCMatrix::solve_this");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lv->dim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match" << endl;
      abort();
    }

  return cmat_solve_lin(rows,0,lv->vector,nrow());
}

void
ReplSCMatrix::schmidt_orthog(SymmSCMatrix *S, int nc)
{
  int i,j,ij;
  int m;

  ReplSymmSCMatrix* lS =
    require_dynamic_cast<ReplSymmSCMatrix*>(S,"ReplSCMatrix::schmidt_orthog");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::schmidt_orthog(): "
           << "dimensions don't match" << endl;
      abort();
    }

#if 0
  cmat_schmidt(rows,lS->matrix,nrow(),nc);
#else
  int me = messagegrp()->me();
  int nproc = messagegrp()->n();
  int nr = nrow();
  
  double vtmp;
  double *v = new double[nr];
  double *cm = new double[nr];

  double **sblock = cmat_new_square_matrix(D1);
  
  int mod = nc%nproc;
  int ncoli = nc/nproc + (mod <= me ? 0 : 1);
  int cstart = (nc/nproc)*me + (mod <= me ? mod : me);
  int cend = cstart+ncoli;

  // copy my columns to a rows of temp matrix
  double **cols = cmat_new_rect_matrix(ncoli, nr);
  for (i=cstart; i < cend; i++)
    for (j=0; j < nr; j++)
      cols[i-cstart][j] = rows[j][i];
    
  for (m=0; m < nc; m++) {
    // who has this column
    for (i=0; i < nproc; i++) {
      int ni = nc/nproc + (mod <= i ? 0 : 1);
      int csi = (nc/nproc)*i + (mod <= i ? mod : i);
      if (m >= csi && m < csi+ni) {
        if (i==me)
          memcpy(cm, cols[m-csi], sizeof(double)*nr);
        messagegrp()->bcast(cm, nr, i);
        break;
      }
    }
    
    memset(v, 0, sizeof(double)*nr);
    
    size_t ij = 0;
    for (i=0; i < nr; i += D1) {
      int ni = nr-i;
      if (ni > D1) ni = D1;
      
      for (j=0; j < nr; j += D1, ij++) {
        if (ij%nproc != me)
          continue;

        int nj = nr-j;
        if (nj > D1) nj = D1;
        
        copy_sym_block(sblock, lS->rows, i, ni, j, nj);
        
        for (int ii=0; ii < ni; ii++)
          for (int jj=0; jj < nj; jj++)
            v[i+ii] += cm[j+jj]*sblock[ii][jj];
      }
    }

    messagegrp()->sum(v, nr);

    for (i=0,vtmp=0.0; i < nr; i++)
      vtmp += v[i]*cm[i];

    if (!vtmp) {
      ExEnv::errn() << "cmat_schmidt: bogus" << endl;
      abort();
    }

    if (vtmp < 1.0e-15)
      vtmp = 1.0e-15;

    vtmp = 1.0/sqrt(vtmp);
    
    for (i=0; i < nr; i++) {
      v[i] *= vtmp;
      cm[i] *= vtmp;
    }

    if (m < nc-1) {
      for (i=m+1; i < nc; i++) {
        if (i < cstart)
          continue;
        if (i >= cend)
          break;
        
        double *ci = cols[i-cstart];
        
        for (j=0,vtmp=0.0; j < nr; j++)
          vtmp += v[j] * ci[j];
        for (j=0; j < nr; j++)
          ci[j] -= vtmp * cm[j];
      }
    }

    // if I own cm then put it back into cols
    if (m >= cstart && m < cend)
      memcpy(cols[m-cstart], cm, sizeof(double)*nr);
  }

  // now collect columns again
  for (i=0; i < nproc; i++) {
    int ni = nc/nproc + (mod <= i ? 0 : 1);
    int csi = (nc/nproc)*i + (mod <= i ? mod : i);
    for (j=0; j < ni; j++) {
      if (i==me) {
        messagegrp()->bcast(cols[j], nr, i);
        for (int k=0; k < nr; k++)
          rows[k][j+csi] = cols[j][k];
      }
      else {
        messagegrp()->bcast(cm, nr, i);
        for (int k=0; k < nr; k++)
          rows[k][j+csi] = cm[k];
      }
    }
  }

  cmat_delete_matrix(sblock);
  cmat_delete_matrix(cols);
  delete[] v;
  delete[] cm;
#endif
}

int
ReplSCMatrix::schmidt_orthog_tol(SymmSCMatrix *S, double tol, double *res)
{
  ReplSymmSCMatrix* lS =
    require_dynamic_cast<ReplSymmSCMatrix*>(S,"ReplSCMatrix::schmidt_orthog_tol");
  
  // make sure that the dimensions match
  if (!rowdim()->equiv(lS->dim())) {
      ExEnv::errn() << indent << "ReplSCMatrix::schmidt_orthog_tol(): " <<
          "dimensions don't match" << endl;
      abort();
    }

  int northog;

  if (messagegrp()->me() == 0) {
      northog = cmat_schmidt_tol(rows,lS->matrix,nrow(),ncol(),tol,res);
    }

  // make sure everybody ends up with the same data
  messagegrp()->bcast(northog);
  messagegrp()->bcast(*res);
  for (int i=0; i<nrow(); i++) {
      messagegrp()->bcast(rows[i],ncol());
    }

  return northog;
}

void
ReplSCMatrix::element_op(const Ref<SCElementOp>& op)
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
ReplSCMatrix::element_op(const Ref<SCElementOp2>& op,
                          SCMatrix* m)
{
  ReplSCMatrix *lm
      = require_dynamic_cast<ReplSCMatrix*>(m,"ReplSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim())) {
      ExEnv::errn() << indent << "ReplSCMatrix: bad element_op" << endl;
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
ReplSCMatrix::element_op(const Ref<SCElementOp3>& op,
                          SCMatrix* m,SCMatrix* n)
{
  ReplSCMatrix *lm
      = require_dynamic_cast<ReplSCMatrix*>(m,"ReplSCMatrix::element_op");
  ReplSCMatrix *ln
      = require_dynamic_cast<ReplSCMatrix*>(n,"ReplSCMatrix::element_op");

  if (!rowdim()->equiv(lm->rowdim()) || !coldim()->equiv(lm->coldim()) ||
      !rowdim()->equiv(ln->rowdim()) || !coldim()->equiv(ln->coldim())) {
      ExEnv::errn() << indent << "ReplSCMatrix: bad element_op" << endl;
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

// from Ed Seidl at the NIH
void
ReplSCMatrix::vprint(const char *title, ostream& os, int prec) const
{
  int ii,jj,kk,nn;
  int i,j;
  int lwidth,width;
  double max=this->maxabs();

  if (messagegrp()->me() != 0) return;

  max = (max==0.0) ? 1.0 : log10(max);
  if (max < 0.0) max=1.0;

  lwidth = prec + 5 + (int) max;
  width = 75/(lwidth+SCFormIO::getindent(os));

  if (title)
    os << endl << indent << title << endl;
  else
    os << endl;

  if (nrow()==0 || ncol()==0) {
    os << indent << "empty matrix" << endl;
    return;
  }

  for (ii=jj=0;;) {
    ii++; jj++;
    kk=width*jj;
    nn = (ncol() > kk) ? kk : ncol();

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
ReplSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          matrix, d1->n()*d2->n());
}

Ref<SCMatrixSubblockIter>
ReplSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << indent << "ReplSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  Ref<SCMatrixBlockList> allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCMatrixRectSubBlock(0, d1->n(), d1->n(),
                                                0, d2->n(), matrix));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          matrix, d1->n()*d2->n());
}

Ref<ReplSCMatrixKit>
ReplSCMatrix::skit()
{
  return dynamic_cast<ReplSCMatrixKit*>(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
