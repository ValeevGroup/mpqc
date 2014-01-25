//
// replsymm.cc
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
#include <math/scmat/disthql.h>
#include <math/scmat/offset.h>
#include <math/scmat/mops.h>
#include <math/scmat/util.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// ReplSymmSCMatrix member functions

static ClassDesc ReplSymmSCMatrix_cd(
  typeid(ReplSymmSCMatrix),"ReplSymmSCMatrix",1,"public SymmSCMatrix",
  0, 0, 0);

static double **
init_symm_rows(double *data, int n)
{
  double** r = new double*[n];
  for (int i=0; i<n; i++) r[i] = &data[(i*(i+1))/2];
  return r;
}

ReplSymmSCMatrix::ReplSymmSCMatrix(const RefSCDimension&a,ReplSCMatrixKit*k):
  SymmSCMatrix(a,k),
  rows(0)
{
  int n = d->n();

  matrix = allocate<double>(n*(n+1)>>1);
  rows = init_symm_rows(matrix,n);

  init_blocklist();
}

void
ReplSymmSCMatrix::before_elemop()
{
  // zero out the blocks not in my block list
  int i, j, index;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  for (i=0, index=0; i<d->blocks()->nblock(); i++) {
      for (j=0; j<=i; j++, index++) {
          if (index%nproc == me) continue;
          for (int ii=d->blocks()->start(i); ii<d->blocks()->fence(i); ii++) {
              for (int jj=d->blocks()->start(j);
                   jj < d->blocks()->fence(j) && jj <= ii;
                   jj++) {
                  matrix[(ii*(ii+1)>>1) + jj] = 0.0;
                }
            }
        }
    }
}

void
ReplSymmSCMatrix::after_elemop()
{
  messagegrp()->sum(matrix, d->n()*(d->n()+1)>>1);
}

void
ReplSymmSCMatrix::init_blocklist()
{
  int i, j, index;
  int nproc = messagegrp()->n();
  int me = messagegrp()->me();
  blocklist = new SCMatrixBlockList;
  for (i=0, index=0; i<d->blocks()->nblock(); i++) {
      for (j=0; j<=i; j++, index++) {
          if (index%nproc != me) continue;
          blocklist->insert(
              new SCMatrixLTriSubBlock(d->blocks()->start(i),
                                       d->blocks()->fence(i),
                                       d->blocks()->start(j),
                                       d->blocks()->fence(j),
                                       matrix));
        }
    }
}

ReplSymmSCMatrix::~ReplSymmSCMatrix()
{
  if (matrix) deallocate(matrix);
  if (rows) delete[] rows;
}

int
ReplSymmSCMatrix::compute_offset(int i,int j) const
{
  if (i<0 || j<0 || i>=d->n() || j>=d->n()) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix: index out of bounds\n";
      abort();
    }
  return ij_offset(i,j);
}

double
ReplSymmSCMatrix::get_element(int i,int j) const
{
  return matrix[compute_offset(i,j)];
}

void
ReplSymmSCMatrix::set_element(int i,int j,double a)
{
  matrix[compute_offset(i,j)] = a;
}

void
ReplSymmSCMatrix::accumulate_element(int i,int j,double a)
{
  matrix[compute_offset(i,j)] += a;
}

SCMatrix *
ReplSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::get_subblock: trying to get too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }

  RefSCDimension dnrow = (nsrow==n()) ? dim().pointer():new SCDimension(nsrow);
  RefSCDimension dncol = (nscol==n()) ? dim().pointer():new SCDimension(nscol);

  SCMatrix * sb = kit()->matrix(dnrow,dncol);
  sb->assign(0.0);

  ReplSCMatrix *lsb =
    require_dynamic_cast<ReplSCMatrix*>(sb, "ReplSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      lsb->rows[i][j] = get_element(i+br,j+bc);

  return sb;
}

SymmSCMatrix *
ReplSymmSCMatrix::get_subblock(int br, int er)
{
  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::get_subblock: trying to get too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") from (" << n() << "," << n() << ")\n";
    abort();
  }

  RefSCDimension dnrow = new SCDimension(nsrow);

  SymmSCMatrix * sb = kit()->symmmatrix(dnrow);
  sb->assign(0.0);

  ReplSymmSCMatrix *lsb =
    require_dynamic_cast<ReplSymmSCMatrix*>(sb, "ReplSymmSCMatrix::get_subblock");

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      lsb->rows[i][j] = get_element(i+br,j+br);

  return sb;
}

void
ReplSymmSCMatrix::assign_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  ReplSCMatrix *lsb = require_dynamic_cast<ReplSCMatrix*>(sb,
                                      "ReplSCMatrix::assign_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::assign_subblock: trying to assign too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") to (" << n() << "," << n() << ")\n";
    abort();
  }

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+bc,lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::assign_subblock(SymmSCMatrix*sb, int br, int er)
{
  ReplSymmSCMatrix *lsb = require_dynamic_cast<ReplSymmSCMatrix*>(sb,
                                      "ReplSymmSCMatrix::assign_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::assign_subblock: trying to assign too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") to (" << n() << "," << n() << ")\n";
    abort();
  }

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::accumulate_subblock(SCMatrix*sb, int br, int er, int bc, int ec)
{
  ReplSCMatrix *lsb = require_dynamic_cast<ReplSCMatrix*>(sb,
                                  "ReplSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;
  int nscol = ec-bc+1;

  if (nsrow > n() || nscol > n()) {
    ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate_subblock: "
         << "trying to accumulate too big a "
         << "subblock (" << nsrow << "," << nscol
         << ") to (" << n() << "," << n() << ")\n";
    abort();
  }

  for (int i=0; i < nsrow; i++)
    for (int j=0; j < nscol; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

void
ReplSymmSCMatrix::accumulate_subblock(SymmSCMatrix*sb, int br, int er)
{
  ReplSCMatrix *lsb = require_dynamic_cast<ReplSCMatrix*>(sb,
                                  "ReplSymmSCMatrix::accumulate_subblock");

  int nsrow = er-br+1;

  if (nsrow > n()) {
    ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate_subblock: trying to "
         << "accumulate too big a "
         << "subblock (" << nsrow << "," << nsrow
         << ") to (" << n() << "," << n() << ")\n";
    abort();
  }

  for (int i=0; i < nsrow; i++)
    for (int j=0; j <= i; j++)
      set_element(i+br,j+br,get_element(i+br,j+br)+lsb->rows[i][j]);
}

SCVector *
ReplSymmSCMatrix::get_row(int i)
{
  if (i >= n()) {
    ExEnv::errn() << indent << "ReplSymmSCMatrix::get_row: trying to get invalid row "
         << i << " max " << n() << endl;
    abort();
  }

  SCVector * v = kit()->vector(dim());

  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSymmSCMatrix::get_row");

  for (int j=0; j < n(); j++)
    lv->set_element(j,get_element(i,j));

  return v;
}

void
ReplSymmSCMatrix::assign_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::assign_row: trying to assign invalid row "
         << i << " max " << n() << endl;
    abort();
  }

  if (v->n() != n()) {
    ExEnv::errn() << indent << "ReplSymmSCMatrix::assign_row: vector is wrong size, "
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }

  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSymmSCMatrix::assign_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,lv->get_element(j));
}

void
ReplSymmSCMatrix::accumulate_row(SCVector *v, int i)
{
  if (i >= n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::accumulate_row: trying to assign invalide row "
         << i << " max " << n() << endl;
    abort();
  }

  if (v->n() != n()) {
    ExEnv::errn() << indent
         << "ReplSymmSCMatrix::accumulate_row: vector is wrong size, "
         << "is " << v->n() << ", should be " << n() << endl;
    abort();
  }

  ReplSCVector *lv =
    require_dynamic_cast<ReplSCVector*>(v, "ReplSymmSCMatrix::accumulate_row");

  for (int j=0; j < n(); j++)
    set_element(i,j,get_element(i,j)+lv->get_element(j));
}

void
ReplSymmSCMatrix::assign_val(double val)
{
  int n = (d->n()*(d->n()+1))/2;
  for (int i=0; i<n; i++) matrix[i] = val;
}

void
ReplSymmSCMatrix::assign_s(SymmSCMatrix*m)
{
  ReplSymmSCMatrix* lm = dynamic_cast<ReplSymmSCMatrix*>(m);
  if (lm && dim()->equiv(lm->dim())) {
      int d = i_offset(n());
      memcpy(matrix, lm->matrix, sizeof(double)*d);
    }
  else
      SymmSCMatrix::assign_s(m);
}

void
ReplSymmSCMatrix::assign_p(const double*m)
{
  int d = i_offset(n());
  memcpy(matrix, m, sizeof(double)*d);
}

void
ReplSymmSCMatrix::assign_pp(const double**m)
{
  for (int i=0; i < n(); i++)
      for (int j=0; j <= i; j++)
          rows[i][j] = m[i][j];
}

void
ReplSymmSCMatrix::convert_p(double*m) const
{
  int d = i_offset(n());
  memcpy(m, matrix, sizeof(double)*d);
}

void
ReplSymmSCMatrix::convert_pp(double**m) const
{
  for (int i=0; i < n(); i++)
      for (int j=0; j <= i; j++)
          m[i][j] = rows[i][j];
}

void
ReplSymmSCMatrix::scale(double s)
{
  int nelem = i_offset(n());
  for (int i=0; i < nelem; i++) matrix[i] *= s;
}

void
ReplSymmSCMatrix::accumulate(const SymmSCMatrix*a)
{
  // make sure that the arguments is of the correct type
  const ReplSymmSCMatrix* la
    = require_dynamic_cast<const ReplSymmSCMatrix*>(a,"ReplSymmSCMatrix::accumulate");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate(SCMatrix*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = (this->n() * (this->n() + 1))/2;
  for (int i=0; i<nelem; i++) matrix[i] += la->matrix[i];
}

double
ReplSymmSCMatrix::invert_this()
{
  if (messagegrp()->n() == 1)
      return cmat_invert(rows,1,n());
  else {
      RefDiagSCMatrix refa = kit()->diagmatrix(d);
      RefSCMatrix refb = kit()->matrix(d,d);
      diagonalize(refa.pointer(),refb.pointer());
      double determ = 1.0;
      for (int i=0; i<dim()->n(); i++) {
          double val = refa->get_element(i);
          determ *= val;
        }
      Ref<SCElementOp> op = new SCElementInvert(1.0e-12);
      refa->element_op(op.pointer());
      assign(0.0);
      accumulate_transform(refb.pointer(), refa.pointer());
      return determ;
    }
}

double
ReplSymmSCMatrix::determ_this()
{
  if (messagegrp()->n() == 1)
      return cmat_determ(rows,1,n());
  else {
      RefDiagSCMatrix refa = kit()->diagmatrix(d);
      RefSCMatrix refb = kit()->matrix(d,d);
      diagonalize(refa.pointer(),refb.pointer());
      double determ = 1.0;
      for (int i=0; i<dim()->n(); i++) {
          double val = refa->get_element(i);
          determ *= val;
        }
      return determ;
    }
}

double
ReplSymmSCMatrix::trace()
{
  double ret=0;
  for (int i=0; i < n(); i++) ret += rows[i][i];
  return ret;
}

double
ReplSymmSCMatrix::solve_this(SCVector*v)
{
  ReplSCVector* lv =
    require_dynamic_cast<ReplSCVector*>(v,"ReplSymmSCMatrix::solve_this");

  // make sure that the dimensions match
  if (!dim()->equiv(lv->dim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::solve_this(SCVector*v): "
           << "dimensions don't match\n";
      abort();
    }

  return cmat_solve_lin(rows,1,lv->vector,n());
}

void
ReplSymmSCMatrix::gen_invert_this(double condition_number_threshold)
{
  RefSCMatrix evecs = kit()->matrix(dim(),dim());
  RefDiagSCMatrix evals = kit()->diagmatrix(dim());

  const char *name = "ReplSymmSCMatrix::gen_invert_this";
  ReplDiagSCMatrix *levals = require_dynamic_cast<ReplDiagSCMatrix*>(evals,name);
  ReplSCMatrix *levecs = require_dynamic_cast<ReplSCMatrix*>(evecs,name);

  this->diagonalize(evals.pointer(), evecs.pointer());
  const double sigma_max = evals->maxabs();
  const double sigma_min_threshold = sigma_max / condition_number_threshold;
  for (int i=0; i < n(); i++) {
      if (fabs(levals->matrix[i]) > sigma_min_threshold)
          levals->matrix[i] = 1.0/levals->matrix[i];
      else
          levals->matrix[i] = 0;
    }

  assign(0.0);
  accumulate_transform(levecs, levals);
}

void
ReplSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b)
{
  int i;

  const char* name = "ReplSymmSCMatrix::diagonalize";
  // make sure that the arguments is of the correct type
  ReplDiagSCMatrix* la = require_dynamic_cast<ReplDiagSCMatrix*>(a,name);
  ReplSCMatrix* lb = require_dynamic_cast<ReplSCMatrix*>(b,name);

  if (!dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) || !dim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent
           << "ReplSymmSCMatrix::diagonalize(DiagSCMatrix*a,SCMatrix*b): "
           << "bad dims\n";
      abort();
    }

  // This sets up the index list of columns to be stored on this node
  int n = dim()->n();
  int me = messagegrp()->me();
  int nproc = messagegrp()->n();

  // if there is one processor or only a few vectors, do serial diagonalization
  if (nproc==1 || n <= 2) {
      double *eigvals = la->matrix;
      double **eigvecs = lb->rows;
      cmat_diag(rows, eigvals, eigvecs, n, 1, 1.0e-15);
    }
  else {
      Ref<MessageGrp> diagmsggrp;

      bool dosplit = n<nproc;

      if (dosplit) {
//           diagmsggrp = messagegrp()->split((me<n)?0:-1);
          std::set<int> diagset;
          for (int i=0; i<n; i++) diagset.insert(i);
          diagmsggrp = messagegrp()->subset(diagset);
      }
      else {
          diagmsggrp = messagegrp();
      }

      if (diagmsggrp) {
          int ndiagproc = diagmsggrp->n();
          int nvec = n/ndiagproc + (me<(n%ndiagproc)?1:0);
          int mvec = n/ndiagproc + ((n%ndiagproc) ? 1 : 0);

          int *ivec = new int[nvec];
          for (i=0; i<nvec; i++)
              ivec[i] = i*ndiagproc + me;

          double *eigvals = new double[n];
          double **eigvecs = cmat_new_rect_matrix(mvec, n);
          double **rect = cmat_new_rect_matrix(mvec, n);

          lb->assign(0.0);

          for (i=0; i < nvec; i++) {
              int c = ivec[i];
              int j;
              for (j=0; j <= c; j++)
                  rect[i][j] = rows[c][j];
              for (; j < n; j++)
                  rect[i][j] = rows[j][c];
          }

          dist_diagonalize(n, nvec, rect[0], eigvals, eigvecs[0], diagmsggrp);

          la->assign(eigvals);
          delete[] eigvals;

          int *tivec = new int [mvec];
          for (i=0; i < ndiagproc; i++) {
              int tnvec;

              if (i==me) {
                  diagmsggrp->bcast(nvec, me);
                  diagmsggrp->bcast(eigvecs[0], n*nvec, me);
                  diagmsggrp->bcast(ivec, nvec, me);
                  tnvec = nvec;
                  memcpy(tivec, ivec, sizeof(int)*nvec);
                  memcpy(rect[0], eigvecs[0], sizeof(double)*n*nvec);
              }
              else {
                  diagmsggrp->bcast(tnvec, i);
                  diagmsggrp->bcast(rect[0], n*tnvec, i);
                  diagmsggrp->bcast(tivec, tnvec, i);
              }

              for (int j=0; j < tnvec; j++) {
                  int c = tivec[j];
                  for (int k=0; k < n; k++)
                      lb->rows[k][c] = rect[j][k];
              }
          }

          delete[] ivec;
          delete[] tivec;
          cmat_delete_matrix(eigvecs);
          cmat_delete_matrix(rect);
      }
      if (dosplit) {
          // now the diagonalization is complete on the first ndiagproc nodes
          // broadcast the results to the remaining nodes.
          std::set<int> bcastset;
          bcastset.insert(0);
          for (int i=n; i<nproc; i++) bcastset.insert(i);
          Ref<MessageGrp> bcastgrp = messagegrp()->subset(bcastset);
//           Ref<MessageGrp> bcastgrp = messagegrp()->split((me==0||me>=n)?0:-1);
          if (bcastgrp) {
              bcastgrp->bcast(la->matrix,n);
              bcastgrp->bcast(lb->matrix,n*n);
          }
      }
    }
}

void
ReplSymmSCMatrix::eigensystem(SymmSCMatrix*s, DiagSCMatrix*a, SCMatrix*b)
{
  int i;

  const char* name = "ReplSymmSCMatrix::eigensystem";
  // make sure that the arguments is of the correct type
  ReplSymmSCMatrix* ls = require_dynamic_cast<ReplSymmSCMatrix*>(s,name);
  ReplDiagSCMatrix* la = require_dynamic_cast<ReplDiagSCMatrix*>(a,name);
  ReplSCMatrix* lb = require_dynamic_cast<ReplSCMatrix*>(b,name);

  if (!dim()->equiv(ls->dim()) ||
      !dim()->equiv(la->dim()) ||
      !dim()->equiv(lb->coldim()) ||
      !dim()->equiv(lb->rowdim())) {
      ExEnv::errn() << indent
           << "ReplSymmSCMatrix::eigensystem: bad dims\n";
      abort();
  }

  // This sets up the index list of columns to be stored on this node
  int n = dim()->n();
  int me = messagegrp()->me();
  int nproc = messagegrp()->n();

  // if there is one processor or only a few vectors, do serial diagonalization
  if (nproc==1) {
    double *eigvals = la->matrix;
    double **eigvecs = lb->rows;
    cmat_eigensystem(this->rows, ls->rows, eigvals, eigvecs, n, 1);
  }
  else {
    throw FeatureNotImplemented("ReplSymmSCMatrix::eigensystem with nproc > 1 is not implemented",
                                __FILE__, __LINE__, this->class_desc());
  }
}

// computes this += a * a.t
void
ReplSymmSCMatrix::accumulate_symmetric_product(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplSCMatrix* la
    = require_dynamic_cast<ReplSCMatrix*>(a,"ReplSymmSCMatrix::"
                                          "accumulate_symmetric_product");

  if (!dim()->equiv(la->rowdim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::"
           << "accumulate_symmetric_product(SCMatrix*a): bad dim\n";
      abort();
    }

  cmat_symmetric_mxm(rows,n(),la->rows,la->ncol(),1);
}

// computes this += a + a.t
void
ReplSymmSCMatrix::accumulate_symmetric_sum(SCMatrix*a)
{
  // make sure that the argument is of the correct type
  ReplSCMatrix* la
    = require_dynamic_cast<ReplSCMatrix*>(a,"ReplSymmSCMatrix::"
                                          "accumulate_symmetric_sum");

  if (!dim()->equiv(la->rowdim()) || !dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::"
           << "accumulate_symmetric_sum(SCMatrix*a): bad dim\n";
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
ReplSymmSCMatrix::accumulate_symmetric_outer_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = require_dynamic_cast<ReplSCVector*>(a,"ReplSymmSCMatrix::"
                                      "accumulate_symmetric_outer_product");

  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::"
           << "accumulate_symmetric_outer_product(SCMatrix*a): bad dim\n";
      abort();
    }

  int n = dim().n();
  double** tdat = this->rows;
  double* adat = la->vector;
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          tdat[i][j] += adat[i]*adat[j];
        }
    }
}

// this += a * b * transpose(a)
void
ReplSymmSCMatrix::accumulate_transform(SCMatrix*a,SymmSCMatrix*b,
                                       SCMatrix::Transform t)
{
  int i,j,k;
  int ii,jj;
  int nc, nr;

  // do the necessary castdowns
  ReplSCMatrix*la
    = require_dynamic_cast<ReplSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  ReplSymmSCMatrix*lb = require_dynamic_cast<ReplSymmSCMatrix*>(
      b,"%s::accumulate_transform", class_name());

  // check the dimensions
  if (t == SCMatrix::NormalTransform) {
    if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate_transform: bad dim\n";
      abort();
    }

    nc = lb->n();
    nr = la->nrow();
  } else {
    if (!dim()->equiv(la->coldim()) || !lb->dim()->equiv(la->rowdim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate_transform: bad dim\n";
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

  // if one processor then minimize the amount of memory used
  if (nproc == 1) {
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
  }

  // this version requires a full temp matrix be kept
  else {
    int me = messagegrp()->me();
    int mod = nr%nproc;
    int njrow = nr/nproc + ((mod <= me) ? 0 : 1);
    int jstart = (nr/nproc)*me + ((mod <= me) ? mod : me);
    int jend = jstart+njrow;

    double **temp = cmat_new_rect_matrix(nr,nc);
    memset(temp[0], 0, sizeof(double)*nr*nc);

    for (i=0; i < nc; i += D1) {
      int ni = nc-i;
      if (ni > D1) ni = D1;

      for (k=0; k < nc; k += D1) {
        int nk = nc-k;
        if (nk > D1) nk = D1;

        copy_sym_block(ablock, lb->rows, i, ni, k, nk);

        for (j=jstart; j < jend; j += D1) {
          int nj = jend-j;
          if (nj > D1) nj = D1;

          if (t == SCMatrix::NormalTransform)
            copy_block(bblock, la->rows, j, nj, k, nk);
          else
            copy_trans_block(bblock, la->rows, j, nj, k, nk);

          memset(cblock[0], 0, sizeof(double)*D1*D1);
          mult_block(ablock, bblock, cblock, ni, nj, nk);

          for (jj=0; jj < nj; jj++)
            for (ii=0; ii < ni; ii++)
              temp[j+jj][i+ii] += cblock[ii][jj];
        }
      }
    }

    for (i=0; i < nproc; i++) {
      njrow = nr/nproc + ((mod <= i) ? 0 : 1);
      jstart = (nr/nproc)*i + ((mod <= i) ? mod : i);
      if (!njrow)
        break;
      messagegrp()->bcast(temp[jstart], njrow*nc, i);
    }

    int ind=0;
    for (i=0; i < nr; i += D1) {
      int ni = nr-i;
      if (ni > D1) ni = D1;

      for (j=0; j <= i; j += D1, ind++) {
        if (ind%nproc != me)
          continue;

        int nj = nr-j;
        if (nj > D1) nj = D1;

        memset(cblock[0], 0, sizeof(double)*D1*D1);

        for (k=0; k < nc; k += D1) {
          int nk = nc-k;
          if (nk > D1) nk = D1;

          if (t == SCMatrix::NormalTransform)
            copy_block(ablock, la->rows, i, ni, k, nk);
          else
            copy_trans_block(ablock, la->rows, i, ni, k, nk);
          copy_block(bblock, temp, j, nj, k, nk);
          mult_block(ablock, bblock, cblock, ni, nj, nk);
        }

        if (i==j) {
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

    ind=0;
    for (i=0; i < nr; i += D1) {
      int ni = nr-i;
      if (ni > D1) ni = D1;

      for (j=0; j <= i; j += D1, ind++) {
        int nj = nr-j;
        if (nj > D1) nj = D1;

        int proc = ind%nproc;

        if (proc==me)
          copy_sym_block(ablock, rows, i, ni, j, nj);

        messagegrp()->bcast(ablock[0], D1*D1, proc);

        if (i==j) {
          for (ii=0; ii < ni; ii++)
            for (jj=0; jj <= ii; jj++)
              rows[i+ii][j+jj] = ablock[ii][jj];
        } else {
          for (ii=0; ii < ni; ii++)
            for (jj=0; jj < nj; jj++)
              rows[i+ii][j+jj] = ablock[ii][jj];
        }
      }
    }

    cmat_delete_matrix(temp);
  }

  cmat_delete_matrix(ablock);
  cmat_delete_matrix(bblock);
  cmat_delete_matrix(cblock);
}

// this += a * b * transpose(a)
void
ReplSymmSCMatrix::accumulate_transform(SCMatrix*a,DiagSCMatrix*b,
                                       SCMatrix::Transform t)
{
  // do the necessary castdowns
  ReplSCMatrix*la
    = require_dynamic_cast<ReplSCMatrix*>(a,"%s::accumulate_transform",
                                      class_name());
  ReplDiagSCMatrix*lb
    = require_dynamic_cast<ReplDiagSCMatrix*>(b,"%s::accumulate_transform",
                                          class_name());

  // check the dimensions
  if (!dim()->equiv(la->rowdim()) || !lb->dim()->equiv(la->coldim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::accumulate_transform: bad dim\n";
      abort();
    }

  cmat_transform_diagonal_matrix(rows,n(),lb->matrix,lb->n(),la->rows,1);
}

void
ReplSymmSCMatrix::accumulate_transform(SymmSCMatrix*a,SymmSCMatrix*b)
{
  SymmSCMatrix::accumulate_transform(a,b);
}

double
ReplSymmSCMatrix::scalar_product(SCVector*a)
{
  // make sure that the argument is of the correct type
  ReplSCVector* la
    = require_dynamic_cast<ReplSCVector*>(a,"ReplSCVector::scalar_product");

  // make sure that the dimensions match
  if (!dim()->equiv(la->dim())) {
      ExEnv::errn() << indent << "ReplSCVector::scale_product(SCVector*a): "
           << "dimensions don't match\n";
      abort();
    }

  int nelem = n();
  double* adat = la->vector;
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
ReplSymmSCMatrix::element_op(const Ref<SCElementOp>& op)
{
  if (op->has_side_effects()) before_elemop();
  scmat_perform_op_on_blocks(op, blocklist);
  if (op->has_side_effects()) after_elemop();
  if (op->has_collect()) op->collect(messagegrp());
}

void
ReplSymmSCMatrix::element_op(const Ref<SCElementOp2>& op,
                              SymmSCMatrix* m)
{
  ReplSymmSCMatrix *lm
      = require_dynamic_cast<ReplSymmSCMatrix*>(m,"ReplSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix: bad element_op\n";
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
ReplSymmSCMatrix::element_op(const Ref<SCElementOp3>& op,
                              SymmSCMatrix* m,SymmSCMatrix* n)
{
  ReplSymmSCMatrix *lm
      = require_dynamic_cast<ReplSymmSCMatrix*>(m,"ReplSymSCMatrix::element_op");
  ReplSymmSCMatrix *ln
      = require_dynamic_cast<ReplSymmSCMatrix*>(n,"ReplSymSCMatrix::element_op");

  if (!dim()->equiv(lm->dim()) || !dim()->equiv(ln->dim())) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix: bad element_op\n";
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
ReplSymmSCMatrix::vprint(const char *title, ostream& os, int prec) const
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
ReplSymmSCMatrix::local_blocks(SCMatrixSubblockIter::Access access)
{
  return new ReplSCMatrixListSubblockIter(access, blocklist,
                                          messagegrp(),
                                          matrix, (d->n()*(d->n()+1))/2);
}

Ref<SCMatrixSubblockIter>
ReplSymmSCMatrix::all_blocks(SCMatrixSubblockIter::Access access)
{
  if (access == SCMatrixSubblockIter::Write) {
      ExEnv::errn() << indent << "ReplSymmSCMatrix::all_blocks: "
           << "Write access permitted for local blocks only"
           << endl;
      abort();
    }
  Ref<SCMatrixBlockList> allblocklist = new SCMatrixBlockList();
  allblocklist->insert(new SCMatrixLTriSubBlock(0, d->n(),
                                                0, d->n(), matrix));
  return new ReplSCMatrixListSubblockIter(access, allblocklist,
                                          messagegrp(),
                                          matrix, (d->n()*(d->n()+1))/2);
}

Ref<ReplSCMatrixKit>
ReplSymmSCMatrix::skit()
{
  return dynamic_cast<ReplSCMatrixKit*>(kit().pointer());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
