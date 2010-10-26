//
// matrix.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/state/stateio.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>

using namespace std;

namespace sc {

/////////////////////////////////////////////////////////////////////////////
// SCDimension reference member functions

RefSCDimension::RefSCDimension() {}

RefSCDimension::RefSCDimension (const RefSCDimension & o):
  Ref<SCDimension> (o) {}

RefSCDimension::RefSCDimension (SCDimension * o): Ref<SCDimension> (o) {}

RefSCDimension::~RefSCDimension () {}

RefSCDimension&
RefSCDimension::operator=(SCDimension* cr)
{
  Ref<SCDimension>::operator=(cr);
  return *this;
}

RefSCDimension&
RefSCDimension::operator<<(const RefBase & c)
{
  Ref<SCDimension>::operator<<(c);
  return *this;
}

RefSCDimension&
RefSCDimension::operator<<(RefCount*a)
{
  Ref<SCDimension>::operator<<(a);
  return *this;
}

RefSCDimension&
RefSCDimension::operator=(const RefSCDimension & c)
{
  Ref<SCDimension>::operator=(c);
  return *this;
}

int
RefSCDimension::n() const
{
  int result;
  if (null()) result = 0;
  else result = pointer()->n();
  return result;
}

RefSCDimension::operator int() const
{
  if (null()) return 0;
  return pointer()->n();
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrix reference member functions
RefSCMatrix::RefSCMatrix() {}

RefSCMatrix::RefSCMatrix (const RefSCMatrix & o): Ref<SCMatrix> (o) {}

RefSCMatrix::RefSCMatrix (SCMatrix * o): Ref<SCMatrix> (o) {}

RefSCMatrix::~RefSCMatrix () {}

RefSCMatrix&
RefSCMatrix::operator=(SCMatrix* cr)
{
  Ref<SCMatrix>::operator=(cr);
  return *this;
}

RefSCMatrix&
RefSCMatrix::operator=(const RefSCMatrix & c)
{
  Ref<SCMatrix>::operator=(c);
  return *this;
}

RefSCMatrix::RefSCMatrix(const RefSCDimension&a,const RefSCDimension&b,
                         const Ref<SCMatrixKit>&k)
{
  assign_pointer(k->matrix(a,b));
}

void
RefSCMatrix::set_element(int i, int j, double a) const
{
  require_nonnull();
  pointer()->set_element(i,j,a);
}

void
RefSCMatrix::accumulate_element(int i, int j, double a) const
{
  require_nonnull();
  pointer()->accumulate_element(i,j,a);
}

double
RefSCMatrix::get_element(int i, int j) const
{
  require_nonnull();
  return pointer()->get_element(i,j);
}

RefSCVector
RefSCMatrix::operator*(const RefSCVector&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector r = kit()->vector(rowdim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(const RefSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(rowdim(),a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(rowdim(),a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(rowdim(),a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator+(const RefSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim(),kit());

  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSCMatrix
RefSCMatrix::operator-(const RefSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim(),kit());

  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefSCMatrix
RefSCMatrix::t() const
{
  require_nonnull();

  RefSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->transpose_this();
  return ret;
}

RefSCMatrix
RefSCMatrix::i() const
{
  require_nonnull();

  RefSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

RefSCMatrix
RefSCMatrix::gi(double condition_number_threshold) const
{
  require_nonnull();

  RefSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->gen_invert_this(condition_number_threshold);
  return ret;
}

int
RefSCMatrix::nrow() const
{
  if (null()) return 0;
  else return pointer()->nrow();
}

int
RefSCMatrix::ncol() const
{
  if (null()) return 0;
  else return pointer()->ncol();
}

RefSCDimension
RefSCMatrix::rowdim() const
{
  if (null()) return 0;
  else return pointer()->rowdim();
}

RefSCDimension
RefSCMatrix::coldim() const
{
  if (null()) return 0;
  else return pointer()->coldim();
}

Ref<SCMatrixKit>
RefSCMatrix::kit() const
{
  if (null()) return 0;
  else return pointer()->kit();
}

SCMatrixdouble
RefSCMatrix::operator()(int i,int j)  const
{
  return SCMatrixdouble(pointer(),i,j);
}

RefSCMatrix
RefSCMatrix::clone() const
{
  RefSCMatrix r = kit()->matrix(rowdim(),coldim());
  return r;
}

RefSCMatrix
RefSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  require_nonnull();

  RefSCMatrix ret = pointer()->get_subblock(br,er,bc,ec);
  return ret;
}

void
RefSCMatrix::assign_subblock(const RefSCMatrix& sb,
                             int br, int er, int bc, int ec, int sbr, int sbc)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->assign_subblock(sb.pointer(),br,er,bc,ec,sbr,sbc);
}

void
RefSCMatrix::accumulate_subblock(const RefSCMatrix& sb,
                                 int br, int er, int bc, int ec,
                                 int sbr, int sbc)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->accumulate_subblock(sb.pointer(),br,er,bc,ec,sbr,sbc);
}

RefSCVector
RefSCMatrix::get_row(int i) const
{
  require_nonnull();

  RefSCVector ret = pointer()->get_row(i);
  return ret;
}

RefSCVector
RefSCMatrix::get_column(int i) const
{
  require_nonnull();

  RefSCVector ret = pointer()->get_column(i);
  return ret;
}

void
RefSCMatrix::assign_row(const RefSCVector& v, int i) const
{
  require_nonnull();
  v.require_nonnull();
  pointer()->assign_row(v.pointer(),i);
}

void
RefSCMatrix::assign_column(const RefSCVector& v, int i) const
{
  require_nonnull();
  v.require_nonnull();
  pointer()->assign_column(v.pointer(),i);
}

void
RefSCMatrix::accumulate_row(const RefSCVector& v, int i) const
{
  require_nonnull();
  v.require_nonnull();
  pointer()->accumulate_row(v.pointer(),i);
}

void
RefSCMatrix::accumulate_column(const RefSCVector& v, int i) const
{
  require_nonnull();
  v.require_nonnull();
  pointer()->accumulate_column(v.pointer(),i);
}

void
RefSCMatrix::accumulate_product(const RefSCMatrix&a,const RefSCMatrix&b) const
{
  require_nonnull();
  pointer()->accumulate_product(a.pointer(),b.pointer());
}

RefSCMatrix
RefSCMatrix::copy() const
{
  if (null()) return 0;
  RefSCMatrix v = kit()->matrix(rowdim(),coldim());
  v.assign(*this);
  return v;
}

void
RefSCMatrix::randomize() const
{
  require_nonnull();
  pointer()->randomize();
}

void
RefSCMatrix::assign(const RefSCMatrix&a) const
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSCMatrix::assign(const double*v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCMatrix::assign(const double**v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCMatrix::convert(double*v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCMatrix::convert(double**v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCMatrix::scale(double a) const
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSCMatrix::assign(double a) const
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSCMatrix::accumulate(const RefSCMatrix&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCMatrix::accumulate(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCMatrix::accumulate(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCMatrix::element_op(const Ref<SCElementOp>&op) const
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSCMatrix::element_op(const Ref<SCElementOp2>& op,
                        const RefSCMatrix& m) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer());
}

void
RefSCMatrix::element_op(const Ref<SCElementOp3>& op,
                        const RefSCMatrix& m,
                        const RefSCMatrix& n) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer(),n.pointer());
}

void
RefSCMatrix::print(ostream& out) const
{
  print(0,out);
}

void
RefSCMatrix::print(const char*title,ostream&out, int precision) const
{
  if (nonnull()) {
      pointer()->print(title,out,precision);
    }
  else {
      if (title) out << endl << title << endl;
      out << "null matrix" << endl;
    }
}

RefSCMatrix
RefSCMatrix::operator *(double a) const
{
  RefSCMatrix r(copy());
  r.scale(a);
  return r;
}

void
RefSCMatrix::svd(const RefSCMatrix &U,
                 const RefDiagSCMatrix &sigma,
                 const RefSCMatrix &V)
{
  require_nonnull();
  RefSCMatrix c = clone();
  c->assign(pointer());
  c->svd_this(U.pointer(), sigma.pointer(), V.pointer());
}

double
RefSCMatrix::solve_lin(const RefSCVector& v) const
{
  require_nonnull();
  RefSCMatrix c = clone();
  c->assign(pointer());
  return c->solve_this(v.pointer());
}

double
RefSCMatrix::determ() const
{
  require_nonnull();
  RefSCMatrix c = clone();
  c->assign(pointer());
  return c->determ_this();
}

double
RefSCMatrix::trace() const
{
  require_nonnull();
  return pointer()->trace();
}

RefSCMatrix
operator *(double a, const RefSCMatrix& v)
{
  return v*a;
}

void
RefSCMatrix::accumulate_outer_product(const RefSCVector& v1,
                                      const RefSCVector&v2) const
{
  require_nonnull();
  pointer()->accumulate_outer_product(v1.pointer(),v2.pointer());
}

void
RefSCMatrix::save(StateOut&s)
{
  if (null()) s.put(0);
  else {
      s.put(1);
      pointer()->save(s);
    }
}

void
RefSCMatrix::restore(StateIn&s)
{
  int have_matrix;
  s.get(have_matrix);
  if (have_matrix && nonnull()) {
      pointer()->restore(s);
    }
  else if (have_matrix) {
      ExEnv::errn() << "RefSCMatrix::restore: matrix not properly initialized" << endl;
      abort();
    }
  else {
      clear();
    }
}

int
RefSCMatrix::nblock() const
{
  BlockedSCMatrix *b = dynamic_cast<BlockedSCMatrix*>(pointer());
  if (b) return b->nblocks();
  return 1;
}

RefSCMatrix
RefSCMatrix::block(int i) const
{
  BlockedSCMatrix *b = dynamic_cast<BlockedSCMatrix*>(pointer());
  if (b) return b->block(i);
  return *this;
}

///////////////////////////////////////////////////////////////////
// RefSymmSCMatrix members

RefSymmSCMatrix::RefSymmSCMatrix()
{
}

RefSymmSCMatrix::RefSymmSCMatrix (const RefSymmSCMatrix & o):
  Ref<SymmSCMatrix> (o)
{
}

RefSymmSCMatrix::RefSymmSCMatrix (SymmSCMatrix * o):
  Ref<SymmSCMatrix> (o)
{
}

RefSymmSCMatrix::~RefSymmSCMatrix ()
{
}

RefSymmSCMatrix&
RefSymmSCMatrix::operator=(SymmSCMatrix* cr)
{
  Ref<SymmSCMatrix>::operator=(cr);
  return *this;
}

RefSymmSCMatrix&
RefSymmSCMatrix::operator=(const RefSymmSCMatrix & c)
{
  Ref<SymmSCMatrix>::operator=(c);
  return *this;
}

RefSymmSCMatrix::RefSymmSCMatrix(const RefSCDimension&a,
                                 const Ref<SCMatrixKit>&k)
{
  assign_pointer(k->symmmatrix(a));
}

void
RefSymmSCMatrix::set_element(int i, int j, double a) const
{
  require_nonnull();
  pointer()->set_element(i,j,a);
}

void
RefSymmSCMatrix::accumulate_element(int i, int j, double a) const
{
  require_nonnull();
  pointer()->accumulate_element(i,j,a);
}

double
RefSymmSCMatrix::get_element(int i, int j) const
{
  require_nonnull();
  return pointer()->get_element(i,j);
}

RefSCMatrix
RefSymmSCMatrix::get_subblock(int br, int er, int bc, int ec)
{
  require_nonnull();

  RefSCMatrix ret = pointer()->get_subblock(br,er,bc,ec);
  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::get_subblock(int br, int er)
{
  require_nonnull();

  RefSymmSCMatrix ret = pointer()->get_subblock(br,er);
  return ret;
}

void
RefSymmSCMatrix::assign_subblock(const RefSCMatrix& sb,
                                 int br, int er, int bc, int ec)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->assign_subblock(sb.pointer(),br,er,bc,ec);
}

void
RefSymmSCMatrix::assign_subblock(const RefSymmSCMatrix& sb, int br, int er)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->assign_subblock(sb.pointer(),br,er);
}

void
RefSymmSCMatrix::accumulate_subblock(const RefSCMatrix& sb,
                                     int br, int er, int bc, int ec)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->accumulate_subblock(sb.pointer(),br,er,bc,ec);
}

void
RefSymmSCMatrix::accumulate_subblock(const RefSymmSCMatrix& sb, int br, int er)
{
  require_nonnull();
  sb.require_nonnull();
  pointer()->accumulate_subblock(sb.pointer(),br,er);
}

RefSCVector
RefSymmSCMatrix::get_row(int i)
{
  require_nonnull();

  RefSCVector ret = pointer()->get_row(i);
  return ret;
}

void
RefSymmSCMatrix::assign_row(const RefSCVector& v, int i)
{
  require_nonnull();
  v.require_nonnull();
  pointer()->assign_row(v.pointer(),i);
}

void
RefSymmSCMatrix::accumulate_row(const RefSCVector& v, int i)
{
  require_nonnull();
  v.require_nonnull();
  pointer()->accumulate_row(v.pointer(),i);
}

void
RefSymmSCMatrix::accumulate_symmetric_product(const RefSCMatrix& a) const
{
  require_nonnull();
  pointer()->accumulate_symmetric_product(a.pointer());
}

void
RefSymmSCMatrix::accumulate_symmetric_sum(const RefSCMatrix& a) const
{
  require_nonnull();
  pointer()->accumulate_symmetric_sum(a.pointer());
}

void
RefSymmSCMatrix::accumulate_transform(const RefSCMatrix& a,
                                      const RefSymmSCMatrix&b,
                                      SCMatrix::Transform t) const
{
  require_nonnull();
  pointer()->accumulate_transform(a.pointer(),b.pointer(),t);
}

void
RefSymmSCMatrix::accumulate_transform(const RefSCMatrix& a,
                                      const RefDiagSCMatrix&b,
                                      SCMatrix::Transform t) const
{
  require_nonnull();
  pointer()->accumulate_transform(a.pointer(),b.pointer(),t);
}

void
RefSymmSCMatrix::accumulate_transform(const RefSymmSCMatrix& a,
                                      const RefSymmSCMatrix&b) const
{
  require_nonnull();
  pointer()->accumulate_transform(a.pointer(),b.pointer());
}

RefSymmSCMatrix
RefSymmSCMatrix::operator+(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim(),kit());

  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::operator-(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim(),kit());

  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::i() const
{
  require_nonnull();

  RefSymmSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::gi(double condition_number_threshold) const
{
  require_nonnull();

  RefSymmSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->gen_invert_this(condition_number_threshold);
  return ret;
}

int
RefSymmSCMatrix::n() const
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefSymmSCMatrix::dim() const
{
  if (null()) return 0;
  else return pointer()->dim();
}

Ref<SCMatrixKit>
RefSymmSCMatrix::kit() const
{
  if (null()) return 0;
  else return pointer()->kit();
}

SymmSCMatrixdouble
RefSymmSCMatrix::operator()(int i,int j) const
{
  return SymmSCMatrixdouble(pointer(),i,j);
}

RefSymmSCMatrix
RefSymmSCMatrix::clone() const
{
  RefSymmSCMatrix r = kit()->symmmatrix(dim());
  return r;
}

RefDiagSCMatrix
RefSymmSCMatrix::eigvals() const
{
  if (null()) return 0;
  RefDiagSCMatrix vals = kit()->diagmatrix(dim());
  RefSCMatrix vecs = kit()->matrix(dim(),dim());
  diagonalize(vals,vecs);
  return vals;
}

RefSCMatrix
RefSymmSCMatrix::eigvecs() const
{
  if (null()) return 0;
  RefDiagSCMatrix vals = kit()->diagmatrix(dim());
  RefSCMatrix vecs = kit()->matrix(dim(),dim());
  diagonalize(vals,vecs);
  return vecs;
}

double
RefSymmSCMatrix::solve_lin(const RefSCVector& v) const
{
  require_nonnull();

  RefSymmSCMatrix ret = clone();
  ret->assign(pointer());
  return ret->solve_this(v.pointer());
}

double
RefSymmSCMatrix::determ() const
{
  require_nonnull();

  RefSymmSCMatrix ret = clone();
  ret->assign(pointer());
  return ret->determ_this();
}

double
RefSymmSCMatrix::trace() const
{
  require_nonnull();
  return pointer()->trace();
}

void
RefSymmSCMatrix::diagonalize(const RefDiagSCMatrix& vals,
                             const RefSCMatrix& vecs) const
{
  require_nonnull();
  pointer()->diagonalize(vals.pointer(),vecs.pointer());
}

RefSymmSCMatrix
RefSymmSCMatrix::copy() const
{
  if (null()) return 0;
  RefSymmSCMatrix v = kit()->symmmatrix(dim());
  v.assign(*this);
  return v;
}

void
RefSymmSCMatrix::randomize() const
{
  require_nonnull();
  pointer()->randomize();
}

void
RefSymmSCMatrix::assign(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSymmSCMatrix::assign(const double*v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSymmSCMatrix::assign(const double**v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSymmSCMatrix::convert(double*v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSymmSCMatrix::convert(double**v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSymmSCMatrix::scale(double a) const
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSymmSCMatrix::scale_diagonal(double a) const
{
  require_nonnull();
  pointer()->scale_diagonal(a);
}

void
RefSymmSCMatrix::assign(double a) const
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSymmSCMatrix::accumulate(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSymmSCMatrix::element_op(const Ref<SCElementOp>&op) const
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSymmSCMatrix::element_op(const Ref<SCElementOp2>&op,
                            const RefSymmSCMatrix&m) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer());
}

void
RefSymmSCMatrix::element_op(const Ref<SCElementOp3>&op,
                            const RefSymmSCMatrix&m,
                            const RefSymmSCMatrix&n) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer(),n.pointer());
}

void
RefSymmSCMatrix::print(ostream& out) const
{
  print(0,out);
}

void
RefSymmSCMatrix::print(const char*title,ostream&out, int precision) const
{
  if (nonnull()) {
      pointer()->print(title,out,precision);
    }
  else {
      if (title) out << endl << title << endl;
      out << "null matrix" << endl;
    }
}

RefSCMatrix
RefSymmSCMatrix::operator*(const RefSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(dim(),a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSymmSCMatrix::operator*(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(dim(),a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSymmSCMatrix::operator*(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(dim(),a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCVector
RefSymmSCMatrix::operator*(const RefSCVector&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector r = kit()->vector(dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSymmSCMatrix
RefSymmSCMatrix::operator *(double a) const
{
  RefSymmSCMatrix r(copy());
  r.scale(a);
  return r;
}

RefSymmSCMatrix
operator *(double a, const RefSymmSCMatrix& v)
{
  return v*a;
}

void
RefSymmSCMatrix::accumulate_symmetric_outer_product(const RefSCVector&v) const
{
  require_nonnull();
  pointer()->accumulate_symmetric_outer_product(v.pointer());
}

double
RefSymmSCMatrix::scalar_product(const RefSCVector&v) const
{
  if (null()) return 0.0;
  return pointer()->scalar_product(v.pointer());
}

void
RefSymmSCMatrix::save(StateOut&s)
{
  if (null()) s.put(0);
  else {
      s.put(1);
      pointer()->save(s);
    }
}

void
RefSymmSCMatrix::restore(StateIn&s)
{
  int have_matrix;
  s.get(have_matrix);
  if (have_matrix && nonnull()) {
      pointer()->restore(s);
    }
  else if (have_matrix) {
      ExEnv::errn() << "RefSymmSCMatrix::restore: "
           << "matrix not properly initialized" << endl;
      abort();
    }
  else {
      clear();
    }
}

int
RefSymmSCMatrix::nblock() const
{
  BlockedSymmSCMatrix *b = dynamic_cast<BlockedSymmSCMatrix*>(pointer());
  if (b) return b->nblocks();
  return 1;
}

RefSymmSCMatrix
RefSymmSCMatrix::block(int i) const
{
  BlockedSymmSCMatrix *b = dynamic_cast<BlockedSymmSCMatrix*>(pointer());
  if (b) return b->block(i);
  return *this;
}


///////////////////////////////////////////////////////////////////
// RefDiagSCMatrix members

RefDiagSCMatrix::RefDiagSCMatrix()
{
}

RefDiagSCMatrix::RefDiagSCMatrix (const RefDiagSCMatrix & o):
  Ref<DiagSCMatrix> (o)
{
}

RefDiagSCMatrix::RefDiagSCMatrix (DiagSCMatrix * o):
  Ref<DiagSCMatrix> (o)
{
}

RefDiagSCMatrix::~RefDiagSCMatrix ()
{
}

RefDiagSCMatrix&
RefDiagSCMatrix::operator=(DiagSCMatrix* cr)
{
  Ref<DiagSCMatrix>::operator=(cr);
  return *this;
}

RefDiagSCMatrix&
RefDiagSCMatrix::operator=(const RefDiagSCMatrix & c)
{
  Ref<DiagSCMatrix>::operator=(c);
  return *this;
}

RefDiagSCMatrix::RefDiagSCMatrix(const RefSCDimension&a,
                                 const Ref<SCMatrixKit>&k)
{
  a.require_nonnull();
  assign_pointer(k->diagmatrix(a));
}

void
RefDiagSCMatrix::set_element(int i, double a) const
{
  require_nonnull();
  pointer()->set_element(i,a);
}

void
RefDiagSCMatrix::accumulate_element(int i, double a) const
{
  require_nonnull();
  pointer()->accumulate_element(i,a);
}

double
RefDiagSCMatrix::get_element(int i) const
{
  require_nonnull();
  return pointer()->get_element(i);
}

RefSCMatrix
RefDiagSCMatrix::operator*(const RefSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(dim(),a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefDiagSCMatrix::operator*(const RefSymmSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = kit()->matrix(dim(),a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator*(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix r = copy();
  Ref<SCElementOp2> op = new SCDestructiveElementProduct;
  r->element_op(op,a.pointer());
  return r;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator+(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim(),kit());

  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator-(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim(),kit());

  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::i() const
{
  require_nonnull();

  RefDiagSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::gi(double condition_number_threshold) const
{
  require_nonnull();

  RefDiagSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->gen_invert_this(condition_number_threshold);
  return ret;
}

int
RefDiagSCMatrix::n() const
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefDiagSCMatrix::dim() const
{
  if (null()) return 0;
  else return pointer()->dim();
}

Ref<SCMatrixKit>
RefDiagSCMatrix::kit() const
{
  if (null()) return 0;
  else return pointer()->kit();
}

DiagSCMatrixdouble
RefDiagSCMatrix::operator()(int i) const
{
  return DiagSCMatrixdouble(pointer(),i,i);
}

RefDiagSCMatrix
RefDiagSCMatrix::clone() const
{
  RefDiagSCMatrix r = kit()->diagmatrix(dim());
  return r;
}

RefDiagSCMatrix
RefDiagSCMatrix::copy() const
{
  if (null()) return 0;
  RefDiagSCMatrix v = kit()->diagmatrix(dim());
  v.assign(*this);
  return v;
}

void
RefDiagSCMatrix::randomize() const
{
  require_nonnull();
  pointer()->randomize();
}

void
RefDiagSCMatrix::assign(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefDiagSCMatrix::assign(const double*v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefDiagSCMatrix::convert(double*v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefDiagSCMatrix::scale(double a) const
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefDiagSCMatrix::assign(double a) const
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefDiagSCMatrix::accumulate(const RefDiagSCMatrix&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefDiagSCMatrix::element_op(const Ref<SCElementOp>&op) const
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefDiagSCMatrix::element_op(const Ref<SCElementOp2>&op,
                            const RefDiagSCMatrix&m) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer());
}

void
RefDiagSCMatrix::element_op(const Ref<SCElementOp3>&op,
                            const RefDiagSCMatrix&m,
                            const RefDiagSCMatrix&n) const
{
  if (nonnull()) pointer()->element_op(op,m.pointer(),n.pointer());
}

double
RefDiagSCMatrix::determ() const
{
  return pointer()->determ_this();
}

double
RefDiagSCMatrix::trace() const
{
  return pointer()->trace();
}

void
RefDiagSCMatrix::print(ostream& out) const
{
  print(0,out);
}

void
RefDiagSCMatrix::print(const char*title,ostream&out, int precision) const
{
  if (nonnull()) {
      pointer()->print(title,out,precision);
    }
  else {
      if (title) out << endl << title << endl;
      out << "null matrix" << endl;
    }
}

RefDiagSCMatrix
RefDiagSCMatrix::operator *(double a) const
{
  RefDiagSCMatrix r(copy());
  r.scale(a);
  return r;
}

RefDiagSCMatrix
operator *(double a, const RefDiagSCMatrix& v)
{
  return v*a;
}

void
RefDiagSCMatrix::save(StateOut&s)
{
  if (null()) s.put(0);
  else {
      s.put(1);
      pointer()->save(s);
    }
}

void
RefDiagSCMatrix::restore(StateIn&s)
{
  int have_matrix;
  s.get(have_matrix);
  if (have_matrix && nonnull()) {
      pointer()->restore(s);
    }
  else if (have_matrix) {
      ExEnv::errn() << "RefDiagSCMatrix::restore: "
           << "matrix not properly initialized" << endl;
      abort();
    }
  else {
      clear();
    }
}

int
RefDiagSCMatrix::nblock() const
{
  BlockedDiagSCMatrix *b = dynamic_cast<BlockedDiagSCMatrix*>(pointer());
  if (b) return b->nblocks();
  return 1;
}

RefDiagSCMatrix
RefDiagSCMatrix::block(int i) const
{
  BlockedDiagSCMatrix *b = dynamic_cast<BlockedDiagSCMatrix*>(pointer());
  if (b) return b->block(i);
  return *this;
}

///////////////////////////////////////////////////////////////////
// RefSCVector members

RefSCVector::RefSCVector()
{
}

RefSCVector::RefSCVector (const RefSCVector & o):
  Ref<SCVector> (o)
{
}

RefSCVector::RefSCVector (SCVector * o):
  Ref<SCVector> (o)
{
}

RefSCVector::~RefSCVector ()
{
}

RefSCVector&
RefSCVector::operator=(SCVector* cr)
{
  Ref<SCVector>::operator=(cr);
  return *this;
}

RefSCVector&
RefSCVector::operator=(const RefSCVector & c)
{
  Ref<SCVector>::operator=(c);
  return *this;
}

RefSCVector::RefSCVector(const RefSCDimension&a,
                         const Ref<SCMatrixKit>&k)
{
  a.require_nonnull();
  assign_pointer(k->vector(a));
}

void
RefSCVector::set_element(int i, double a) const
{
  require_nonnull();
  pointer()->set_element(i,a);
}

void
RefSCVector::accumulate_element(int i, double a) const
{
  require_nonnull();
  pointer()->accumulate_element(i,a);
}

double
RefSCVector::get_element(int i) const
{
  require_nonnull();
  return pointer()->get_element(i);
}

RefSCVector
RefSCVector::operator+(const RefSCVector&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector ret(dim(),kit());

  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSCVector
RefSCVector::operator-(const RefSCVector&a) const
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector ret(dim(),kit());

  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

int
RefSCVector::n() const
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefSCVector::dim() const
{
  if (null()) return 0;
  else return pointer()->dim();
}

Ref<SCMatrixKit>
RefSCVector::kit() const
{
  if (null()) return 0;
  else return pointer()->kit();
}

SCVectordouble
RefSCVector::operator()(int i) const
{
  return SCVectordouble(pointer(),i);
}

SCVectordouble
RefSCVector::operator[](int i) const
{
  return SCVectordouble(pointer(),i);
}

RefSCVector
RefSCVector::clone() const
{
  RefSCVector r = kit()->vector(dim());
  return r;
}

RefSCVector
RefSCVector::copy() const
{
  if (null()) return 0;
  RefSCVector v = kit()->vector(dim());
  v.assign(*this);
  return v;
}

double
RefSCVector::dot(const RefSCVector&a) const
{
  require_nonnull();
  return pointer()->scalar_product(a.pointer());
}

double
RefSCVector::scalar_product(const RefSCVector&a) const
{
  require_nonnull();
  return pointer()->scalar_product(a.pointer());
}

void
RefSCVector::randomize() const
{
  require_nonnull();
  pointer()->randomize();
}

void
RefSCVector::assign(const RefSCVector&a) const
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSCVector::assign(const double*v) const
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCVector::convert(double*v) const
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCVector::scale(double a) const
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSCVector::assign(double a) const
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSCVector::accumulate(const RefSCVector&a) const
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCVector::accumulate_product(const RefSymmSCMatrix&a, const RefSCVector&b)
{
  require_nonnull();
  pointer()->accumulate_product(a.pointer(), b.pointer());
}

void
RefSCVector::accumulate_product(const RefSCMatrix&a, const RefSCVector&b)
{
  require_nonnull();
  pointer()->accumulate_product(a.pointer(), b.pointer());
}

void
RefSCVector::element_op(const Ref<SCElementOp>&op) const
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSCVector::element_op(const Ref<SCElementOp2>&op,
                        const RefSCVector&v) const
{
  if (nonnull()) pointer()->element_op(op,v.pointer());
}

void
RefSCVector::element_op(const Ref<SCElementOp3>&op,
                        const RefSCVector&v,
                        const RefSCVector&w) const
{
  if (nonnull()) pointer()->element_op(op,v.pointer(),w.pointer());
}

void
RefSCVector::print(ostream& out) const
{
  print(0,out);
}

void
RefSCVector::print(const char*title,ostream&out, int precision) const
{
  if (nonnull()) {
      pointer()->print(title,out,precision);
    }
  else {
      if (title) out << endl << title << endl;
      out << "null matrix" << endl;
    }
}

RefSCVector
RefSCVector::operator *(double a) const
{
  RefSCVector r(copy());
  r.scale(a);
  return r;
}

RefSCVector
operator *(double a, const RefSCVector& v)
{
  return v*a;
}

void
RefSCVector::normalize() const
{
  require_nonnull();
  pointer()->normalize();
}

RefSymmSCMatrix
RefSCVector::symmetric_outer_product() const
{
  RefSymmSCMatrix result(dim(),kit());
  result.assign(0.0);
  result.accumulate_symmetric_outer_product(pointer());
  return result;
}

RefSCMatrix
RefSCVector::outer_product(const RefSCVector&v) const
{
  RefSCMatrix result(dim(),v.dim(),kit());
  result.assign(0.0);
  result.accumulate_outer_product(*this,v);
  return result;
}

double
RefSCVector::maxabs() const
{
  if (null()) return 0.0;
  return pointer()->maxabs();
}


void
RefSCVector::save(StateOut&s)
{
  if (null()) s.put(0);
  else {
      s.put(1);
      pointer()->save(s);
    }
}

void
RefSCVector::restore(StateIn&s)
{
  int have_matrix;
  s.get(have_matrix);
  if (have_matrix && nonnull()) {
      pointer()->restore(s);
    }
  else if (have_matrix) {
      ExEnv::errn() << "RefSCVector::restore: vector not properly initialized" << endl;
      abort();
    }
  else {
      clear();
    }
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
