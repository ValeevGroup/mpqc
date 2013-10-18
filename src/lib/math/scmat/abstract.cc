//
// abstract.cc
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

#include <util/misc/consumableresources.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCMatrix classes.
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// SCMatrixKit members

static ClassDesc SCMatrixKit_cd(
  typeid(SCMatrixKit),"SCMatrixKit",1,"public DescribedClass",
  0, 0, 0);

SCMatrixKit::SCMatrixKit()
{
  grp_ = MessageGrp::get_default_messagegrp();
}

SCMatrixKit::SCMatrixKit(const Ref<KeyVal>& keyval)
{
  grp_ << keyval->describedclassvalue("messagegrp");
  if (grp_.null()) grp_ = MessageGrp::get_default_messagegrp();
}

SCMatrixKit::~SCMatrixKit()
{
}

SCMatrix*
SCMatrixKit::restore_matrix(StateIn& s,
                            const RefSCDimension& d1,
                            const RefSCDimension& d2)
{
  SCMatrix *r = matrix(d1,d2);
  r->restore(s);
  return r;
}

SymmSCMatrix*
SCMatrixKit::restore_symmmatrix(StateIn& s, const RefSCDimension& d)
{
  SymmSCMatrix *r = symmmatrix(d);
  r->restore(s);
  return r;
}

DiagSCMatrix*
SCMatrixKit::restore_diagmatrix(StateIn& s, const RefSCDimension& d)
{
  DiagSCMatrix *r = diagmatrix(d);
  r->restore(s);
  return r;
}

SCVector*
SCMatrixKit::restore_vector(StateIn& s, const RefSCDimension& d)
{
  SCVector *r = vector(d);
  r->restore(s);
  return r;
}

Ref<MessageGrp>
SCMatrixKit::messagegrp() const
{
  return grp_;
}

/////////////////////////////////////////////////////////////////////////
// SCMatrix members

static ClassDesc SCMatrix_cd(
  typeid(SCMatrix),"SCMatrix",1,"public DescribedClass",
  0, 0, 0);

SCMatrix::SCMatrix(const RefSCDimension&a1, const RefSCDimension&a2,
                   SCMatrixKit*kit):
  d1(a1),
  d2(a2),
  kit_(kit)
{
}

SCMatrix::~SCMatrix()
{
}

void
SCMatrix::save(StateOut&s)
{
  int nr = nrow();
  int nc = ncol();
  s.put(nr);
  s.put(nc);
  int has_subblocks = 0;
  s.put(has_subblocks);
  for (int i=0; i<nr; i++) {
      for (int j=0; j<nc; j++) {
          s.put(get_element(i,j));
        }
    }
}

void
SCMatrix::restore(StateIn& s)
{
  int nrt, nr = nrow();
  int nct, nc = ncol();
  s.get(nrt);
  s.get(nct);
  if (nrt != nr || nct != nc) {
      ExEnv::errn() << "SCMatrix::restore(): bad dimensions" << endl;
      abort();
    }
  int has_subblocks;
  s.get(has_subblocks);
  if (!has_subblocks) {
      for (int i=0; i<nr; i++) {
          for (int j=0; j<nc; j++) {
              double tmp;
              s.get(tmp);
              set_element(i,j, tmp);
            }
        }
    }
  else {
      ExEnv::errn() << "SCMatrix::restore(): matrix has subblocks--cannot restore"
           << endl;
      abort();
    }
}

double
SCMatrix::maxabs() const
{
  Ref<SCElementMaxAbs> op = new SCElementMaxAbs();
  Ref<SCElementOp> abop = op.pointer();
  ((SCMatrix *)this)->element_op(abop);
  return op->result();
}

void
SCMatrix::randomize()
{
  Ref<SCElementOp> op = new SCElementRandomize();
  this->element_op(op);
}

void
SCMatrix::assign_val(double a)
{
  Ref<SCElementOp> op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCMatrix::scale(double a)
{
  Ref<SCElementOp> op = new SCElementScale(a);
  this->element_op(op);
}

void
SCMatrix::scale_diagonal(double a)
{
  Ref<SCElementOp> op = new SCElementScaleDiagonal(a);
  this->element_op(op);
}

void
SCMatrix::shift_diagonal(double a)
{
  Ref<SCElementOp> op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SCMatrix::assign_r(SCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCMatrix::assign_p(const double*a)
{
  int i;
  int nr = nrow();
  int nc = ncol();
  // some compilers need the following cast
  const double **v = (const double**) new double*[nr];
  for (i=0; i<nr; i++) {
      v[i] = &a[i*nc];
    }
  assign(v);
  delete[] v;
}

void
SCMatrix::assign_pp(const double**a)
{
  int i;
  int j;
  int nr;
  int nc;
  nr = nrow();
  nc = ncol();
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SCMatrix::convert_r(SCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SCMatrix::convert_accumulate(SCMatrix*a)
{
  Ref<SCElementOp> op = new SCElementAccumulateSCMatrix(a);
  element_op(op);
}

void
SCMatrix::convert_p(double*a) const
{
  int i;
  int nr = nrow();
  int nc = ncol();
  double **v = new double*[nr];
  for (i=0; i<nr; i++) {
      v[i] = &a[i*nc];
    }
  convert(v);
  delete[] v;
}

void
SCMatrix::convert_pp(double**a) const
{
  int i, j;
  int nr, nc;
  nr = nrow();
  nc = ncol();
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          a[i][j] = get_element(i,j);
        }
    }
}

void
SCMatrix::accumulate_product_sr(SymmSCMatrix*a,SCMatrix*b)
{
  SCMatrix *t = b->copy();
  t->transpose_this();
  SCMatrix *t2 = this->copy();
  t2->transpose_this();
  t2->accumulate_product(t,a);
  delete t;
  t2->transpose_this();
  assign(t2);
  delete t2;
}

void
SCMatrix::accumulate_product_dr(DiagSCMatrix*a,SCMatrix*b)
{
  SCMatrix *t = b->copy();
  t->transpose_this();
  SCMatrix *t2 = kit()->matrix(coldim(),rowdim());
  t2->assign(0.0);
  t2->accumulate_product(t,a);
  delete t;
  t2->transpose_this();
  accumulate(t2);
  delete t2;
}

void
SCMatrix::print(ostream&o) const
{
  vprint(0, o, 10);
}

void
SCMatrix::print(const char *t, ostream&o, int i) const
{
  vprint(t, o, i);
}

SCMatrix*
SCMatrix::clone()
{
  return kit()->matrix(rowdim(),coldim());
}

SCMatrix*
SCMatrix::copy()
{
  SCMatrix* result = clone();
  result->assign(this);
  return result;
}

void
SCMatrix::gen_invert_this(double condition_number_threshold)
{
  int i;

  // Compute the singular value decomposition of this
  RefSCMatrix U(rowdim(),rowdim(),kit());
  RefSCMatrix V(coldim(),coldim(),kit());
  RefSCDimension min;
  if (coldim().n()<rowdim().n()) min = coldim();
  else min = rowdim();
  int nmin = min.n();
  RefDiagSCMatrix sigma(min,kit());
  svd_this(U.pointer(),sigma.pointer(),V.pointer());

  // compute the epsilon rank of this
  const double sigma_max = sigma->maxabs();
  const double sigma_min_threshold = sigma_max / condition_number_threshold;
  int rank = 0;
  for (i=0; i<nmin; i++) {
      if (sigma(i) > sigma_min_threshold) rank++;
    }

  RefSCDimension drank = new SCDimension(rank);
  RefDiagSCMatrix sigma_i(drank,kit());
  for (i=0; i<rank; i++) {
      sigma_i(i) = 1.0/sigma(i);
    }
  RefSCMatrix Ur(rowdim(), drank, kit());
  RefSCMatrix Vr(coldim(), drank, kit());
  Ur.assign_subblock(U,0, rowdim().n()-1, 0, drank.n()-1, 0, 0);
  Vr.assign_subblock(V,0, coldim().n()-1, 0, drank.n()-1, 0, 0);
  assign((Vr * sigma_i * Ur.t()).t());
  transpose_this();
}

void
SCMatrix::svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V)
{
  ExEnv::errn() << indent << class_name() << ": SVD not implemented\n";
  abort();
}

void
SCMatrix::accumulate_product_rs(SCMatrix*a,SymmSCMatrix*b)
{
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

void
SCMatrix::accumulate_product_ss(SymmSCMatrix*a,SymmSCMatrix*b)
{
  SCMatrix *arect = kit()->matrix(a->dim(),a->dim());
  arect->assign(0.0);
  arect->accumulate(a);
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(arect,brect);
  delete arect;
  delete brect;
}

void
SCMatrix::accumulate_product_sd(SymmSCMatrix*a,DiagSCMatrix*b)
{
  SCMatrix *arect = kit()->matrix(a->dim(),a->dim());
  arect->assign(0.0);
  arect->accumulate(a);
  accumulate_product(arect,b);
  delete arect;
}

void
SCMatrix::accumulate_product_ds(DiagSCMatrix*a,SymmSCMatrix*b)
{
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

void
SCMatrix::accumulate_product_rd(SCMatrix*a,DiagSCMatrix*b)
{
  SCMatrix *brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);
  accumulate_product(a,brect);
  delete brect;
}

Ref<MessageGrp>
SCMatrix::messagegrp() const
{
  return kit_->messagegrp();
}

/////////////////////////////////////////////////////////////////////////
// SymmSCMatrix member functions

static ClassDesc SymmSCMatrix_cd(
  typeid(SymmSCMatrix),"SymmSCMatrix",1,"public DescribedClass",
  0, 0, 0);

SymmSCMatrix::SymmSCMatrix(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
{
}

SymmSCMatrix::~SymmSCMatrix()
{
}

void
SymmSCMatrix::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; j++) {
          s.put(get_element(i,j));
        }
    }
}

void
SymmSCMatrix::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      ExEnv::errn() << "SymmSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      for (int j=0; j<=i; j++) {
          double tmp;
          s.get(tmp);
          set_element(i,j, tmp);
        }
    }
}

double
SymmSCMatrix::maxabs() const
{
  Ref<SCElementMaxAbs> op = new SCElementMaxAbs();
  Ref<SCElementOp> abop = op.pointer();
  ((SymmSCMatrix*)this)->element_op(abop);
  return op->result();
}

void
SymmSCMatrix::randomize()
{
  Ref<SCElementOp> op = new SCElementRandomize();
  this->element_op(op);
}

void
SymmSCMatrix::assign_val(double a)
{
  Ref<SCElementOp> op = new SCElementAssign(a);
  this->element_op(op);
}

void
SymmSCMatrix::assign_p(const double*a)
{
  int i;
  int nr = n();
  // some compilers need the following cast
  const double **v = (const double **) new double*[nr];
  int ioff= 0;
  for (i=0; i<nr; i++) {
      ioff += i;
      v[i] = &a[ioff];
//    ioff += i;
    }
  assign(v);
  delete[] v;
}

void
SymmSCMatrix::assign_pp(const double**a)
{
  int i;
  int j;
  int nr;
  nr = n();
  for (i=0; i<nr; i++) {
      for (j=0; j<=i; j++) {
          set_element(i,j,a[i][j]);
        }
    }
}

void
SymmSCMatrix::convert_s(SymmSCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SymmSCMatrix::convert_accumulate(SymmSCMatrix*a)
{
  Ref<SCElementOp> op = new SCElementAccumulateSymmSCMatrix(a);
  element_op(op);
}

void
SymmSCMatrix::convert_p(double*a) const
{
  int i;
  int nr = n();
  double **v = new double*[nr];
  int ioff= 0;
  for (i=0; i<nr; i++) {
      ioff += i;
      v[i] = &a[ioff];
//    ioff += i;
    }
  convert(v);
  delete[] v;
}

void
SymmSCMatrix::convert_pp(double**a) const
{
  int i;
  int j;
  int nr;
  nr = n();
  for (i=0; i<nr; i++) {
      for (j=0; j<=i; j++) {
          a[i][j] = get_element(i,j);
        }
    }
}

void
SymmSCMatrix::scale(double a)
{
  Ref<SCElementOp> op = new SCElementScale(a);
  this->element_op(op);
}

void
SymmSCMatrix::scale_diagonal(double a)
{
  Ref<SCElementOp> op = new SCElementScaleDiagonal(a);
  this->element_op(op);
}

void
SymmSCMatrix::shift_diagonal(double a)
{
  Ref<SCElementOp> op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SymmSCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SymmSCMatrix::assign_s(SymmSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SymmSCMatrix::print(ostream&o) const
{
  vprint(0, o, 10);
}

void
SymmSCMatrix::print(const char *t, ostream&o, int i) const
{
  vprint(t, o, i);
}

void
SymmSCMatrix::vprint(const char* title, ostream& out, int i) const
{
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
}

SymmSCMatrix*
SymmSCMatrix::clone()
{
  return kit()->symmmatrix(dim());
}

SymmSCMatrix*
SymmSCMatrix::copy()
{
  SymmSCMatrix* result = clone();
  result->assign(this);
  return result;
}

void
SymmSCMatrix::accumulate_symmetric_product(SCMatrix *a)
{
  RefSCMatrix at = a->copy();
  at->transpose_this();
  RefSCMatrix m = kit()->matrix(a->rowdim(),a->rowdim());
  m->assign(0.0);
  m->accumulate_product(a, at.pointer());
  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SCMatrix *a, SymmSCMatrix *b,
                                   SCMatrix::Transform t)
{
  RefSCMatrix brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix res;

  if (t == SCMatrix::TransposeTransform) { // return A^T*B*A
      RefSCMatrix at = a->copy();
      at->transpose_this();

      RefSCMatrix tmp = at->clone();
      tmp->assign(0.0);

      tmp->accumulate_product(at.pointer(), brect.pointer());
      brect = 0;
      at = 0;

      res = kit()->matrix(dim(),dim());
      res->assign(0.0);
      res->accumulate_product(tmp.pointer(), a);
    }
  else { // return A *B* A^T
      RefSCMatrix tmp = a->clone();
      tmp->assign(0.0);

      tmp->accumulate_product(a,brect);
      brect = 0;

      RefSCMatrix at = a->copy();
      at->transpose_this();

      res = kit()->matrix(dim(),dim());
      res->assign(0.0);
      res->accumulate_product(tmp.pointer(), at.pointer());
      at = 0;
    }

  scale(2.0);
  accumulate_symmetric_sum(res.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SCMatrix *a, DiagSCMatrix *b,
                                   SCMatrix::Transform t)
{
  RefSCMatrix m = kit()->matrix(a->rowdim(),a->rowdim());
  RefSCMatrix brect = kit()->matrix(b->dim(),b->dim());
  brect->assign(0.0);
  brect->accumulate(b);

  RefSCMatrix tmp = a->clone();
  tmp->assign(0.0);

  RefSCMatrix res;

  if (t == SCMatrix::TransposeTransform) {
      RefSCMatrix at = a->copy();
      at->transpose_this();

      tmp->accumulate_product(at.pointer(), brect.pointer());
      brect = 0;
      at = 0;

      res = kit()->matrix(dim(),dim());
      res->assign(0.0);
      res->accumulate_product(tmp.pointer(), a);
    }
  else {
      tmp->accumulate_product(a, brect.pointer());
      brect = 0;

      RefSCMatrix at = a->copy();
      at->transpose_this();

      res = kit()->matrix(dim(),dim());
      res->assign(0.0);
      res->accumulate_product(tmp.pointer(), at.pointer());
      at = 0;
    }

  tmp = 0;

  scale(2.0);
  accumulate_symmetric_sum(res.pointer());
  scale(0.5);
}

void
SymmSCMatrix::accumulate_transform(SymmSCMatrix *a, SymmSCMatrix *b)
{
  RefSCMatrix m = kit()->matrix(a->dim(),a->dim());
  m->assign(0.0);
  m->accumulate(a);
  accumulate_transform(m.pointer(),b);
}

void
SymmSCMatrix::accumulate_symmetric_outer_product(SCVector *v)
{
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate_outer_product(v,v);

  scale(2.0);
  accumulate_symmetric_sum(m.pointer());
  scale(0.5);
}

double
SymmSCMatrix::scalar_product(SCVector *v)
{
  RefSCVector v2 = kit()->vector(dim());
  v2->assign(0.0);
  v2->accumulate_product(this,v);
  return v2->scalar_product(v);
}

Ref<MessageGrp>
SymmSCMatrix::messagegrp() const
{
  return kit_->messagegrp();
}

/////////////////////////////////////////////////////////////////////////
// DiagSCMatrix member functions

static ClassDesc DiagSCMatrix_cd(
  typeid(DiagSCMatrix),"DiagSCMatrix",1,"public DescribedClass",
  0, 0, 0);

DiagSCMatrix::DiagSCMatrix(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
{
}

DiagSCMatrix::~DiagSCMatrix()
{
}

void
DiagSCMatrix::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      s.put(get_element(i));
    }
}

void
DiagSCMatrix::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      ExEnv::errn() << "DiagSCMatrix::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      double tmp;
      s.get(tmp);
      set_element(i, tmp);
    }
}

double
DiagSCMatrix::maxabs() const
{
  Ref<SCElementMaxAbs> op = new SCElementMaxAbs();
  Ref<SCElementOp> abop = op.pointer();
  ((DiagSCMatrix*)this)->element_op(abop);
  return op->result();
}

void
DiagSCMatrix::randomize()
{
  Ref<SCElementOp> op = new SCElementRandomize();
  this->element_op(op);
}

void
DiagSCMatrix::assign_val(double a)
{
  Ref<SCElementOp> op = new SCElementAssign(a);
  this->element_op(op);
}

void
DiagSCMatrix::assign_p(const double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
DiagSCMatrix::convert_d(DiagSCMatrix*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
DiagSCMatrix::convert_accumulate(DiagSCMatrix*a)
{
  Ref<SCElementOp> op = new SCElementAccumulateDiagSCMatrix(a);
  element_op(op);
}

void
DiagSCMatrix::convert_p(double*a) const
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
DiagSCMatrix::scale(double a)
{
  Ref<SCElementOp> op = new SCElementScale(a);
  this->element_op(op);
}

void
DiagSCMatrix::assign_d(DiagSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
DiagSCMatrix::print(ostream&o) const
{
  vprint(0, o, 10);
}

void
DiagSCMatrix::print(const char *t, ostream&o, int i) const
{
  vprint(t, o, i);
}

void
DiagSCMatrix::vprint(const char* title, ostream& out, int i) const
{
  RefSCMatrix m = kit()->matrix(dim(),dim());
  m->assign(0.0);
  m->accumulate(this);
  m->print(title, out, i);
}

DiagSCMatrix*
DiagSCMatrix::clone()
{
  return kit()->diagmatrix(dim());
}

DiagSCMatrix*
DiagSCMatrix::copy()
{
  DiagSCMatrix* result = clone();
  result->assign(this);
  return result;
}

Ref<MessageGrp>
DiagSCMatrix::messagegrp() const
{
  return kit_->messagegrp();
}

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCVector classes.
/////////////////////////////////////////////////////////////////////////

static ClassDesc SCVector_cd(
  typeid(SCVector),"SCVector",1,"public DescribedClass",
  0, 0, 0);

SCVector::SCVector(const RefSCDimension&a, SCMatrixKit *kit):
  d(a),
  kit_(kit)
{
}

SCVector::~SCVector()
{
}

void
SCVector::save(StateOut&s)
{
  int nr = n();
  s.put(nr);
  for (int i=0; i<nr; i++) {
      s.put(get_element(i));
    }
}

void
SCVector::write_xml(boost::property_tree::ptree& pt)
{
  pt.put("SCVector.<xmlattr>.n", n());
  pt.put("SCVector.data.<xmlattr>.compressed", true);

  double* data = allocate<double>(n());
  this->convert(data);
  XMLDataStream<double> xds(data, n());
  pt.put("SCVector.data", xds);
}


void
SCVector::restore(StateIn& s)
{
  int nrt, nr = n();
  s.get(nrt);
  if (nrt != nr) {
      ExEnv::errn() << "SCVector::restore(): bad dimension" << endl;
      abort();
    }
  for (int i=0; i<nr; i++) {
      double tmp;
      s.get(tmp);
      set_element(i, tmp);
    }
}

double
SCVector::maxabs() const
{
  Ref<SCElementMaxAbs> op = new SCElementMaxAbs();
  Ref<SCElementOp> abop = op.pointer();
  ((SCVector*)this)->element_op(abop);
  return op->result();
}

void
SCVector::randomize()
{
  Ref<SCElementOp> op = new SCElementRandomize();
  this->element_op(op);
}

void
SCVector::assign_val(double a)
{
  Ref<SCElementOp> op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCVector::assign_p(const double*a)
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      set_element(i,a[i]);
    }
}

void
SCVector::convert_v(SCVector*a)
{
  assign(0.0);
  convert_accumulate(a);
}

void
SCVector::convert_accumulate(SCVector*a)
{
  Ref<SCElementOp> op = new SCElementAccumulateSCVector(a);
  element_op(op);
}

void
SCVector::convert_p(double*a) const
{
  int i;
  int nr = n();
  for (i=0; i<nr; i++) {
      a[i] = get_element(i);
    }
}

void
SCVector::scale(double a)
{
  Ref<SCElementOp> op = new SCElementScale(a);
  this->element_op(op);
}

void
SCVector::assign_v(SCVector*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCVector::print(ostream&o) const
{
  vprint(0, o, 10);
}

void
SCVector::print(const char *t, ostream&o, int i) const
{
  vprint(t, o, i);
}

void
SCVector::normalize()
{
  double norm = sqrt(scalar_product(this));
  if (norm > 1.e-20) norm = 1.0/norm;
  else {
      ExEnv::errn() << indent
           << "SCVector::normalize: tried to normalize tiny vector\n";
      abort();
    }
  scale(norm);
}

SCVector*
SCVector::clone()
{
  return kit()->vector(dim());
}

SCVector*
SCVector::copy()
{
  SCVector* result = clone();
  result->assign(this);
  return result;
}

void
SCVector::accumulate_product_sv(SymmSCMatrix *m, SCVector *v)
{
  RefSCMatrix mrect = kit()->matrix(m->dim(),m->dim());
  mrect->assign(0.0);
  mrect->accumulate(m);
  accumulate_product(mrect.pointer(), v);
}

Ref<MessageGrp>
SCVector::messagegrp() const
{
  return kit_->messagegrp();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
