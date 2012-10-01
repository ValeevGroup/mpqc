//
// elemop.cc
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

#include <stdexcept>

#include <stdlib.h>
#include <cmath>

#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/elemop.h>
#include <math/scmat/abstract.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// SCElementOp member functions

static ClassDesc SCElementOp_cd(
  typeid(SCElementOp),"SCElementOp",1,"public SavableState",
  0, 0, 0);

SCElementOp::SCElementOp()
{
}

SCElementOp::~SCElementOp()
{
}

int
SCElementOp::has_collect()
{
  return 0;
}

void
SCElementOp::defer_collect(int)
{
}

int
SCElementOp::has_side_effects()
{
  return 0;
}

void
SCElementOp::collect(const Ref<MessageGrp>&)
{
}

bool
SCElementOp::threadsafe()
{
  return false;
}

bool
SCElementOp::cloneable()
{
  return false;
}

Ref<SCElementOp>
SCElementOp::clone()
{
  throw std::runtime_error("SCElementOp::clone: not implemented");
}

void
SCElementOp::collect(const Ref<SCElementOp> &)
{
  throw std::runtime_error("SCElementOp::collect(const Ref<SCElementOp> &): "
                           "not implemented");
}

void
SCElementOp::process_base(SCMatrixBlock* a)
{
  if (dynamic_cast<SCMatrixRectBlock*>(a))
    process_spec_rect(dynamic_cast<SCMatrixRectBlock*>(a));
  else if (dynamic_cast<SCMatrixLTriBlock*>(a))
    process_spec_ltri(dynamic_cast<SCMatrixLTriBlock*>(a));
  else if (dynamic_cast<SCMatrixDiagBlock*>(a))
    process_spec_diag(dynamic_cast<SCMatrixDiagBlock*>(a));
  else if (dynamic_cast<SCVectorSimpleBlock*>(a))
    process_spec_vsimp(dynamic_cast<SCVectorSimpleBlock*>(a));
  else if (dynamic_cast<SCMatrixRectSubBlock*>(a))
    process_spec_rectsub(dynamic_cast<SCMatrixRectSubBlock*>(a));
  else if (dynamic_cast<SCMatrixLTriSubBlock*>(a))
    process_spec_ltrisub(dynamic_cast<SCMatrixLTriSubBlock*>(a));
  else if (dynamic_cast<SCMatrixDiagSubBlock*>(a))
    process_spec_diagsub(dynamic_cast<SCMatrixDiagSubBlock*>(a));
  else if (dynamic_cast<SCVectorSimpleSubBlock*>(a))
    process_spec_vsimpsub(dynamic_cast<SCVectorSimpleSubBlock*>(a));
  else
    a->process(this);
}

// If specializations of SCElementOp do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp must handle since it is pure virtual.

void
SCElementOp::process_spec_rect(SCMatrixRectBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter&r=*i;
  process(r);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
}
void
SCElementOp::process_spec_ltri(SCMatrixLTriBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process_spec_diag(SCMatrixDiagBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process_spec_vsimp(SCVectorSimpleBlock* a)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process_spec_rectsub(SCMatrixRectSubBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixRectSubBlockIter(a);
  SCMatrixBlockIter&r=*i;
  process(r);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
}
void
SCElementOp::process_spec_ltrisub(SCMatrixLTriSubBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixLTriSubBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process_spec_diagsub(SCMatrixDiagSubBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixDiagSubBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process_spec_vsimpsub(SCVectorSimpleSubBlock* a)
{
  SCMatrixBlockIter*i = new SCVectorSimpleSubBlockIter(a);
  process(*i);
  delete i;
}

/////////////////////////////////////////////////////////////////////////////
// SCElementOp2 member functions

static ClassDesc SCElementOp2_cd(
  typeid(SCElementOp2),"SCElementOp2",1,"public SavableState",
  0, 0, 0);

SCElementOp2::SCElementOp2()
{
}

SCElementOp2::~SCElementOp2()
{
}

int
SCElementOp2::has_collect()
{
  return 0;
}

void
SCElementOp2::defer_collect(int)
{
}

int
SCElementOp2::has_side_effects()
{
  return 0;
}

int
SCElementOp2::has_side_effects_in_arg()
{
  return 0;
}

void
SCElementOp2::collect(const Ref<MessageGrp>&)
{
}

void
SCElementOp2::process_base(SCMatrixBlock* a, SCMatrixBlock* b)
{
  a->process(this, b);
}

// If specializations of SCElementOp2 do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp2 must handle since it is pure virtual.

void
SCElementOp2::process_spec_rect(SCMatrixRectBlock* a,SCMatrixRectBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixRectBlockIter(b);
  process(*i,*j);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
  delete j;
}
void
SCElementOp2::process_spec_ltri(SCMatrixLTriBlock* a,SCMatrixLTriBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixLTriBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}
void
SCElementOp2::process_spec_diag(SCMatrixDiagBlock* a,SCMatrixDiagBlock* b)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixDiagBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}
void
SCElementOp2::process_spec_vsimp(SCVectorSimpleBlock* a,SCVectorSimpleBlock* b)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  SCMatrixBlockIter*j = new SCVectorSimpleBlockIter(b);
  process(*i,*j);
  delete i;
  delete j;
}

/////////////////////////////////////////////////////////////////////////////
// SCElementOp3 member functions

static ClassDesc SCElementOp3_cd(
  typeid(SCElementOp3),"SCElementOp3",1,"public SavableState",
  0, 0, 0);

SCElementOp3::SCElementOp3()
{
}

SCElementOp3::~SCElementOp3()
{
}

int
SCElementOp3::has_collect()
{
  return 0;
}

void
SCElementOp3::defer_collect(int)
{
}

int
SCElementOp3::has_side_effects()
{
  return 0;
}

int
SCElementOp3::has_side_effects_in_arg1()
{
  return 0;
}

int
SCElementOp3::has_side_effects_in_arg2()
{
  return 0;
}

void
SCElementOp3::collect(const Ref<MessageGrp>&)
{
}

void
SCElementOp3::process_base(SCMatrixBlock* a,
                           SCMatrixBlock* b,
                           SCMatrixBlock* c)
{
  a->process(this, b, c);
}

// If specializations of SCElementOp3 do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp3 must handle since it is pure virtual.

void
SCElementOp3::process_spec_rect(SCMatrixRectBlock* a,
                                SCMatrixRectBlock* b,
                                SCMatrixRectBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixRectBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixRectBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process_spec_ltri(SCMatrixLTriBlock* a,
                                SCMatrixLTriBlock* b,
                                SCMatrixLTriBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixLTriBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixLTriBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process_spec_diag(SCMatrixDiagBlock* a,
                                SCMatrixDiagBlock* b,
                                SCMatrixDiagBlock* c)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  SCMatrixBlockIter*j = new SCMatrixDiagBlockIter(b);
  SCMatrixBlockIter*k = new SCMatrixDiagBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}
void
SCElementOp3::process_spec_vsimp(SCVectorSimpleBlock* a,
                                 SCVectorSimpleBlock* b,
                                 SCVectorSimpleBlock* c)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  SCMatrixBlockIter*j = new SCVectorSimpleBlockIter(b);
  SCMatrixBlockIter*k = new SCVectorSimpleBlockIter(c);
  process(*i,*j,*k);
  delete i;
  delete j;
  delete k;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScale members

static ClassDesc SCElementScale_cd(
  typeid(SCElementScale),"SCElementScale",1,"public SCElementOp",
  0, 0, create<SCElementScale>);
SCElementScale::SCElementScale(double a):scale(a) {}
SCElementScale::SCElementScale(StateIn&s):
  SCElementOp(s)
{
  s.get(scale);
}
void
SCElementScale::save_data_state(StateOut&s)
{
  s.put(scale);
}
SCElementScale::~SCElementScale() {}
void
SCElementScale::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(scale*i.get());
    }
}

int
SCElementScale::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScalarProduct members

static ClassDesc SCElementScalarProduct_cd(
  typeid(SCElementScalarProduct),"SCElementScalarProduct",1,"public SCElementOp2",
  0, 0, create<SCElementScalarProduct>);

SCElementScalarProduct::SCElementScalarProduct():
  deferred_(0), product(0.0)
{
}

SCElementScalarProduct::SCElementScalarProduct(StateIn&s):
  SCElementOp2(s)
{
  s.get(product);
  s.get(deferred_);
}

void
SCElementScalarProduct::save_data_state(StateOut&s)
{
  s.put(product);
  s.put(deferred_);
}

SCElementScalarProduct::~SCElementScalarProduct()
{
}

void
SCElementScalarProduct::process(SCMatrixBlockIter&i,
                                SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
      product += i.get()*j.get();
    }
}

int
SCElementScalarProduct::has_collect()
{
  return 1;
}

void
SCElementScalarProduct::defer_collect(int h)
{
  deferred_=h;
}

void
SCElementScalarProduct::collect(const Ref<MessageGrp>&grp)
{
  if (!deferred_)
    grp->sum(product);
}

double
SCElementScalarProduct::result()
{
  return product;
}

/////////////////////////////////////////////////////////////////////////
// SCElementDestructiveProduct members

static ClassDesc SCElementDestructiveProduct_cd(
  typeid(SCElementDestructiveProduct),"SCElementDestructiveProduct",1,"public SCElementOp2",
  0, 0, create<SCElementDestructiveProduct>);
SCElementDestructiveProduct::SCElementDestructiveProduct() {}
SCElementDestructiveProduct::SCElementDestructiveProduct(StateIn&s):
  SCElementOp2(s)
{
}
void
SCElementDestructiveProduct::save_data_state(StateOut&s)
{
}
SCElementDestructiveProduct::~SCElementDestructiveProduct() {}
void
SCElementDestructiveProduct::process(SCMatrixBlockIter&i,
                                     SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
      i.set(i.get()*j.get());
    }
}

int
SCElementDestructiveProduct::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementDAXPY members

static ClassDesc SCElementDAXPY_cd(
  typeid(SCElementDAXPY),"SCElementDAXPY",1,"public SCElementOp2",
  0, 0, create<SCElementDAXPY>);
SCElementDAXPY::SCElementDAXPY(double a) : a_(a) {}
SCElementDAXPY::SCElementDAXPY(StateIn&s):
  SCElementOp2(s)
{
  s.get(a_);
}
void
SCElementDAXPY::save_data_state(StateOut&s)
{
  s.put(a_);
}
SCElementDAXPY::~SCElementDAXPY() {}
void
SCElementDAXPY::process(SCMatrixBlockIter&i,
                        SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
    i.set(i.get() + a_*j.get());
  }
}

int
SCElementDAXPY::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementInvert members

static ClassDesc SCElementInvert_cd(
  typeid(SCElementInvert),"SCElementInvert",1,"public SCElementOp",
  0, 0, create<SCElementInvert>);
SCElementInvert::SCElementInvert(double threshold):
  threshold_(threshold),
  nbelowthreshold_(0),
  deferred_(0)
{}
SCElementInvert::SCElementInvert(StateIn&s):
  SCElementOp(s)
{
  s.get(threshold_);
  s.get(nbelowthreshold_);
  s.get(deferred_);
}
void
SCElementInvert::save_data_state(StateOut&s)
{
  s.put(threshold_);
  s.put(nbelowthreshold_);
  s.put(deferred_);
}
SCElementInvert::~SCElementInvert() {}
void
SCElementInvert::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      double val = i.get();
      if (fabs(val) > threshold_) val = 1.0/val;
      else { val = 0.0; nbelowthreshold_++; }
      i.set(val);
    }
}

int
SCElementInvert::has_side_effects()
{
  return 1;
}
int
SCElementInvert::has_collect()
{
  return 1;
}
void
SCElementInvert::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementInvert::collect(const Ref<MessageGrp>&msg)
{
  if (!deferred_)
    msg->sum(nbelowthreshold_);
}
void
SCElementInvert::collect(const Ref<SCElementOp>&op)
{
  throw std::runtime_error(
      "SCElementInvert::collect(const Ref<SCElementOp> &): not implemented");
}

/////////////////////////////////////////////////////////////////////////
// SCElementSquareRoot members

static ClassDesc SCElementSquareRoot_cd(
  typeid(SCElementSquareRoot),"SCElementSquareRoot",1,"public SCElementOp",
  0, 0, create<SCElementSquareRoot>);
SCElementSquareRoot::SCElementSquareRoot() {}
SCElementSquareRoot::SCElementSquareRoot(double a) {}
SCElementSquareRoot::SCElementSquareRoot(StateIn&s):
  SCElementOp(s)
{
}
void
SCElementSquareRoot::save_data_state(StateOut&s)
{
}
SCElementSquareRoot::~SCElementSquareRoot() {}
void
SCElementSquareRoot::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      double val = i.get();
      if (val > 0.0) i.set(sqrt(i.get()));
      else i.set(0.0);
    }
}

int
SCElementSquareRoot::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementMaxAbs members

static ClassDesc SCElementMaxAbs_cd(
  typeid(SCElementMaxAbs),"SCElementMaxAbs",1,"public SCElementOp",
  0, 0, create<SCElementMaxAbs>);

SCElementMaxAbs::SCElementMaxAbs():deferred_(0), r(0.0) {}
SCElementMaxAbs::SCElementMaxAbs(StateIn&s):
  SCElementOp(s)
{
  s.get(r);
  s.get(deferred_);
}
void
SCElementMaxAbs::save_data_state(StateOut&s)
{
  s.put(r);
  s.put(deferred_);
}
SCElementMaxAbs::~SCElementMaxAbs() {}
void
SCElementMaxAbs::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (fabs(i.get()) > r) r = fabs(i.get());
    }
}
double
SCElementMaxAbs::result()
{
  return r;
}
int
SCElementMaxAbs::has_collect()
{
  return 1;
}
void
SCElementMaxAbs::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementMaxAbs::collect(const Ref<MessageGrp>&msg)
{
  if (!deferred_)
    msg->max(r);
}
void
SCElementMaxAbs::collect(const Ref<SCElementOp>&op)
{
  throw std::runtime_error(
      "SCElementMaxAbs::collect(const Ref<SCElementOp> &): not implemented");
}

/////////////////////////////////////////////////////////////////////////
// SCElementKNorm members

static ClassDesc SCElementKNorm_cd(
  typeid(SCElementKNorm),"SCElementKNorm",1,"public SCElementOp",
  0, 0, create<SCElementKNorm>);

SCElementKNorm::SCElementKNorm(unsigned int k):deferred_(0), r_(0.0), k_(k) {}
SCElementKNorm::SCElementKNorm(StateIn&s):
  SCElementOp(s)
{
  s.get(k_);
  s.get(r_);
  s.get(deferred_);
}
void
SCElementKNorm::save_data_state(StateOut&s)
{
  s.put(r_);
  s.put(deferred_);
}
SCElementKNorm::~SCElementKNorm() {}
void
SCElementKNorm::process(SCMatrixBlockIter&i)
{
  switch (k_) {
    case 1u:
      r_ += _process1(i); break;
    case 2u:
      r_ += _process2(i); break;
    default:
      r_ += _process(i, k_);
  }
}
double
SCElementKNorm::_process(SCMatrixBlockIter& i, int k) {
  double result = 0.0;
  for (i.reset(); i; ++i) {
    result += std::pow(std::abs(i.get()),static_cast<double>(k));
  }
  return result;
}
double
SCElementKNorm::_process1(SCMatrixBlockIter& i) {
  double result = 0.0;
  for (i.reset(); i; ++i) {
    result += std::abs(i.get());
  }
  return result;
}
double
SCElementKNorm::_process2(SCMatrixBlockIter& i) {
  double result = 0.0;
  for (i.reset(); i; ++i) {
    const double v = i.get();
    result += v * v;
  }
  return result;
}
double
SCElementKNorm::result()
{
  return r_;
}
int
SCElementKNorm::has_collect()
{
  return 1;
}
void
SCElementKNorm::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementKNorm::collect(const Ref<MessageGrp>&msg)
{
  if (!deferred_)
    msg->sum(r_);
  r_ = std::pow(r_,1.0/k_);
}
void
SCElementKNorm::collect(const Ref<SCElementOp>&op)
{
  throw std::runtime_error(
      "SCElementKNorm::collect(const Ref<SCElementOp> &): not implemented");
}

/////////////////////////////////////////////////////////////////////////
// SCElementMin members

static ClassDesc SCElementMinAbs_cd(
  typeid(SCElementMinAbs),"SCElementMinAbs",1,"public SCElementOp",
  0, 0, create<SCElementMinAbs>);
SCElementMinAbs::SCElementMinAbs(double rinit):deferred_(0), r(rinit) {}
SCElementMinAbs::SCElementMinAbs(StateIn&s):
  SCElementOp(s)
{
  s.get(r);
  s.get(deferred_);
}
void
SCElementMinAbs::save_data_state(StateOut&s)
{
  s.put(r);
  s.put(deferred_);
}
SCElementMinAbs::~SCElementMinAbs() {}
void
SCElementMinAbs::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (fabs(i.get()) < r) r = fabs(i.get());
    }
}
double
SCElementMinAbs::result()
{
  return r;
}
int
SCElementMinAbs::has_collect()
{
  return 1;
}
void
SCElementMinAbs::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementMinAbs::collect(const Ref<MessageGrp>&msg)
{
  if (!deferred_)
    msg->min(r);
}
void
SCElementMinAbs::collect(const Ref<SCElementOp>&op)
{
  throw std::runtime_error(
      "SCElementMinAbs::collect(const Ref<SCElementOp> &): not implemented");
}

/////////////////////////////////////////////////////////////////////////
// SCElementSum members

static ClassDesc SCElementSum_cd(
  typeid(SCElementSum),"SCElementSum",1,"public SCElementOp",
  0, 0, create<SCElementSum>);
SCElementSum::SCElementSum():deferred_(0), r(0.0) {}
SCElementSum::SCElementSum(StateIn&s):
  SCElementOp(s)
{
  s.get(r);
  s.get(deferred_);
}
void
SCElementSum::save_data_state(StateOut&s)
{
  s.put(r);
  s.put(deferred_);
}
SCElementSum::~SCElementSum() {}
void
SCElementSum::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      r += i.get();
    }
}
double
SCElementSum::result()
{
  return r;
}
int
SCElementSum::has_collect()
{
  return 1;
}
void
SCElementSum::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementSum::collect(const Ref<MessageGrp>&msg)
{
  if (!deferred_)
    msg->sum(r);
}
void
SCElementSum::collect(const Ref<SCElementOp>&op)
{
  throw std::runtime_error(
      "SCElementSum::collect(const Ref<SCElementOp> &): not implemented");
}

/////////////////////////////////////////////////////////////////////////
// SCElementAssign members

static ClassDesc SCElementAssign_cd(
  typeid(SCElementAssign),"SCElementAssign",1,"public SCElementOp",
  0, 0, create<SCElementAssign>);
SCElementAssign::SCElementAssign(double a):assign(a) {}
SCElementAssign::SCElementAssign(StateIn&s):
  SCElementOp(s)
{
  s.get(assign);
}
void
SCElementAssign::save_data_state(StateOut&s)
{
  s.put(assign);
}
SCElementAssign::~SCElementAssign() {}
void
SCElementAssign::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(assign);
    }
}

int
SCElementAssign::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementRandomize members

static ClassDesc SCElementRandomize_cd(
  typeid(SCElementRandomize),"SCElementRandomize",1,"public SCElementOp",
  0, 0, create<SCElementRandomize>);
SCElementRandomize::SCElementRandomize() {}
SCElementRandomize::SCElementRandomize(StateIn&s):
  SCElementOp(s)
{
}
void
SCElementRandomize::save_data_state(StateOut&s)
{
}
SCElementRandomize::~SCElementRandomize() {}
void
SCElementRandomize::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
#ifdef HAVE_DRAND48
      i.set(drand48()*(drand48()<0.5?1.0:-1.0));
#else
      int r=rand();
      double dr = (double) r / 32767.0;
      i.set(dr*(dr<0.5?1.0:-1.0));
#endif
    }
}

int
SCElementRandomize::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementShiftDiagonal members

static ClassDesc SCElementShiftDiagonal_cd(
  typeid(SCElementShiftDiagonal),"SCElementShiftDiagonal",1,"public SCElementOp",
  0, 0, create<SCElementShiftDiagonal>);
SCElementShiftDiagonal::SCElementShiftDiagonal(double a):shift_diagonal(a) {}
SCElementShiftDiagonal::SCElementShiftDiagonal(StateIn&s):
  SCElementOp(s)
{
  s.get(shift_diagonal);
}
void
SCElementShiftDiagonal::save_data_state(StateOut&s)
{
  s.put(shift_diagonal);
}
SCElementShiftDiagonal::~SCElementShiftDiagonal() {}
void
SCElementShiftDiagonal::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (i.i() == i.j()) i.set(shift_diagonal+i.get());
    }
}

int
SCElementShiftDiagonal::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementScaleDiagonal members

static ClassDesc SCElementScaleDiagonal_cd(
  typeid(SCElementScaleDiagonal),"SCElementScaleDiagonal",1,"public SCElementOp",
  0, 0, create<SCElementScaleDiagonal>);
SCElementScaleDiagonal::SCElementScaleDiagonal(double a):scale_diagonal(a) {}
SCElementScaleDiagonal::SCElementScaleDiagonal(StateIn&s):
  SCElementOp(s)
{
  s.get(scale_diagonal);
}
void
SCElementScaleDiagonal::save_data_state(StateOut&s)
{
  s.put(scale_diagonal);
}
SCElementScaleDiagonal::~SCElementScaleDiagonal() {}
void
SCElementScaleDiagonal::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (i.i() == i.j()) i.set(scale_diagonal*i.get());
    }
}

int
SCElementScaleDiagonal::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////
// SCElementDot members

static ClassDesc SCElementDot_cd(
  typeid(SCElementDot),"SCElementDot",1,"public SCElementOp",
  0, 0, create<SCElementDot>);

SCElementDot::SCElementDot(double**a, double**b, int n):
  avects(a),
  bvects(b),
  length(n)
{
}

SCElementDot::SCElementDot(StateIn&s)
{
  ExEnv::errn() << indent << "SCElementDot does not permit StateIn CTOR\n";
  abort();
}

void
SCElementDot::save_data_state(StateOut&s)
{
  ExEnv::errn() << indent << "SCElementDot does not permit save_data_state\n";
  abort();
}

int
SCElementDot::has_side_effects()
{
  return 1;
}

void
SCElementDot::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      double tmp = i.get();
      double* a = avects[i.i()];
      double* b = bvects[i.j()];
      for (int j = length; j; j--, a++, b++) {
          tmp += *a * *b;
        }
      i.accum(tmp);
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementAccumulateSCMatrix members

static ClassDesc SCElementAccumulateSCMatrix_cd(
  typeid(SCElementAccumulateSCMatrix),"SCElementAccumulateSCMatrix",1,"public SCElementOp",
  0, 0, 0);

SCElementAccumulateSCMatrix::SCElementAccumulateSCMatrix(SCMatrix*a):
  m(a)
{
}

int
SCElementAccumulateSCMatrix::has_side_effects()
{
  return 1;
}

void
SCElementAccumulateSCMatrix::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.accum(m->get_element(i.i(), i.j()));
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementAccumulateSymmSCMatrix members

static ClassDesc SCElementAccumulateSymmSCMatrix_cd(
  typeid(SCElementAccumulateSymmSCMatrix),"SCElementAccumulateSymmSCMatrix",1,"public SCElementOp",
  0, 0, 0);

SCElementAccumulateSymmSCMatrix::SCElementAccumulateSymmSCMatrix(
    SymmSCMatrix*a):
  m(a)
{
}

int
SCElementAccumulateSymmSCMatrix::has_side_effects()
{
  return 1;
}

void
SCElementAccumulateSymmSCMatrix::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.accum(m->get_element(i.i(), i.j()));
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementAccumulateDiagSCMatrix members

static ClassDesc SCElementAccumulateDiagSCMatrix_cd(
  typeid(SCElementAccumulateDiagSCMatrix),"SCElementAccumulateDiagSCMatrix",1,"public SCElementOp",
  0, 0, 0);

SCElementAccumulateDiagSCMatrix::SCElementAccumulateDiagSCMatrix(
    DiagSCMatrix*a):
  m(a)
{
}

int
SCElementAccumulateDiagSCMatrix::has_side_effects()
{
  return 1;
}

void
SCElementAccumulateDiagSCMatrix::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.accum(m->get_element(i.i()));
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementAccumulateSCVector members

static ClassDesc SCElementAccumulateSCVector_cd(
  typeid(SCElementAccumulateSCVector),"SCElementAccumulateSCVector",1,"public SCElementOp",
  0, 0, 0);

SCElementAccumulateSCVector::SCElementAccumulateSCVector(SCVector*a):
  m(a)
{
}

int
SCElementAccumulateSCVector::has_side_effects()
{
  return 1;
}

void
SCElementAccumulateSCVector::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.accum(m->get_element(i.i()));
    }
}

/////////////////////////////////////////////////////////////////////////////

SCElement SCElement::null(-1,-1,0.0);

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
