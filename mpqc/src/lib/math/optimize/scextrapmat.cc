//
// scextrapmat.cc
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

#include <util/state/stateio.h>
#include <math/scmat/elemop.h>
#include <math/scmat/blocked.h>
#include <math/optimize/scextrapmat.h>

using namespace sc;

static ClassDesc SymmSCMatrixSCExtrapData_cd(
  typeid(SymmSCMatrixSCExtrapData),"SymmSCMatrixSCExtrapData",1,"public SCExtrapData",
  0, 0, create<SymmSCMatrixSCExtrapData>);

SymmSCMatrixSCExtrapData::SymmSCMatrixSCExtrapData(StateIn& s) :
  SCExtrapData(s)
{
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m = k->symmmatrix(dim);
  m.restore(s);
}

SymmSCMatrixSCExtrapData::SymmSCMatrixSCExtrapData(const RefSymmSCMatrix& mat)
{
  m = mat;
}

void
SymmSCMatrixSCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);
  SavableState::save_state(m.dim().pointer(),s);

  int blocked = (dynamic_cast<BlockedSymmSCMatrix*>(m.pointer())) ? 1 : 0;
  s.put(blocked);
  
  m.save(s);
}

void
SymmSCMatrixSCExtrapData::zero()
{
  m.assign(0.0);
}

SCExtrapData*
SymmSCMatrixSCExtrapData::copy()
{
  return new SymmSCMatrixSCExtrapData(m.copy());
}

void
SymmSCMatrixSCExtrapData::accumulate_scaled(double scale,
                                            const Ref<SCExtrapData>& data)
{
  SymmSCMatrixSCExtrapData* a
      = require_dynamic_cast<SymmSCMatrixSCExtrapData*>(
          data.pointer(), "SymmSCMatrixSCExtrapData::accumulate_scaled");

  RefSymmSCMatrix am = a->m.copy();
  am.scale(scale);
  m.accumulate(am);
}

///////////////////////////////////////////////////////////////////////////

static ClassDesc SymmSCMatrix2SCExtrapData_cd(
  typeid(SymmSCMatrix2SCExtrapData),"SymmSCMatrix2SCExtrapData",1,"public SCExtrapData",
  0, 0, create<SymmSCMatrix2SCExtrapData>);

SymmSCMatrix2SCExtrapData::SymmSCMatrix2SCExtrapData(StateIn&s) :
  SCExtrapData(s)
{
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m1 = k->symmmatrix(dim);
  m2 = k->symmmatrix(dim);
  m1.restore(s);
  m2.restore(s);
}

SymmSCMatrix2SCExtrapData::SymmSCMatrix2SCExtrapData(
    const RefSymmSCMatrix& mat1,
    const RefSymmSCMatrix& mat2)
{
  m1 = mat1;
  m2 = mat2;
}

void
SymmSCMatrix2SCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);
  SavableState::save_state(m1.dim().pointer(),s);

  int blocked = (dynamic_cast<BlockedSymmSCMatrix*>(m1.pointer())) ? 1 : 0;
  s.put(blocked);
  
  m1.save(s);
  m2.save(s);
}

void
SymmSCMatrix2SCExtrapData::zero()
{
  m1.assign(0.0);
  m2.assign(0.0);
}

SCExtrapData*
SymmSCMatrix2SCExtrapData::copy()
{
  return new SymmSCMatrix2SCExtrapData(m1.copy(), m2.copy());
}

void
SymmSCMatrix2SCExtrapData::accumulate_scaled(double scale,
                                             const Ref<SCExtrapData>& data)
{
  SymmSCMatrix2SCExtrapData* a
      = require_dynamic_cast<SymmSCMatrix2SCExtrapData*>(
          data.pointer(), "SymmSCMatrix2SCExtrapData::accumulate_scaled");

  RefSymmSCMatrix am = a->m1.copy();
  am.scale(scale);
  m1.accumulate(am);
  am = 0;

  am = a->m2.copy();
  am.scale(scale);
  m2.accumulate(am);
}

///////////////////////////////////////////////////////////////////////////

static ClassDesc SymmSCMatrix4SCExtrapData_cd(
  typeid(SymmSCMatrix4SCExtrapData),"SymmSCMatrix4SCExtrapData",1,"public SCExtrapData",
  0, 0, create<SymmSCMatrix4SCExtrapData>);

SymmSCMatrix4SCExtrapData::SymmSCMatrix4SCExtrapData(StateIn&s) :
  SCExtrapData(s)
{
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m1 = k->symmmatrix(dim);
  m2 = k->symmmatrix(dim);
  m3 = k->symmmatrix(dim);
  m4 = k->symmmatrix(dim);
  m1.restore(s);
  m2.restore(s);
  m3.restore(s);
  m4.restore(s);
}

SymmSCMatrix4SCExtrapData::SymmSCMatrix4SCExtrapData(
    const RefSymmSCMatrix& mat1,
    const RefSymmSCMatrix& mat2,
    const RefSymmSCMatrix& mat3,
    const RefSymmSCMatrix& mat4)
{
  m1 = mat1;
  m2 = mat2;
  m3 = mat3;
  m4 = mat4;
}

void
SymmSCMatrix4SCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);
  SavableState::save_state(m1.dim().pointer(),s);

  int blocked = (dynamic_cast<BlockedSymmSCMatrix*>(m1.pointer())) ? 1 : 0;
  s.put(blocked);
  
  m1.save(s);
  m2.save(s);
  m3.save(s);
  m4.save(s);
}

void
SymmSCMatrix4SCExtrapData::zero()
{
  m1.assign(0.0);
  m2.assign(0.0);
  m3.assign(0.0);
  m4.assign(0.0);
}

SCExtrapData*
SymmSCMatrix4SCExtrapData::copy()
{
  return new SymmSCMatrix4SCExtrapData(m1.copy(), m2.copy(),
                                       m3.copy(), m4.copy());
}

void
SymmSCMatrix4SCExtrapData::accumulate_scaled(double scale,
                                             const Ref<SCExtrapData>& data)
{
  SymmSCMatrix4SCExtrapData* a
      = require_dynamic_cast<SymmSCMatrix4SCExtrapData*>(
          data.pointer(), "SymmSCMatrix4SCExtrapData::accumulate_scaled");

  RefSymmSCMatrix am = a->m1.copy();
  am.scale(scale);
  m1.accumulate(am);
  am = 0;

  am = a->m2.copy();
  am.scale(scale);
  m2.accumulate(am);

  am = a->m3.copy();
  am.scale(scale);
  m3.accumulate(am);

  am = a->m4.copy();
  am.scale(scale);
  m4.accumulate(am);
}

///////////////////////////////////////////////////////////////////////////

static ClassDesc SymmSCMatrixNSCExtrapData_cd(
  typeid(SymmSCMatrixNSCExtrapData),"SymmSCMatrixNSCExtrapData",1,"public SCExtrapData",
  0, 0, create<SymmSCMatrixNSCExtrapData>);

SymmSCMatrixNSCExtrapData::SymmSCMatrixNSCExtrapData(StateIn&s) :
  SCExtrapData(s)
{
  s.get(n_);
  
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m = new RefSymmSCMatrix[n_];
  
  for (int i=0; i < n_; i++) {
    m[i] = k->symmmatrix(dim);
    m[i].restore(s);
  }
}

SymmSCMatrixNSCExtrapData::SymmSCMatrixNSCExtrapData(int n,
                                                     RefSymmSCMatrix *mats)
{
  n_=n;
  m = new RefSymmSCMatrix[n_];
  for (int i=0; i < n_; i++)
    m[i] = mats[i];
}

void
SymmSCMatrixNSCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);

  s.put(n_);
  SavableState::save_state(m[0].dim().pointer(),s);

  int blocked = (dynamic_cast<BlockedSymmSCMatrix*>(m[0].pointer())) ? 1 : 0;
  s.put(blocked);
  
  for (int i=0; i < n_; i++)
    m[i].save(s);
}

void
SymmSCMatrixNSCExtrapData::zero()
{
  for (int i=0; i < n_; i++)
    m[i].assign(0.0);
}

SCExtrapData*
SymmSCMatrixNSCExtrapData::copy()
{
  RefSymmSCMatrix *m2 = new RefSymmSCMatrix[n_];
  for (int i=0; i < n_; i++)
    m2[i] = m[i].copy();
  
  SCExtrapData *ret = new SymmSCMatrixNSCExtrapData(n_, m2);
  delete[] m2;

  return ret;
}

void
SymmSCMatrixNSCExtrapData::accumulate_scaled(double scale,
                                             const Ref<SCExtrapData>& data)
{
  SymmSCMatrixNSCExtrapData* a
      = require_dynamic_cast<SymmSCMatrixNSCExtrapData*>(
          data.pointer(), "SymmSCMatrixNSCExtrapData::accumulate_scaled");

  for (int i=0; i < n_; i++) {
    RefSymmSCMatrix am = a->m[i].copy();
    am.scale(scale);
    m[i].accumulate(am);
  }
}

///////////////////////////////////////////////////////////////////////////

static ClassDesc SymmSCMatrixSCExtrapError_cd(
  typeid(SymmSCMatrixSCExtrapError),"SymmSCMatrixSCExtrapError",1,"public SCExtrapError",
  0, 0, create<SymmSCMatrixSCExtrapError>);

SymmSCMatrixSCExtrapError::SymmSCMatrixSCExtrapError(StateIn& s) :
  SCExtrapError(s)
{
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m = k->symmmatrix(dim);
  m.restore(s);
}

SymmSCMatrixSCExtrapError::SymmSCMatrixSCExtrapError(
    const RefSymmSCMatrix& mat)
{
  m = mat;
}

void
SymmSCMatrixSCExtrapError::save_data_state(StateOut& s)
{
  SCExtrapError::save_data_state(s);
  SavableState::save_state(m.dim().pointer(),s);

  int blocked = (dynamic_cast<BlockedSymmSCMatrix*>(m.pointer())) ? 1 : 0;
  s.put(blocked);
  
  m.save(s);
}

double
SymmSCMatrixSCExtrapError::error()
{
  return m->maxabs();
}

double
SymmSCMatrixSCExtrapError::scalar_product(const Ref<SCExtrapError>& arg)
{
  SymmSCMatrixSCExtrapError* a
      = require_dynamic_cast<SymmSCMatrixSCExtrapError*>(
          arg.pointer(), "SymmSCMatrixSCExtrapError::scalar_product");
  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  m->element_op(sp.pointer(), a->m.pointer());
  return sp->result();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
