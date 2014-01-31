//
// intparams.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <float.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/intparams.h>

using namespace sc;

///////////////////////////////////////////////////////////////////////

EfieldDotVectorData::~EfieldDotVectorData()
{
}

void
EfieldDotVectorData::set_position(double*p)
{
  position[0] = p[0];
  position[1] = p[1];
  position[2] = p[2];
}

void
EfieldDotVectorData::set_vector(double*v)
{
  vector[0] = v[0];
  vector[1] = v[1];
  vector[2] = v[2];
}

///////////////////////////////////////////////////////////////////////

DipoleData::~DipoleData()
{
}

void
DipoleData::set_origin(double*o)
{
  origin[0] = o[0];
  origin[1] = o[1];
  origin[2] = o[2];
}

///////////////////////////////////////////////////////////////////////

PointChargeData::PointChargeData(int ncharges,
                                 const double *const*positions,
                                 const double *charges,
                                 int copy_data)
{
  ncharges_ = ncharges;
  if (copy_data) {
    alloced_positions_ = new double*[ncharges];
    alloced_charges_ = new double[ncharges];
    memcpy(alloced_charges_, charges, sizeof(double)*ncharges);
    double *tmp = new double[ncharges*3];
    for (int i=0; i<ncharges; i++) {
      alloced_positions_[i] = tmp;
      for (int j=0; j<3; j++) {
        *tmp++ = positions[i][j];
      }
    }
    positions_ = alloced_positions_;
    charges_ = alloced_charges_;
  }
  else {
    charges_ = charges;
    alloced_charges_ = 0;
    alloced_positions_ = 0;
    charges_ = charges;
    positions_ = positions;
  }
}

PointChargeData::~PointChargeData()
{
  if (alloced_positions_) {
    delete[] alloced_positions_[0];
    delete[] alloced_positions_;
  }
  delete[] alloced_charges_;
}

/////////////////////////////

static ClassDesc IntParams_cd(
  typeid(IntParams),"IntParams",1,"virtual public SavableState",
  0, 0, 0);

IntParams::IntParams(unsigned int nparams) : nparams_(nparams) {}

IntParams::IntParams(StateIn& si) : SavableState(si)
{
  si.get(nparams_);
}

IntParams::~IntParams() {}

void
IntParams::save_data_state(StateOut& so)
{
  so.put(nparams_);
}

unsigned int
IntParams::nparams() const
{
  return nparams_;
}

/////////////////////////////

static ClassDesc IntParamsVoid_cd(
  typeid(IntParamsVoid),"IntParamsVoid",1,"public IntParams",
  create<IntParamsVoid>, 0, create<IntParamsVoid>);

IntParamsVoid::IntParamsVoid() : IntParams(0) {}

IntParamsVoid::IntParamsVoid(StateIn& si) : IntParams(si) {}

IntParamsVoid::~IntParamsVoid() {}

void
IntParamsVoid::save_data_state(StateOut& so)
{
  IntParams::save_data_state(so);
}

bool
IntParamsVoid::equiv(const IntParams& other) const {
  return (downcast<IntParamsVoid>(other) != 0);
}

/////////////////////////////

static ClassDesc IntParamsOrigin_cd(
  typeid(IntParamsOrigin),"IntParamsOrigin",1,"public IntParams",
  create<IntParamsOrigin>, 0, create<IntParamsOrigin>);

IntParamsOrigin::IntParamsOrigin() : IntParams(0), O_(3, 0.0) {
}

IntParamsOrigin::IntParamsOrigin(const double (&O)[3]) : IntParams(0), O_(3) {
  std::copy(O, O+3, O_.begin());
}

IntParamsOrigin::IntParamsOrigin(StateIn& si) : IntParams(si) {
  si.get(O_);
}

IntParamsOrigin::~IntParamsOrigin() {}

void
IntParamsOrigin::save_data_state(StateOut& so)
{
  IntParams::save_data_state(so);
  so.put(O_);
}

const double*
IntParamsOrigin::r() const {
  return &O_[0];
}

double
IntParamsOrigin::r(unsigned int xyz) const {
  return O_.at(xyz);
}

bool
IntParamsOrigin::equiv(const IntParams& other) const {
  return (downcast<IntParamsOrigin>(other) != 0);
}

/////////////////////////////

static ClassDesc IntParamsG12_cd(
  typeid(IntParamsG12),"IntParamsG12",1,"public IntParams",
  0, 0, create<IntParamsG12>);

IntParamsG12::IntParamsG12(const ContractedGeminal& bra) :
  IntParams(2), bra_(bra), ket_(null_geminal)
{
  if (bra_.size() == 0)
    throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal contractions of zero length",__FILE__,__LINE__);

  typedef ContractedGeminal::const_iterator citer;
  citer end = bra_.end();
  for(citer i=bra_.begin(); i<end; i++) {
    if ( (*i).first < 0.0)
      throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
}

IntParamsG12::IntParamsG12(const ContractedGeminal& bra,
                           const ContractedGeminal& ket) :
  IntParams(2), bra_(bra), ket_(ket)
{
  if (bra_.size() == 0 ||
      ket_.size() == 0)
    throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal contractions of zero length",__FILE__,__LINE__);

  typedef ContractedGeminal::const_iterator citer;
  citer end = bra_.end();
  for(citer i=bra_.begin(); i<end; i++) {
    if ( (*i).first < 0.0)
      throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
  end = ket_.end();
  for(citer i=ket_.begin(); i<end; i++) {
    if ( (*i).first < 0.0)
      throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
}

IntParamsG12::IntParamsG12(StateIn& si) : SavableState(si)
{
  si.get(bra_);
  si.get(ket_);
}

IntParamsG12::~IntParamsG12() {}

void
IntParamsG12::save_data_state(StateOut& so)
{
  IntParams::save_data_state(so);
  so.put(bra_);
  so.put(ket_);
}

const IntParamsG12::ContractedGeminal&
IntParamsG12::bra() const { return bra_; }

const IntParamsG12::ContractedGeminal&
IntParamsG12::ket() const { return ket_; }

IntParamsG12::ContractedGeminal
IntParamsG12::zero_exponent_geminal(1,std::make_pair(0.0,1.0));

double
IntParamsG12::null_exponent(DBL_MAX);

IntParamsG12::ContractedGeminal
IntParamsG12::null_geminal(1,std::make_pair(null_exponent,1.0));

IntParamsG12::PrimitiveGeminal
IntParamsG12::product(const PrimitiveGeminal& A,
                      const PrimitiveGeminal& B) {
  return std::make_pair(A.first+B.first,A.second*B.second);
}

IntParamsG12::ContractedGeminal
IntParamsG12::product(const ContractedGeminal& A,
                      const ContractedGeminal& B) {
  const unsigned int na = A.size();
  const unsigned int nb = B.size();
  ContractedGeminal result;
  for(unsigned int a=0; a<na; ++a) {
    for(unsigned int b=0; b<nb; ++b) {
      result.push_back( product(A[a],B[b]) );
    }
  }
  return result;
}

bool
IntParamsG12::equiv(const IntParams& other) const {
  const IntParamsG12* otherptr = downcast<IntParamsG12>(other);
  if (otherptr) { // compare contents
    return (otherptr->bra_ == bra_ &&
            otherptr->ket_ == ket_);
  }
  return false;
}
