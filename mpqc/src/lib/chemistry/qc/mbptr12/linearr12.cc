//
// linearr12.cc
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/linearr12.h>

using namespace sc;
using namespace LinearR12;

CorrelationFactor::CorrelationFactor(const CorrelationFactorID& id, const CorrelationParameters& params) :
  id_(id), params_(params)
{
  init();
}

void
CorrelationFactor::init()
{
  switch (id_) {
  case NullCorrFactor:
    init2("NONE",&Integral::grt,1);
    tbint_type_eri(TwoBodyInt::eri);
    break;
  case R12CorrFactor:
    init2("R12",&Integral::grt,4);
    tbint_type_eri(TwoBodyInt::eri);
    tbint_type_f12(TwoBodyInt::r12);
    tbint_type_t1f12(TwoBodyInt::r12t1);
    tbint_type_t2f12(TwoBodyInt::r12t2);
    break;
  case G12CorrFactor:
    init2("G12",&Integral::g12,6);
    tbint_type_eri(TwoBodyInt::eri);
    tbint_type_f12(TwoBodyInt::r12_0_g12);
    tbint_type_t1f12(TwoBodyInt::t1g12);
    tbint_type_t2f12(TwoBodyInt::t2g12);
    tbint_type_f12eri(TwoBodyInt::r12_m1_g12);
    tbint_type_f12f12(TwoBodyInt::r12_0_g12);
    tbint_type_f12t1f12(TwoBodyInt::g12t1g12);
    break;
  default:
    throw FeatureNotImplemented("CorrelationFactor::CorrelationFactor() -- correlation factor with this id not found",__FILE__,__LINE__);
  }
}

void
CorrelationFactor::init2(const std::string& label,
                         IntegralCallback callback,
                         unsigned int num_tbint_types)
{
  label_ = label;
  callback_ = callback;
  num_tbint_types_ = num_tbint_types;
  tbint_type_eri_ = -1;
  tbint_type_f12_ = -1;
  tbint_type_t1f12_ = -1;
  tbint_type_t2f12_ = -1;
  tbint_type_f12eri_ = -1;
  tbint_type_f12f12_ = -1;
  tbint_type_f12t1f12_ = -1;
}

CorrelationFactor::~CorrelationFactor()
{
}

const LinearR12::CorrelationFactorID&
CorrelationFactor::id() const
{
  return id_;
}

unsigned int
CorrelationFactor::nfunctions() const
{
  return params_.size();
}

const LinearR12::CorrelationFactor::ContractedGeminal&
CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const LinearR12::CorrelationFactor::PrimitiveGeminal&
CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}

const std::string&
CorrelationFactor::label() const
{
  return label_;
}

const CorrelationFactor::IntegralCallback&
CorrelationFactor::callback() const
{
  return callback_;
}

unsigned int
CorrelationFactor::num_tbint_types() const
{
  return num_tbint_types_;
}

void
CorrelationFactor::print(std::ostream& os) const
{
  using std::endl;
  os << indent << "CorrelationFactor:" << endl;
  os << incindent;
  const int nfunc = nfunctions();
  for(int f=0; f<nfunc; f++) {
    os << indent << "Function " << f << ":" << endl << incindent;
    os << indent << "Functional form: " << label() << endl;
    if (id() == G12CorrFactor) {
      os << indent << "[ Exponent Coefficient] = [ ";
      const int nprim = nprimitives(f);
      for(int p=0; p<nprim; p++) {
        const PrimitiveGeminal& prim = primitive(f,p);
        os << "[" << prim.first << " " << prim.second << "] ";
      }
      os << " ]" << endl;
    }
    os << decindent;
  }
  os << decindent;
}

