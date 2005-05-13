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

Ref<CorrelationFactor>
CorrelationFactor::Instance(CorrelationFactorID id)
{
  Ref<CorrelationFactor> cf;
  switch (id) {
  case R12CorrFactor:
    cf = new CorrelationFactor(id,"R12",&Integral::grt,4);
    cf->tbint_type_eri(TwoBodyInt::eri);
    cf->tbint_type_f12(TwoBodyInt::r12);
    cf->tbint_type_t1f12(TwoBodyInt::r12t1);
    cf->tbint_type_t2f12(TwoBodyInt::r12t2);
    break;
  case G12CorrFactor:
    cf = new CorrelationFactor(id,"G12",&Integral::g12,6);
    cf->tbint_type_eri(TwoBodyInt::eri);
    cf->tbint_type_f12(TwoBodyInt::r12_0_g12);
    cf->tbint_type_t1f12(TwoBodyInt::t1g12);
    cf->tbint_type_t2f12(TwoBodyInt::t2g12);
    cf->tbint_type_f12eri(TwoBodyInt::r12_m1_g12);
    cf->tbint_type_f12f12(TwoBodyInt::r12_0_g12);
    cf->tbint_type_f12t1f12(TwoBodyInt::g12t1g12);
  default:
    throw FeatureNotImplemented("CorrelationFactor::Instance() -- correlation factor with this id not found",__FILE__,__LINE__);
  }
  return cf;
}

CorrelationFactor::CorrelationFactor(CorrelationFactorID id, const std::string& label,
                                     IntegralCallback callback, unsigned int num_tbint_types) :
  id_(id), label_(label), callback_(callback), num_tbint_types_(num_tbint_types),
  tbint_type_eri_(-1),
  tbint_type_f12_(-1),
  tbint_type_t1f12_(-1),
  tbint_type_t2f12_(-1),
  tbint_type_f12eri_(-1),
  tbint_type_f12f12_(-1),
  tbint_type_f12t1f12_(-1)
{
}

CorrelationFactor::~CorrelationFactor()
{
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

