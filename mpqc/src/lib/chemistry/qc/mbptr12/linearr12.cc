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
#include <util/class/scexception.h>
#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/linearr12.h>

using namespace sc;
using namespace LinearR12;

CorrelationFactor::CorrelationFactor(const std::string& label, const CorrelationParameters& params) :
  label_(label), params_(params)
{
}

CorrelationFactor::~CorrelationFactor()
{
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
    
    // It doesn't make sense to print out parameters for some correlation factors
    Ref<G12CorrelationFactor> g12ptr;  g12ptr << const_cast<CorrelationFactor*>(this);
    if (g12ptr.nonnull()) {
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

int
LinearR12::CorrelationFactor::tbint_type_eri() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_t2f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t2f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_f12eri() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_f12f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
LinearR12::CorrelationFactor::tbint_type_f12t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbintdescr(f) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

Ref<TwoBodyIntDescr>
CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int fbra, unsigned int fket) const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbintdescr(f,g) -- should not be called for this CorrelationFactor",__FILE__,__LINE__);
}

////

LinearR12::NullCorrelationFactor::NullCorrelationFactor() :
  CorrelationFactor("NONE")
{
}

int
LinearR12::NullCorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

double
LinearR12::NullCorrelationFactor::value(unsigned int c, double r12) const
{
  return 0;
}

////

LinearR12::R12CorrelationFactor::R12CorrelationFactor() :
  CorrelationFactor("R12", CorrelationParameters(1,ContractedGeminal(1,std::make_pair(0.0,1.0))))
{
}

int
LinearR12::R12CorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
LinearR12::R12CorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12);
}

int
LinearR12::R12CorrelationFactor::tbint_type_t1f12() const
{
  return static_cast<int>(TwoBodyInt::r12t1);
}

int
LinearR12::R12CorrelationFactor::tbint_type_t2f12() const
{
  return static_cast<int>(TwoBodyInt::r12t2);
}

Ref<TwoBodyIntDescr>
LinearR12::R12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  return new TwoBodyIntDescrR12(IF);
}

double
LinearR12::R12CorrelationFactor::value(unsigned int c, double r12) const
{
  return r12;
}

////

LinearR12::G12CorrelationFactor::G12CorrelationFactor(const CorrelationParameters& params) :
  CorrelationFactor("G12",params)
{
}

int
LinearR12::G12CorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
LinearR12::G12CorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
LinearR12::G12CorrelationFactor::tbint_type_t1f12() const
{
  return static_cast<int>(TwoBodyInt::t1g12);
}

int
LinearR12::G12CorrelationFactor::tbint_type_t2f12() const
{
  return static_cast<int>(TwoBodyInt::t2g12);
}

int
LinearR12::G12CorrelationFactor::tbint_type_f12eri() const
{
  return static_cast<int>(TwoBodyInt::r12_m1_g12);
}

int
LinearR12::G12CorrelationFactor::tbint_type_f12f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
LinearR12::G12CorrelationFactor::tbint_type_f12t1f12() const
{
  return static_cast<int>(TwoBodyInt::g12t1g12);
}

Ref<TwoBodyIntDescr>
LinearR12::G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f), IntParamsG12::zero_exponent_geminal);
  return new TwoBodyIntDescrG12(IF,params);
}

Ref<TwoBodyIntDescr>
LinearR12::G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
                                            unsigned int fbra,
                                            unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12(IF,params);
}

double
LinearR12::G12CorrelationFactor::value(unsigned int c, double r12) const
{
  double val = 0.0;
  const unsigned int nprims = nprimitives(c);
  for(unsigned int p=0; p<nprims; p++) {
    const PrimitiveGeminal& prim = primitive(c,p);
    const double expon = prim.first;
    const double coef = prim.second;
    val += coef*exp(-expon*r12*r12);
  }
  return val;
}

