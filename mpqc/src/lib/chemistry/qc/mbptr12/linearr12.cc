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

CorrelationFactor::CorrelationFactor(const std::string& label) :
  label_(label)
{
}

CorrelationFactor::~CorrelationFactor()
{
}

unsigned int
CorrelationFactor::nfunctions() const
{
  return 1;
}

unsigned int
CorrelationFactor::nprimitives(unsigned int c) const
{
  return 1;
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
    print_params(os,f);
    os << decindent;
  }
  os << decindent;
}

void
CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
}

double
CorrelationFactor::value(unsigned int c, double r12, double r1, double r2) const
{
  return value(c,r12);
}

int
CorrelationFactor::tbint_type_eri() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_t2f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_t2f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_f12eri() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12eri() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_f12f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_f12t1f12() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12t1f12() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
}

int
CorrelationFactor::tbint_type_f12f12_anti() const
{
  throw ProgrammingError("LinearR12::CorrelationFactor::tbint_type_f12f12_anti() -- invalid type of integrals for the given CorrelationFactor",__FILE__,__LINE__);
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

NullCorrelationFactor::NullCorrelationFactor() :
  CorrelationFactor("NONE")
{
}

int
NullCorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

double
NullCorrelationFactor::value(unsigned int c, double r12) const
{
  return 0.0;
}

////

R12CorrelationFactor::R12CorrelationFactor() :
  CorrelationFactor("R12")
{
}

int
R12CorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
R12CorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12);
}

int
R12CorrelationFactor::tbint_type_t1f12() const
{
  return static_cast<int>(TwoBodyInt::r12t1);
}

int
R12CorrelationFactor::tbint_type_t2f12() const
{
  return static_cast<int>(TwoBodyInt::r12t2);
}

Ref<TwoBodyIntDescr>
R12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  return new TwoBodyIntDescrR12(IF);
}

double
R12CorrelationFactor::value(unsigned int c, double r12) const
{
  return r12;
}

////

G12CorrelationFactor::G12CorrelationFactor(const CorrelationParameters& params) :
  CorrelationFactor("G12"), params_(params)
{
}

unsigned int
G12CorrelationFactor::nfunctions() const
{
  return params_.size();
}

const G12CorrelationFactor::ContractedGeminal&
G12CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
G12CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const G12CorrelationFactor::PrimitiveGeminal&
G12CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}


int
G12CorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
G12CorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
G12CorrelationFactor::tbint_type_t1f12() const
{
  return static_cast<int>(TwoBodyInt::t1g12);
}

int
G12CorrelationFactor::tbint_type_t2f12() const
{
  return static_cast<int>(TwoBodyInt::t2g12);
}

int
G12CorrelationFactor::tbint_type_f12eri() const
{
  return static_cast<int>(TwoBodyInt::r12_m1_g12);
}

int
G12CorrelationFactor::tbint_type_f12f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
G12CorrelationFactor::tbint_type_f12t1f12() const
{
  return static_cast<int>(TwoBodyInt::g12t1g12);
}

Ref<TwoBodyIntDescr>
G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12(IF,params);
}

Ref<TwoBodyIntDescr>
G12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
                                            unsigned int fbra,
                                            unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12(IF,params);
}

double
G12CorrelationFactor::value(unsigned int c, double r12) const
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

void
G12CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[" << prim.first << " " << prim.second << "] ";
  }
  os << " ]" << endl;
}


////

G12NCCorrelationFactor::G12NCCorrelationFactor(const CorrelationParameters& params) :
  CorrelationFactor("G12"), params_(params)
{
}

unsigned int
G12NCCorrelationFactor::nfunctions() const
{
  return params_.size();
}

const G12NCCorrelationFactor::ContractedGeminal&
G12NCCorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
G12NCCorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const G12NCCorrelationFactor::PrimitiveGeminal&
G12NCCorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}


int
G12NCCorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
G12NCCorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
G12NCCorrelationFactor::tbint_type_f12eri() const
{
  return static_cast<int>(TwoBodyInt::r12_m1_g12);
}

int
G12NCCorrelationFactor::tbint_type_f12f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_g12);
}

int
G12NCCorrelationFactor::tbint_type_f12t1f12() const
{
  return static_cast<int>(TwoBodyInt::g12t1g12);
}

int
G12NCCorrelationFactor::tbint_type_f12f12_anti() const
{
  return static_cast<int>(TwoBodyInt::anti_g12g12);
}

Ref<TwoBodyIntDescr>
G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(f));
  return new TwoBodyIntDescrG12NC(IF,params);
}

Ref<TwoBodyIntDescr>
G12NCCorrelationFactor::tbintdescr(const Ref<Integral>& IF,
					      unsigned int fbra,
					      unsigned int fket) const
{
  Ref<IntParamsG12> params = new IntParamsG12(function(fbra),function(fket));
  return new TwoBodyIntDescrG12NC(IF,params);
}

double
G12NCCorrelationFactor::value(unsigned int c, double r12) const
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

void
G12NCCorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ Exponent Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[" << prim.first << " " << prim.second << "] ";
  }
  os << " ]" << endl;
}


////

GenG12CorrelationFactor::GenG12CorrelationFactor(const CorrelationParameters& params) :
  CorrelationFactor("GenG12"), params_(params)
{
  if (params_.size() != 1)
    throw ProgrammingError("GenG12CorrelationFactor::GenG12CorrelationFactor() -- only 1 general Geminal correlation factor can now be handled",__FILE__,__LINE__);
}

unsigned int
GenG12CorrelationFactor::nfunctions() const
{
  return params_.size();
}


const GenG12CorrelationFactor::ContractedGeminal&
GenG12CorrelationFactor::function(unsigned int c) const
{
  return params_.at(c);
}

unsigned int
GenG12CorrelationFactor::nprimitives(unsigned int c) const
{
  return params_.at(c).size();
}

const GenG12CorrelationFactor::PrimitiveGeminal&
GenG12CorrelationFactor::primitive(unsigned int c, unsigned int p) const
{
  return params_.at(c).at(p);
}

int
GenG12CorrelationFactor::tbint_type_eri() const
{
  return static_cast<int>(TwoBodyInt::eri);
}

int
GenG12CorrelationFactor::tbint_type_f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_gg12);
}

int
GenG12CorrelationFactor::tbint_type_f12eri() const
{
  return static_cast<int>(TwoBodyInt::r12_m1_gg12);
}

int
GenG12CorrelationFactor::tbint_type_f12f12() const
{
  return static_cast<int>(TwoBodyInt::r12_0_gg12);
}

int
GenG12CorrelationFactor::tbint_type_f12t1f12() const
{
  return static_cast<int>(TwoBodyInt::gg12t1gg12);
}

Ref<TwoBodyIntDescr>
GenG12CorrelationFactor::tbintdescr(const Ref<Integral>& IF, unsigned int f) const
{
  Ref<IntParamsGenG12> params = new IntParamsGenG12(function(f));
  return new TwoBodyIntDescrGenG12(IF,params);
}

Ref<TwoBodyIntDescr>
GenG12CorrelationFactor::tbintdescr(const Ref<Integral>& IF,
					       unsigned int fbra,
					       unsigned int fket) const
{
  Ref<IntParamsGenG12> params = new IntParamsGenG12(function(fbra),function(fket));
  return new TwoBodyIntDescrGenG12(IF,params);
}

double
GenG12CorrelationFactor::value(unsigned int c, double r12) const
{
  throw ProgrammingError("GenG12CorrelationFactor::value(c,r12) -- not defined for general Geminal correlation factors",__FILE__,__LINE__);
}

double
GenG12CorrelationFactor::value(unsigned int c, double r12, double r1, double r2) const
{
  double val = 0.0;
  const unsigned int nprims = nprimitives(c);
  for(unsigned int p=0; p<nprims; p++) {
    const PrimitiveGeminal& prim = primitive(c,p);
    const double alpha = prim.first.first;
    const double gamma = prim.first.second;
    const double coef = prim.second;
    val += coef*exp( - alpha*(r1*r1 + r2*r2) - gamma*r12*r12 );
  }
  return val;
}

void
GenG12CorrelationFactor::print_params(std::ostream& os, unsigned int f) const
{
  using std::endl;

  os << indent << "[ [Alpha Gamma] Coefficient] = [ ";
  const int nprim = nprimitives(f);
  for(int p=0; p<nprim; p++) {
    const PrimitiveGeminal& prim = primitive(f,p);
    os << "[ [" << prim.first.first << " " << prim.first.second << "] " << prim.second << "] ";
  }
  os << " ]" << endl;
}
