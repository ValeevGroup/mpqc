//
// intparams.cc
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

#include <float.h>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/intparams.h>

using namespace sc;

/////////////////////////////

IntParams::IntParams(unsigned int nparams) :
  nparams_(nparams)
{}

IntParams::~IntParams()
{}

unsigned int
IntParams::nparams() const
{
  return nparams_;
}

/////////////////////////////

IntParamsVoid::IntParamsVoid() : IntParams(0) {}

IntParamsVoid::~IntParamsVoid() {}

/////////////////////////////

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

IntParamsG12::~IntParamsG12()
{}

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

/////////////////////////////

IntParamsGenG12::IntParamsGenG12(const ContractedGeminal& bra) :
  IntParams(2), bra_(bra), ket_(null_geminal)
{
  if (bra_.size() == 0)
    throw ProgrammingError("IntParamsGenG12::IntParamsGenG12() -- geminal contractions of zero length",__FILE__,__LINE__);

#if 0  
  typedef ContractedGeminal::const_iterator citer;
  citer end = bra_.end();
  for(citer i=bra_.begin(); i<end; i++) {
    if ( (*i).second < 0.0)
      throw ProgrammingError("IntParamsGenG12::IntParamsGenG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
#endif
}

IntParamsGenG12::IntParamsGenG12(const ContractedGeminal& bra,
                           const ContractedGeminal& ket) :
  IntParams(2), bra_(bra), ket_(ket)
{
  if (bra != ket)
    throw ProgrammingError("IntParamsGenG12::IntParamsGenG12(bra,ket) -- bra != ket",__FILE__,__LINE__);

  if (bra_.size() == 0)
    throw ProgrammingError("IntParamsGenG12::IntParamsGenG12() -- geminal contractions of zero length",__FILE__,__LINE__);

#if 0  
  typedef ContractedGeminal::const_iterator citer;
  citer end = bra_.end();
  for(citer i=bra_.begin(); i<end; i++) {
    if ( (*i).first < 0.0)
      throw ProgrammingError("IntParamsGenG12::IntParamsGenG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
#endif
}

IntParamsGenG12::~IntParamsGenG12()
{}

const IntParamsGenG12::ContractedGeminal&
IntParamsGenG12::bra() const { return bra_; }

const IntParamsGenG12::ContractedGeminal&
IntParamsGenG12::ket() const { return ket_; }

IntParamsGenG12::ContractedGeminal
IntParamsGenG12::zero_exponent_geminal(1,std::make_pair(std::make_pair(0.0,0.0),1.0));

double
IntParamsGenG12::null_exponent(DBL_MAX);

IntParamsGenG12::ContractedGeminal
IntParamsGenG12::null_geminal(1,std::make_pair(std::make_pair(null_exponent,null_exponent),1.0));
