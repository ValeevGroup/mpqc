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

#include <util/misc/scexception.h>
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

IntParamsG12::IntParamsG12(double g1, double g2) :
  IntParams(2), gamma1_(g1), gamma2_(g2)
{
  if (gamma1_ < 0.0 || gamma2_ < 0.0) {
    throw ProgrammingError("IntParamsG12::IntParamsG12() -- geminal parameters must be nonnegative",__FILE__,__LINE__);
  }
}

IntParamsG12::~IntParamsG12()
{}

double
IntParamsG12::gamma1() const { return gamma1_; }

double
IntParamsG12::gamma2() const { return gamma2_; }
