//
// spin.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/spin.h>
#include <util/class/scexception.h>
#include <ctype.h>

namespace sc {
  unsigned int nspincases1(bool spin_polarized) { return spin_polarized ? 2 : 1; }
  unsigned int nspincases2(bool spin_polarized) { return spin_polarized ? 3 : 1; }
  unsigned int npurespincases2(){ return(NPureSpinCases2); }
  SpinCase1 case1(SpinCase2 S) { return S==BetaBeta ? Beta : Alpha; }
  SpinCase1 case2(SpinCase2 S) { return S==AlphaAlpha ? Alpha : Beta; }
  std::string to_string(SpinCase1 S) {
    switch(S) {
      case Alpha:
        return std::string("alpha");
      case Beta:
        return std::string("beta");
      case AnySpinCase1:
        return std::string("");
    }
    throw ProgrammingError("to_string() -- invalid argument",__FILE__,__LINE__);
  }

  std::string to_string(SpinCase2 S) {
    switch(S) {
      case AlphaAlpha:
        return std::string("alpha-alpha");
      case AlphaBeta:
        return std::string("alpha-beta");
      case BetaBeta:
        return std::string("beta-beta");
      case AnySpinCase2:
        return std::string("");
    }
    throw ProgrammingError("to_string() -- invalid argument",__FILE__,__LINE__);
  }

  std::string to_string(PureSpinCase2 S) {
    switch(S) {
      case Singlet:
        return std::string("singlet");
      case Triplet:
        return std::string("triplet");
      case AnyPureSpinCase2:
        return std::string("");
    }
    throw ProgrammingError("to_string() -- invalid argument",__FILE__,__LINE__);
  }

  SpinCase1 to_spincase1(const std::string& key) {
    for(int s=AnySpinCase1; s!=InvalidSpinCase1; ++s) {
      SpinCase1 sc = static_cast<SpinCase1>(s);
      if (key == to_string(sc))
        return sc;
    }
    throw ProgrammingError("to_spincase1() -- invalid argument",__FILE__,__LINE__);
  }

  PureSpinCase2 to_purespincase2(const std::string& key) {
    for(int s=AnyPureSpinCase2; s!=InvalidPureSpinCase2; ++s) {
      PureSpinCase2 sc = static_cast<PureSpinCase2>(s);
      if (key == to_string(sc))
        return sc;
    }
    throw ProgrammingError("to_purespincase2() -- invalid argument",__FILE__,__LINE__);
  }

  SpinCase2 to_spincase2(const std::string& key) {
    for(int s=AnySpinCase2; s!=InvalidSpinCase2; ++s) {
      SpinCase2 sc = static_cast<SpinCase2>(s);
      if (key == to_string(sc))
        return sc;
    }
    throw ProgrammingError("to_spincase2() -- invalid argument",__FILE__,__LINE__);
  }

  std::string prepend_spincase(SpinCase2 S, const std::string& R, bool lowercase)
  {
    std::string prefix(to_string(S));
    if (!lowercase)
      prefix[0] = toupper(prefix[0]);
    return prefix + ' ' + R;
  }

  std::string prepend_spincase(PureSpinCase2 S, const std::string& R, bool lowercase)
  {
    std::string prefix(to_string(S));
    if (!lowercase)
      prefix[0] = toupper(prefix[0]);
    return prefix + ' ' + R;
  }

  std::string prepend_spincase(SpinCase1 S, const std::string& R, bool lowercase)
  {
    std::string prefix(to_string(S));
    if (!lowercase)
      prefix[0] = toupper(prefix[0]);
    return prefix + ' ' + R;
  }

}

