//
// spin.cpp
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

#include "mpqc/chemistry/qc/wfn/spin.h"

#include <cassert>
#include <ctype.h>

#include "mpqc/util/misc/assert.h"
#include "mpqc/util/misc/exception.h"

namespace mpqc {
constexpr unsigned int nspincases1(bool spin_polarized) {
  return spin_polarized ? 2 : 1;
}
constexpr unsigned int nspincases2(bool spin_polarized) {
  return spin_polarized ? 3 : 1;
}
constexpr unsigned int npurespincases2() { return (NPureSpinCases2); }
constexpr SpinCase1 case1(SpinCase2 S) {
  return S == SpinCase2::BetaBeta ? SpinCase1::Beta : SpinCase1::Alpha;
}
constexpr SpinCase1 case2(SpinCase2 S) {
  return S == SpinCase2::AlphaAlpha ? SpinCase1::Alpha : SpinCase1::Beta;
}
constexpr SpinCase2 case12(SpinCase1 S1, SpinCase1 S2) {
  return (S1 == SpinCase1::Alpha && S2 == SpinCase1::Alpha)
             ? SpinCase2::AlphaAlpha
             : ((S1 == SpinCase1::Alpha && S2 == SpinCase1::Beta)
                    ? SpinCase2::AlphaBeta
                    : ((S1 == SpinCase1::Beta && S2 == SpinCase1::Alpha)
                           ? SpinCase2::AlphaBeta
                           : ((S1 == SpinCase1::Beta && S2 == SpinCase1::Beta)
                                  ? SpinCase2::BetaBeta
                                  : SpinCase2::Invalid)));
}
constexpr SpinCase1 other(SpinCase1 S) {
  return (S == SpinCase1::Alpha)
             ? SpinCase1::Beta
             : ((S == SpinCase1::Beta) ? SpinCase1::Alpha : SpinCase1::Invalid);
}
constexpr const char* to_string(SpinCase1 S) {
  return (S == SpinCase1::Alpha)
             ? "alpha"
             : ((S == SpinCase1::Beta)
                    ? "beta"
                    : (S == SpinCase1::Any ? "" : "invalid"));
}

constexpr const char* to_string(SpinCase2 S) {
  return (S == SpinCase2::AlphaAlpha)
             ? "alpha-alpha"
             : ((S == SpinCase2::AlphaBeta)
                    ? "alpha-beta"
                    : ((S == SpinCase2::BetaBeta)
                           ? "beta-beta"
                           : ((S == SpinCase2::Any) ? "" : "invalid")));
}

constexpr const char* to_string(PureSpinCase2 S) {
  return (S == PureSpinCase2::Singlet)
             ? "singlet"
             : ((S == PureSpinCase2::Triplet)
                    ? "triplet"
                    : ((S == PureSpinCase2::Any) ? "" : "invalid"));
}

SpinCase1 to_spincase1(std::string key) {
  for (int s = static_cast<int>(SpinCase1::Any); s != static_cast<int>(SpinCase1::Invalid); ++s) {
    SpinCase1 sc = static_cast<SpinCase1>(s);
    if (key == to_string(sc)) return sc;
  }
  throw ProgrammingError("to_spincase1() -- invalid argument", __FILE__,
                         __LINE__);
}

PureSpinCase2 to_purespincase2(std::string key) {
  for (int s = static_cast<int>(PureSpinCase2::Any);
       s != static_cast<int>(PureSpinCase2::Invalid); ++s) {
    PureSpinCase2 sc = static_cast<PureSpinCase2>(s);
    if (key == to_string(sc)) return sc;
  }
  throw ProgrammingError("to_purespincase2() -- invalid argument", __FILE__,
                         __LINE__);
}

SpinCase2 to_spincase2(std::string key) {
  for (int s = static_cast<int>(SpinCase2::Any);
       s != static_cast<int>(SpinCase2::Invalid); ++s) {
    SpinCase2 sc = static_cast<SpinCase2>(s);
    if (key == to_string(sc)) return sc;
  }
  throw ProgrammingError("to_spincase2() -- invalid argument", __FILE__,
                         __LINE__);
}

std::string prepend_spincase(SpinCase2 S, const std::string& R,
                             bool lowercase) {
  std::string prefix(to_string(S));
  if (!lowercase) prefix[0] = toupper(prefix[0]);
  return prefix + ' ' + R;
}

std::string prepend_spincase(PureSpinCase2 S, const std::string& R,
                             bool lowercase) {
  std::string prefix(to_string(S));
  if (!lowercase) prefix[0] = toupper(prefix[0]);
  return prefix + ' ' + R;
}

std::string prepend_spincase(SpinCase1 S, const std::string& R,
                             bool lowercase) {
  std::string prefix(to_string(S));
  if (!lowercase) prefix[0] = toupper(prefix[0]);
  return prefix + ' ' + R;
}
}
