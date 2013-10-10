//
// spin.h
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

#include <stdexcept>
#include <util/misc/scexception.h>

#ifndef _chemistry_qc_mbptr12_spin_h
#define _chemistry_qc_mbptr12_spin_h

namespace sc {

  typedef enum { NSpinCases1 = 2, NSpinCases2 = 3} NSpinCases;
  typedef enum { NPureSpinCases2 = 2 } NPureSpinCases;
  /// Returns the number of unique spin cases (1 or 2)
  unsigned int nspincases1(bool spin_polarized);
  /// Returns the number of unique combinations of 2 spin cases (1 or 3)
  unsigned int nspincases2(bool spin_polarized);
  /// Returns the number of pure 2 spin cases
  unsigned int npurespincases2();
  enum { ToLowerCase = true,
     ToUpperCase = false };

  /// @addtogroup ChemistryElectronicStructureSpin
  /// @{

  typedef enum { AnySpinCase1 = -1, Alpha = 0, Beta = 1, InvalidSpinCase1 = 2} SpinCase1;
  typedef enum { AnySpinCase2 = -1, AlphaBeta = 0, AlphaAlpha = 1, BetaBeta = 2, InvalidSpinCase2 = 3} SpinCase2;
  typedef enum { AnyPureSpinCase2 = -1, Singlet = 0, Triplet = 1, InvalidPureSpinCase2 = 2} PureSpinCase2;
  /// returns the first spin case of the 2-spin S
  SpinCase1 case1(SpinCase2 S);
  /// returns the second spin case of the 2-spin S
  SpinCase1 case2(SpinCase2 S);
  /// combines 2 spins to give 1 2-spin
  SpinCase2 case12(SpinCase1 S1, SpinCase1 S2);
  /// given 1-spin return the other 1-spin
  SpinCase1 other(SpinCase1 S);

  std::string to_string(SpinCase1 S);
  std::string to_string(SpinCase2 S);
  std::string to_string(PureSpinCase2 S);
  SpinCase1 to_spincase1(const std::string& key);
  SpinCase2 to_spincase2(const std::string& key);
  PureSpinCase2 to_purespincase2(const std::string& key);

  /// Prepend string representation of S to R and return
  std::string prepend_spincase(SpinCase1 S, const std::string& R, bool lowercase = false);
  /// Prepend string representation of S to R and return
  std::string prepend_spincase(SpinCase2 S, const std::string& R, bool lowercase = false);
  /// Prepend string representation of S to R and return
  std::string prepend_spincase(PureSpinCase2 S, const std::string& R, bool lowercase = false);

  /// @}
  // end of addtogroup ChemistryElectronicStructure

};

#endif

