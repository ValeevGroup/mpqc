//
// spin.h
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

#ifdef __GNUC__
#pragma interface
#endif

#include <stdexcept>
#include <util/misc/scexception.h>

#ifndef _chemistry_qc_mbptr12_spin_h
#define _chemistry_qc_mbptr12_spin_h

namespace sc {
  
  typedef enum { NSpinCases1 = 2, NSpinCases2 = 3} NSpinCases;
  typedef enum { Alpha = 0, Beta = 1} SpinCase1;
  typedef enum { AlphaBeta = 0, AlphaAlpha = 1, BetaBeta = 2} SpinCase2;
  /// Returns the number of unique combinations of 2 spin cases
  int nspincases2(bool spin_polarized);
  /// returns the first spin case of the 2-spin S
  SpinCase1 case1(SpinCase2 S);
  /// returns the second spin case of the 2-spin S
  SpinCase1 case2(SpinCase2 S);
  
};

#endif

