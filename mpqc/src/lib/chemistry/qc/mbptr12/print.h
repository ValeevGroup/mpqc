//
// print.h
//
// Copyright (C) 2006 Edward Valeev
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_print_h
#define _chemistry_qc_mbptr12_print_h

namespace sc {

  /// Default print thresholds
  class DefaultPrintThresholds {
  public:
    /// Only essential results + events
    static const unsigned int terse = 0;
    /// Additional diagnostics, e.g., transform progress, etc.
    static const unsigned int diagnostics = 1;
    /// Additional fine-level diagnostics, e.g. detailed timings, etc.
    static const unsigned int fine = 2;
    /// Print N^0 quantities
    static const unsigned int N0 = terse;
    /// Print most N^0 quantities
    static const unsigned int mostN0 = diagnostics;
    /// Print most N quantities
    static const unsigned int mostN = 3;
    /// Print essential N^2 quantities
    static const unsigned int N2 = 2;
    /// Print most N^2 quantities
    static const unsigned int mostN2 = 3;
    /// Print all N^2 quantities
    static const unsigned int allN2 = 5;
    /// Print essential o^4 quantities
    static const unsigned int O4 = 3;
    /// Print most o^4 quantities
    static const unsigned int mostO4 = 4;
    /// Print all o^4 quantities
    static const unsigned int allO4 = 6;
    /// Print essential o^2N^2 quantities
    static const unsigned int O2N2 = 5;
    /// Print most o^2N^2 quantities
    static const unsigned int mostO2N2 = 6;
    /// Print all o^2N^2 quantities
    static const unsigned int allO2N2 = 7;
    /// Print essential N^4 quantities, e.g., fully transformed integrals of any kind
    static const unsigned int N4 = 6;
    /// Print most N^4 quantities, e.g., partially transformed integrals
    static const unsigned int mostN4 = 7;
    /// Print all N^4 quantities, e.g., partially transformed integrals
    static const unsigned int allN4 = 10;
    /// Print essential o^6 quantities
    static const unsigned int O6 = N4;
    /// Print most o^6 quantities
    static const unsigned int mostO6 = mostN4;
  };

};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
