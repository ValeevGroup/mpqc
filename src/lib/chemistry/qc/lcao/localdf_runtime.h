//
// localdf_runtime.h
//
// Copyright (C) 2013 MPQC Developers
//
// Author: David Hollman <dhollman@vt.edu>
// Maintainer: DSH
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

#ifndef _mpqc_src_lib_chemistry_qc_df_localdf_runtime_h
#define _mpqc_src_lib_chemistry_qc_df_localdf_runtime_h

#include <chemistry/qc/lcao/df.h>
#include <chemistry/qc/lcao/tbint_runtime.h>

namespace sc {

  /** Parsed representation of a string key that represents fitting of a product of space1 and space2 into fspace
      Coulomb fitting kernel_key is the default. */
  class ParsedLocalDensityFittingKey : ParsedDensityFittingKey {
    public:
      ParsedLocalDensityFittingKey(const std::string& key)
        : ParsedDensityFittingKey(key) {
      };

  };

}

class LocalDensityFittingRuntime : virtual public SavableState {
  public:
    typedef LocalDensityFittingRuntime this_type;

    // TODO Finish moving stuff over to this new class


};

#endif
