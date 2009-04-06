//
// union.h
//
// Copyright (C) 2009 Edward Valeev
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_basis_union_h
#define _mpqc_src_lib_chemistry_qc_basis_union_h

#include <chemistry/qc/basis/basis.h>

namespace sc {

  /** The UnionBasisSet class is a union of 2 GaussianBasisSet objects
  */
  class UnionBasisSet: public GaussianBasisSet {

    public:

      /** The KeyVal constructor.  It does not call the basis class
          KeyVal constructor.

          <dl>

          <dt><tt>basis1</tt><dd> The first GaussianBasisSet object in the union.
          <dt><tt>basis2</tt><dd> The second GaussianBasisSet object in the union.

          </dl>

          */
      UnionBasisSet(const Ref<KeyVal>&);
      UnionBasisSet(StateIn&);

      void save_data_state(StateOut&);

    private:

      // set to 1 to debug
      static int debug() { return 0; }

      static ClassDesc class_desc_;

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
