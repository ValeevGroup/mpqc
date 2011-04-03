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

#ifndef _mpqc_src_lib_chemistry_qc_basis_union_h
#define _mpqc_src_lib_chemistry_qc_basis_union_h

#include <chemistry/qc/basis/basis.h>

namespace sc {

  /** UnionBasisSet constructs a union of two GaussianBasisSet objects.
      The union basis  on atom i includes the basis
      functions of A centered on i followed by the basis functions of B
      centered on i. Duplicate basis functions/shells are eliminated.
      The Molecule object for the two
      basis sets must be identical.
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
      UnionBasisSet(const Ref<GaussianBasisSet>& bs1,
                    const Ref<GaussianBasisSet>& bs2);

      void save_data_state(StateOut&);

      /// return basis1
      const Ref<GaussianBasisSet>& basis1() const;
      /// return basis2
      const Ref<GaussianBasisSet>& basis2() const;

      enum Basis12 { Basis1=1, Basis2=2, Basis1_and_Basis2=12 };
      /// reports in which basis  shell s of the union basis is found
      Basis12 shell_to_basis(int s) const;
      /// reports in which basis  function f of the union basis is found
      Basis12 function_to_basis(int f) const;

    private:

      Ref<GaussianBasisSet> basis1_;
      Ref<GaussianBasisSet> basis2_;

      // maps shell to the basis from which it came from
      std::vector<int> shell_to_basis_;
      // maps function to the basis from which it came from
      std::vector<int> function_to_basis_;

      /// computes A + B
      void sum(const Ref<GaussianBasisSet>& A,
               const Ref<GaussianBasisSet>& B);

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
