//
// cart.h
//
// Copyright (C) 2011 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_basis_cart_h
#define _mpqc_src_lib_chemistry_qc_basis_cart_h

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>

namespace sc {

  /// CartesianBasisSet is obtained from the parent basis by converting spherical harmonic shells to cartesian
  /// counterparts
  class CartesianBasisSet : public GaussianBasisSet {
    public:
      /** A KeyVal constructor is used to generate a CartesianBasisSet
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>basis</tt><td>GaussianBasisSet<td>none<td>parent basis
          <tr><td><tt>integral</tt><td>Integral<td>defined by environment<td>Integral factory

          </table>
       */
      CartesianBasisSet(const Ref<KeyVal>& kv);
      CartesianBasisSet(const Ref<GaussianBasisSet>& basis,
                        Ref<Integral> integral = Integral::get_default_integral());
      CartesianBasisSet(StateIn&);
      virtual ~CartesianBasisSet();
      void save_data_state(StateOut&);

      const Ref<GaussianBasisSet>& parent() const;

    private:
      static ClassDesc class_desc_;

      Ref<Integral> integral_;  //< this defines the transformation from cartesian to spherical
      Ref<GaussianBasisSet> parent_;  //< parent basis

      void convert(const Ref<GaussianBasisSet>& parent,
                   const Ref<Integral>& integral);
  };


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
