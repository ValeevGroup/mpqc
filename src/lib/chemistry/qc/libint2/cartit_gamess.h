//
// cartit_gamess.h
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

#ifndef _chemistry_qc_libint2_cartitgamess_h
#define _chemistry_qc_libint2_cartitgamess_h

#include <chemistry/qc/basis/cartiter.h>
#include <cgshellinfo.h>
#include <chemistry/qc/libint2/macros_gamess.h>

namespace sc {

  class CartesianIterGAMESS: public CartesianIter {

      static const int lmax = LIBINT2_CARTGAUSS_MAX_AM;
      static ::libint2::CGShellInfo< ::libint2::CGShellOrderingData< ::libint2::CGShellOrdering_GAMESS,lmax> > ordering_data_;

    public:
      CartesianIterGAMESS(int l) :
        CartesianIter(l) {
      }

      void start() {
        bfn_ = 0;
        ordering_data_.cartindex_to_ijk(l_, bfn_, a_, b_, c_);
      }

      void next() {
        ++bfn_;
        if (bfn_ >= INT_NCART(l_)) // if iterated over all basis functions
          a_ = -1; // set a_ to an invalid value
        else
          ordering_data_.cartindex_to_ijk(l_, bfn_, a_, b_, c_);
      }

      operator int() {
        return (a_ >= 0); // check that a_ is valid
      }
  };

  class RedundantCartesianIterGAMESS: public RedundantCartesianIter {

      static const int lmax = LIBINT2_CARTGAUSS_MAX_AM;
      static ::libint2::CGShellInfo< ::libint2::CGShellOrderingData< ::libint2::CGShellOrdering_GAMESS,lmax> > ordering_data_;

    public:
      RedundantCartesianIterGAMESS(int l) :
        RedundantCartesianIter(l) {
      }

      int bfn() {
        return ordering_data_.cartindex(l(), a(), b());
      }
  };

  class RedundantCartesianSubIterGAMESS: public RedundantCartesianSubIter {

      static const int lmax = LIBINT2_CARTGAUSS_MAX_AM;
      static ::libint2::CGShellInfo< ::libint2::CGShellOrderingData< ::libint2::CGShellOrdering_GAMESS,lmax> > ordering_data_;

      int bfn_;

    public:
      RedundantCartesianSubIterGAMESS(int l) :
        RedundantCartesianSubIter(l) {
      }

      void start(int aa, int bb, int cc) {
        RedundantCartesianSubIter::start(aa, bb, cc);
        bfn_ = ordering_data_.cartindex(l(), a(), b());
      }

      int bfn() const { return bfn_; }
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
