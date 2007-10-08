//
// pairiter.impl.h
//
// Copyright (C) 2007 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_pairiterimpl_h
#define _chemistry_qc_mbptr12_pairiterimpl_h

#include <assert.h>
#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/pairiter.h>

namespace sc {
  namespace fastpairiter {
    
    template <>
    MOPairIter<Symm>::MOPairIter(int nI, int nJ) :
      nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0)
      {
        assert(nI == nJ);
        nIJ_ = nI_ * (nI_ + 1)/2;
        init();
      }

    template <>
    MOPairIter<AntiSymm>::MOPairIter(int nI, int nJ) :
      nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0)
      {
        assert(nI == nJ);
        nIJ_ = nI_ * (nI_ - 1)/2;
      }

    template <>
    MOPairIter<ASymm>::MOPairIter(int nI, int nJ) :
      nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0)
      {
        nIJ_ = nI_ * nJ_;
      }

    template <PairSymm PSymm>
    MOPairIter<PSymm>::~MOPairIter() {
    }

    template <>
    void MOPairIter<AntiSymm>::init() {
      I_ = 1; J_ = 0; IJ_ = 0;
    }
    
    template <PairSymm PSymm>
    void MOPairIter<PSymm>::init() {
      I_ = 0; J_ = 0; IJ_ = 0;
    }

    template <PairSymm PSymm>
    void MOPairIter<PSymm>::start() {
      init();
    }
    
    template <>
    void MOPairIter<AntiSymm>::next() {
      ++J_;
      if (J_ == I_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }

    template <>
    void MOPairIter<ASymm>::next() {
      ++J_;
      if (J_ == nJ_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }

    template <>
    void MOPairIter<Symm>::next() {
      ++J_;
      if (J_ > I_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }

    template <>
    int MOPairIter<AntiSymm>::ij(int i, int j) const {
      return i*(i-1)/2 + j;
    }

    template <>
    int MOPairIter<Symm>::ij(int i, int j) const {
      return i*(i+1)/2 + j;
    }

    template <>
    int MOPairIter<ASymm>::ij(int i, int j) const {
      return i*nJ_ + j;
    }

    template <PairSymm PSymm>
    MOPairIter<PSymm>::operator int() const {
      return (IJ_ < nIJ_);
    }

    template <>
    RefSCDimension pairdim<AntiSymm>(int nI, int nJ) {
      assert(nI == nJ);
      return new SCDimension(nI * (nI-1)/2);
    }

    template <>
    RefSCDimension pairdim<ASymm>(int nI, int nJ) {
      return new SCDimension(nI * nJ);
    }

    template <>
    RefSCDimension pairdim<Symm>(int nI, int nJ) {
      assert(nI == nJ);
      return new SCDimension(nI * (nI+1)/2);
    }

  } // namespace fastpairiter
}

#endif

