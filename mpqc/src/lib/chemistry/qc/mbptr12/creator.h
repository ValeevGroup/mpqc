//
// creator.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_creator_h
#define _chemistry_qc_mbptr12_creator_h

#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/moints_runtime.h>

namespace sc {

  /** RangeCreator<T> is Functor which can be used up to n times to create objects
   of type T. operator() returns the objects, or 0 when done. Thus T must be Comparable to
   an int or Constructable from an int. */
  template<typename T>
  class RangeCreator {
    public:
      RangeCreator(unsigned int n) :
        n_(n), ncreated_(0) {
      }
      /// returns a new object T, or null() if done. \sa null()
      virtual T operator()() =0;
      /// returns the null object. Default is to return T(0).
      virtual T null() const {
        return T(0);
      }

    protected:
      bool can_create() const {
        return ncreated_ < n_;
      }
      void next() {
        ++ncreated_;
      }
      unsigned int ncreated() const {
        return ncreated_;
      }

    private:
      unsigned int n_;
      unsigned int ncreated_;
  };

  /** Creates new DistArray4 using TwoBodyFourCenterMOIntsRuntime and a vector of transform keys */
  class DistArray4Creator: public RangeCreator<Ref<DistArray4> > {
    public:
      typedef Ref<DistArray4> ObjT;

      DistArray4Creator(const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_rtime,
                        const std::vector<std::string>& tform_keys);
      /// Implementation of RangeCreator::operator()
      ObjT operator()();

    private:
      const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_rtime_;
      const std::vector<std::string>& tform_keys_;
  };

  /** Creates TwoBodyIntDescr for correlation factor C */
  class TwoBodyIntDescrCreator: public RangeCreator<Ref<TwoBodyIntDescr> > {
    public:
      typedef Ref<TwoBodyIntDescr> ObjT;

      TwoBodyIntDescrCreator(const Ref<LinearR12::CorrelationFactor>& corrfactor,
                             const Ref<Integral>& integral,
                             bool CorrFunctionInBra = false,
                             bool CorrFunctionInKet = false);
      /// Implementation of RangeCreator::operator()
      ObjT operator()();

    private:
      Ref<LinearR12::CorrelationFactor> corrfactor_;
      Ref<Integral> integral_;
      bool CorrFunctionInBraKet_;
      unsigned int nf12bra_;
      unsigned int nf12ket_;
      unsigned int braindex_;
      unsigned int ketindex_;

      void increment_indices();
  };

  /** Creates R12TwoBodyIntKey for the given CorrelationFactor */
  class R12TwoBodyIntKeyCreator: public RangeCreator<std::string> {
    public:
      typedef std::string ObjT;

      R12TwoBodyIntKeyCreator(const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_rtime_,
                              const Ref<OrbitalSpace>& bra1,
                              const Ref<OrbitalSpace>& ket1,
                              const Ref<OrbitalSpace>& bra2,
                              const Ref<OrbitalSpace>& ket2,
                              const Ref<LinearR12::CorrelationFactor>& corrfactor,
                              bool CorrFunctionInBra = false,
                              bool CorrFunctionInKet = false,
                              std::string layout_key = std::string(TwoBodyIntLayout::b1b2_k1k2));
      /// Implementation of RangeCreator::operator()
      ObjT operator()();

      ObjT null() const;

    private:
      Ref<TwoBodyFourCenterMOIntsRuntime> moints_rtime_;
      Ref<LinearR12::CorrelationFactor> corrfactor_;
      const Ref<OrbitalSpace>& bra1_;
      const Ref<OrbitalSpace>& bra2_;
      const Ref<OrbitalSpace>& ket1_;
      const Ref<OrbitalSpace>& ket2_;
      bool CorrFunctionInBra_;
      bool CorrFunctionInKet_;
      std::string layout_key_;
      unsigned int nf12bra_;
      unsigned int nf12ket_;
      unsigned int braindex_;
      unsigned int ketindex_;

      void increment_indices();
  };

}

#endif

