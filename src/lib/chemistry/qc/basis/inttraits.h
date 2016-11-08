//
// inttraits.h
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

#ifndef _mpqc_src_lib_chemistry_qc_basis_inttraits_h
#define _mpqc_src_lib_chemistry_qc_basis_inttraits_h

#include <cassert>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intparams.h>

namespace sc {

  // auxiliary type functions
  namespace detail {

    template <> struct SEvalCreator<2u> {
      static Ref< OneBodyIntEvalType<2u>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsVoid>& params) {
        return factory->overlap();
      }
    };
  };

  /// Traits of a set of one-body integrals
  template <int NumCenters, OneBodyOperSet::type Type> struct OneBodyIntTraits {
    /// the type of the NBodyInt object that evaluates this set
    typedef typename OneBodyIntEvalType<NumCenters>::value EvalType;
    /** the type of IntParams object needed to initialize the evaluator
        for computing this set of integrals */
    typedef typename OneBodyIntParamsType<Type>::value ParamsType;
    /// the type of the NBodyInt object that evaluates this set
    typedef OneBodyOperSetProperties<Type> TypeMap;
    /// number of integral types
    static const int size = TypeMap::size;
    /// creates an Eval object
    static Ref<EvalType> eval(const Ref<Integral>& factory,
                              const Ref<ParamsType>& params) {
      typedef typename detail::OneBodyEvalCreator<NumCenters,Type>::value OneBodyEvalCreator;
      return OneBodyEvalCreator::eval(factory,params);
    }
    /// maps index of the integral type within this set to TwoBodyOper::type
    static OneBodyOper::type intset(unsigned int t) {
      MPQC_ASSERT(t < size);
      return TypeMap::value[t];
    }
    /// inverse of the above intset
    static unsigned int intset(OneBodyOper::type t) {
      for(unsigned int i=0; i<size; ++i)
        if (TypeMap::value[i] == t)
          return i;
      abort();   // should be unreachable if input is valid
    }

  };

  // auxiliary type functions
  namespace {

    template <int NumCenters> struct ERIEvalCreator {
      static Ref< typename TwoBodyIntEvalType<NumCenters>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsVoid>& params) {
        return factory->coulomb<NumCenters>();
      }
    };

    template <int NumCenters> struct R12EvalCreator {
      static Ref< typename TwoBodyIntEvalType<NumCenters>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsVoid>& params) {
        return factory->grt<NumCenters>();
      }
    };

    template <int NumCenters> struct G12EvalCreator {
      static Ref< typename TwoBodyIntEvalType<NumCenters>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsG12>& params) {
        return factory->g12<NumCenters>(params);
      }
    };

    template <int NumCenters> struct G12NCEvalCreator {
      static Ref< typename TwoBodyIntEvalType<NumCenters>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsG12>& params) {
        return factory->g12nc<NumCenters>(params);
      }
    };

    template <int NumCenters> struct G12DKHEvalCreator {
      static Ref< typename TwoBodyIntEvalType<NumCenters>::value >
      eval(const Ref<Integral>& factory, const Ref<IntParamsG12>& params) {
        return factory->g12dkh<NumCenters>(params);
      }
    };

    template <int NumCenters, TwoBodyOperSet::type Type> struct TwoBodyEvalCreator;
    template <int NumCenters> struct TwoBodyEvalCreator<NumCenters,TwoBodyOperSet::ERI> {
      typedef ERIEvalCreator<NumCenters> value;
    };
    template <int NumCenters> struct TwoBodyEvalCreator<NumCenters,TwoBodyOperSet::R12> {
      typedef R12EvalCreator<NumCenters> value;
    };
    template <int NumCenters> struct TwoBodyEvalCreator<NumCenters,TwoBodyOperSet::G12> {
      typedef G12EvalCreator<NumCenters> value;
    };
    template <int NumCenters> struct TwoBodyEvalCreator<NumCenters,TwoBodyOperSet::G12NC> {
      typedef G12NCEvalCreator<NumCenters> value;
    };
    template <int NumCenters> struct TwoBodyEvalCreator<NumCenters,TwoBodyOperSet::G12DKH> {
      typedef G12DKHEvalCreator<NumCenters> value;
    };
  };

  /// Traits of a set of two-body integrals
  template <int NumCenters, TwoBodyOperSet::type Type> struct TwoBodyIntTraits {
    /// the type of the NBodyInt object that evaluates this set
    typedef typename TwoBodyIntEvalType<NumCenters>::value EvalType;
    /** the type of IntParams object needed to initialize the evaluator
        for computing this set of integrals */
    typedef typename TwoBodyIntParamsType<Type>::value ParamsType;
    /// the type of the NBodyInt object that evaluates this set
    typedef TwoBodyOperSetProperties<Type> TypeMap;
    /// number of integral types
    static const int size = TypeMap::size;
    /// creates an Eval object
    static Ref<EvalType> eval(const Ref<Integral>& factory,
                              const Ref<ParamsType>& params) {
      typedef typename detail::TwoBodyEvalCreator<NumCenters,Type>::value TwoBodyEvalCreator;
      return TwoBodyEvalCreator::eval(factory,params);
    }
    /// maps index of the integral type within this set to TwoBodyOper::type
    static TwoBodyOper::type intset(unsigned int t) {
      MPQC_ASSERT(t < size);
      return TypeMap::value[t];
    }
    /// inverse of the above intset
    static unsigned int intset(TwoBodyOper::type t) {
      for(unsigned int i=0; i<size; ++i)
        if (TypeMap::value[i] == t)
          return i;
      abort();   // should be unreachable if input is valid
    }

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
