//
// intdescr.h
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

#ifndef _chemistry_qc_basis_intdescr_h
#define _chemistry_qc_basis_intdescr_h

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/intparams.h>
#include <chemistry/qc/basis/inttraits.h>

namespace sc {

  template <typename IntEval> struct IntEvalToOperSetType;
  template <> struct IntEvalToOperSetType<TwoBodyInt> {
    typedef TwoBodyOperSet::type value;
  };
  template <> struct IntEvalToOperSetType<TwoBodyThreeCenterInt> {
    typedef TwoBodyOperSet::type value;
  };
  template <> struct IntEvalToOperSetType<TwoBodyTwoCenterInt> {
    typedef TwoBodyOperSet::type value;
  };

  /** IntegralSetDescr contains all information necessary to construct an IntEval object that computes
      a particular set of integrals using an Integral factory. This means that it must provide all parameters
      necessary to construct such set (i.e. geminal exponents, etc.).
      It also allows to programmatically (at runtime) to examine the set (i.e. ask what
      is the first integral type?).

      @tparam IntEval the base evaluator type (e.g. TwoBodyInt) for all integrals in the integral set. */
  template <typename IntEval>
  class IntegralSetDescr : public RefCount {
    public:
    IntegralSetDescr() {}

    /// the factory used to create integral evaluator
    virtual const Ref<Integral>& factory() const =0;
    /// call appropriate method to produce TwoBodyInt corresponding to this Set
    virtual Ref<IntEval> inteval() const =0;
    /// how many integral sets
    virtual unsigned int num_sets() const =0;
    /// optional parameters that determine the operator set (e.g., geminal exponents, etc.)
    virtual Ref<IntParams> params() const =0;
    /// the type of the operator set
    virtual typename IntEvalToOperSetType<IntEval>::value operset() const =0;

    /// Maps integral set t to its index in this set
    virtual unsigned int intset(TwoBodyOper::type t) const =0;
    /// Maps integral set t to its TwoBodyOper::type
    virtual TwoBodyOper::type intset(unsigned int t) const =0;

    /// iterator over Property classes for each Integral in the set (Property describes symmetry, labels, etc.)
    //PropertyIterator propiter();
  };

  template <int NumCenters> struct TwoBodyIntType;
  template <> struct TwoBodyIntType<4> {
    typedef TwoBodyInt value;
  };
  template <> struct TwoBodyIntType<3> {
    typedef TwoBodyThreeCenterInt value;
  };
  template <> struct TwoBodyIntType<2> {
    typedef TwoBodyTwoCenterInt value;
  };

  /// Implements descriptors for various two-body evaluators
  template <int NumCenters, TwoBodyOperSet::type TwoBodyIntSet>
    class TwoBodyNCenterIntDescr : public IntegralSetDescr< typename TwoBodyIntType<NumCenters>::value > {
      public:
        typedef TwoBodyIntTraits<NumCenters,TwoBodyIntSet> TraitsType;
        typedef typename TraitsType::ParamsType ParamsType;
        typedef typename TwoBodyIntType<NumCenters>::value EvalType;

        static const unsigned int num_intsets = TraitsType::size;
        TwoBodyNCenterIntDescr(const Ref<Integral>& IF,
                               const Ref<ParamsType>& params = Ref<ParamsType>(dynamic_cast<ParamsType*>(new IntParamsVoid))) :
                                 factory_(IF),
                                 params_(params) { }
        TwoBodyNCenterIntDescr(const Ref<Integral>& IF,
                               const Ref<IntParams>& params) :
                                 factory_(IF),
                                 params_(0) {
          params_ << params;
          MPQC_ASSERT(params_);
        }
        ~TwoBodyNCenterIntDescr() {}

        /// which factory is used
        const Ref<Integral>& factory() const { return factory_; }
        // implementation of TwoBodyIntDescr::inteval()
        Ref<EvalType> inteval() const {
          return TraitsType::eval(factory_, params_);
        }
        // implementation of TwoBodyIntDescr::params()
        Ref<IntParams> params() const {
          return params_;
        }
        // implementation of TwoBodyIntDescr::operset()
        TwoBodyOperSet::type operset() const { return TwoBodyIntSet; }
        // implementation of TwoBodyIntDescr::num_sets()
        unsigned int num_sets() const { return num_intsets; }
        // Implementation of TwoBodyIntDescr::intset()
        unsigned int intset(TwoBodyOper::type t) const {
          return intSet(t);
        }
        // Implementation of TwoBodyIntDescr::intset()
        TwoBodyOper::type intset(unsigned int t) const {
          return intSet(t);
        }
        /// Static version of TwoBodyIntDescr::intset()
        static unsigned int intSet(TwoBodyOper::type t) {
          return TraitsType::intset(t);
        }
        /// Static version of TwoBodyIntDescr::intset()
        static TwoBodyOper::type intSet(unsigned int t) {
          return TraitsType::intset(t);
        }

      private:
        /// which factory is used
        Ref<Integral> factory_;
        /// the parameters of the operator set
        Ref<ParamsType> params_;

    };


  typedef IntegralSetDescr<TwoBodyInt> TwoBodyIntDescr;
  typedef IntegralSetDescr<TwoBodyThreeCenterInt> TwoBodyThreeCenterIntDescr;
  typedef IntegralSetDescr<TwoBodyTwoCenterInt> TwoBodyTwoCenterIntDescr;
  template <int NumCenters, int NumParticles> struct NCentersToDescr;
  template <> struct NCentersToDescr<4,2> {
    typedef TwoBodyIntDescr value;
  };
  template <> struct NCentersToDescr<3,2> {
    typedef TwoBodyThreeCenterIntDescr value;
  };
  template <> struct NCentersToDescr<2,2> {
    typedef TwoBodyTwoCenterIntDescr value;
  };

  typedef TwoBodyNCenterIntDescr<4,TwoBodyOperSet::ERI> TwoBodyIntDescrERI;
  typedef TwoBodyNCenterIntDescr<3,TwoBodyOperSet::ERI> TwoBodyThreeCenterIntDescrERI;
  typedef TwoBodyNCenterIntDescr<2,TwoBodyOperSet::ERI> TwoBodyTwoCenterIntDescrERI;
  typedef TwoBodyNCenterIntDescr<4,TwoBodyOperSet::R12> TwoBodyIntDescrR12;
  typedef TwoBodyNCenterIntDescr<3,TwoBodyOperSet::R12> TwoBodyThreeCenterIntDescrR12;
  typedef TwoBodyNCenterIntDescr<2,TwoBodyOperSet::R12> TwoBodyTwoCenterIntDescrR12;
  typedef TwoBodyNCenterIntDescr<4,TwoBodyOperSet::G12> TwoBodyIntDescrG12;
  typedef TwoBodyNCenterIntDescr<3,TwoBodyOperSet::G12> TwoBodyThreeCenterIntDescrG12;
  typedef TwoBodyNCenterIntDescr<2,TwoBodyOperSet::G12> TwoBodyTwoCenterIntDescrG12;
  typedef TwoBodyNCenterIntDescr<4,TwoBodyOperSet::G12NC> TwoBodyIntDescrG12NC;
  typedef TwoBodyNCenterIntDescr<3,TwoBodyOperSet::G12NC> TwoBodyThreeCenterIntDescrG12NC;
  typedef TwoBodyNCenterIntDescr<2,TwoBodyOperSet::G12NC> TwoBodyTwoCenterIntDescrG12NC;
  typedef TwoBodyNCenterIntDescr<4,TwoBodyOperSet::G12DKH> TwoBodyIntDescrG12DKH;
  typedef TwoBodyNCenterIntDescr<3,TwoBodyOperSet::G12DKH> TwoBodyThreeCenterIntDescrG12DKH;
  typedef TwoBodyNCenterIntDescr<2,TwoBodyOperSet::G12DKH> TwoBodyTwoCenterIntDescrG12DKH;

  struct IntDescrFactory {
    template <int NumCenters>
      static Ref< typename NCentersToDescr<NumCenters,2>::value >
      make(const Ref<Integral>& integral,
           TwoBodyOperSet::type type,
           const Ref<IntParams>& params) {
        typedef typename NCentersToDescr<NumCenters,2>::value ReturnType;
        switch (type) {
          case TwoBodyOperSet::ERI: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::ERI> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::R12: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::R12> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::G12: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::G12> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::G12NC: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::G12NC> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::G12DKH: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::G12DKH> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::R12_0_G12: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::R12_0_G12> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::R12_m1_G12: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::R12_m1_G12> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::G12_T1_G12: {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::G12_T1_G12> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          case TwoBodyOperSet::DeltaFunction : {
            typedef TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::DeltaFunction> ConcreteType;
            return new ConcreteType(integral,params);
          } break;
          default:
            MPQC_ASSERT(false);
        }
        return Ref< typename NCentersToDescr<NumCenters,2>::value >(); // dummy return statement to pacify picky compilers
      }
  };

}

#endif

