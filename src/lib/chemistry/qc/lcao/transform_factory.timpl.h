//
// transform_factory.timpl.h
//
// Copyright (C) 2008 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_transformfactorytimpl_h
#define _mpqc_src_lib_chemistry_qc_lcao_transformfactorytimpl_h

#include <chemistry/qc/lcao/transform_factory.h>
#include <chemistry/qc/lcao/transform_tbint.h>

namespace sc {

  namespace detail {

    // these classes are needed for type-indirection
    template <bool DF> struct MakeTwoBodyTransform;
    template <> struct MakeTwoBodyTransform<false> {
      static const int i = 1;
      template <typename TransformType> static
      Ref<TwoBodyMOIntsTransform>
      evaluate(MOIntsTransformFactory* factory,
          const std::string& name,
          const Ref<TwoBodyIntDescr>& descrarg) {
        Ref<TwoBodyMOIntsTransform> result;
        const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? factory->tbintdescr() : descrarg);
        result = new TransformType(name,factory,descr,
                                   factory->space1_,
                                   factory->space2_,
                                   factory->space3_,
                                   factory->space4_);

        if (factory->top_mole_)
          result->set_top_mole(factory->top_mole_);
        result->set_debug(factory->debug());

        return result;
      }
    };
    template <> struct MakeTwoBodyTransform<true> {
      static const int i = 1;
      template <typename TransformType> static
      Ref<TwoBodyMOIntsTransform>
      evaluate(MOIntsTransformFactory* factory,
          const std::string& name,
          const Ref<TwoBodyIntDescr>& descrarg) {
        Ref<TwoBodyMOIntsTransform> result;
        const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? factory->tbintdescr() : descrarg);
        result = new TransformType(name,factory->df_info(),descr,
                                   factory->space1_,
                                   factory->space2_,
                                   factory->space3_,
                                   factory->space4_);

        if (factory->top_mole_)
          result->set_top_mole(factory->top_mole_);
        result->set_debug(factory->debug());

        return result;
      }
    };

    template <typename TransformType> struct NeedDF {
      static const bool value = false;
    };
    template <> struct NeedDF<TwoBodyMOIntsTransform_ixjy_df> {
      static const bool value = true;
    };

    // convert a transform type to the DF-based transform. Default is to not use a DF-based transform
    template <typename TransformType> struct ToDensityFittingType {
      typedef TransformType value;
    };
    template <> struct ToDensityFittingType<TwoBodyMOIntsTransform_iRjS> {
      typedef TwoBodyMOIntsTransform_ixjy_df value;
    };
#if 1
    template <> struct ToDensityFittingType<TwoBodyMOIntsTransform_ixjy> {
      typedef TwoBodyMOIntsTransform_ixjy_df value;
    };
    template <> struct ToDensityFittingType<TwoBodyMOIntsTransform_ikjy> {
      typedef TwoBodyMOIntsTransform_ixjy_df value;
    };
#endif

    // compare types
    template <typename A, typename B> struct EqualTypes {
      static const bool value = false;
    };
    template <typename A> struct EqualTypes<A,A> {
      static const bool value = true;
    };
  }

  template <typename TransformType> Ref<TwoBodyMOIntsTransform>
    MOIntsTransformFactory::twobody_transform(const std::string& name,
                                              const Ref<TwoBodyIntDescr>& descrarg)
  {
    // is density fitting is possible, try calling for a DF-based transform (unless already asked for a DF-based transform)
    if (df_info() != 0 &&
        ! detail::EqualTypes<TransformType, TwoBodyMOIntsTransform_ixjy_df>::value &&
        ! detail::EqualTypes<TransformType, typename detail::ToDensityFittingType<TransformType>::value>::value) {
      return this->twobody_transform< typename detail::ToDensityFittingType<TransformType>::value >(name, descrarg);
    }

    typedef detail::MakeTwoBodyTransform< detail::NeedDF<TransformType>::value > TformMaker;
    Ref<TwoBodyMOIntsTransform> result = TformMaker::template evaluate<TransformType>(this,name,descrarg);
    return result;
  }

  template <typename TransformType> Ref<TwoBodyThreeCenterMOIntsTransform>
    MOIntsTransformFactory::twobody_transform(const std::string& name,
                      const Ref<TwoBodyThreeCenterIntDescr>& descrarg) {
    Ref<TwoBodyThreeCenterMOIntsTransform> result = new TransformType(name,this,descrarg,space1_,space2_,space3_);
#if 0
    if (top_mole_)
      result->set_top_mole(top_mole_);
#endif
    return result;
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
