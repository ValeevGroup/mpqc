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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_transformfactorytimpl_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_transformfactorytimpl_h

#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/transform_tbint.h>

namespace sc {

  template <typename TransformType> Ref<TwoBodyMOIntsTransform>
    MOIntsTransformFactory::twobody_transform(const std::string& name,
                                              const Ref<TwoBodyIntDescr>& descrarg)
  {
    Ref<TwoBodyMOIntsTransform> result;
    const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);
    result = new TransformType(name,this,descr,space1_,space2_,space3_,space4_);

    if (top_mole_.nonnull())
      result->set_top_mole(top_mole_);
    result->set_debug(debug());
    reserve_memory(result->memory());

    return result;
  }

  // code goes here

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
