//
// moints_runtime.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_mointsruntime_h
#define _mpqc_src_lib_chemistry_qc_lcao_mointsruntime_h

#include <chemistry/qc/lcao/tbint_runtime.h>
#include <chemistry/qc/lcao/df_runtime.h>

namespace sc {

  /** MOIntsRuntime provides runtime support for computing 1-body and 2-body (2-, 3-, and 4-center) MO-basis integrals
      (with or without density fitting).
    */
  class MOIntsRuntime : virtual public SavableState {
    public:
      typedef MOIntsRuntime this_type;

      /// give density fitting basis to enable density fitting, when possible
      MOIntsRuntime(const Ref<MOIntsTransformFactory>& factory,
                    const Ref<DensityFittingParams>& dfparams = 0);
      ~MOIntsRuntime();
      MOIntsRuntime(StateIn&);
      void save_data_state(StateOut&);

      /// obsoletes this object
      void obsolete();

      /// factory for creating AO->MO transforms for 2-body integrals
      const Ref<MOIntsTransformFactory>& factory() const { return factory_; }
      /// density fitting parameters
      const Ref<DensityFittingParams>& dfparams() const { return dfparams_; }
      /// runtime for density fitting matrices. Returns null if density fitting basis was not given.
      const Ref<DensityFittingRuntime>& runtime_df() const { return runtime_df_; }
      /// runtime for 2-body 2-center integrals
      const Ref<TwoBodyTwoCenterMOIntsRuntime>& runtime_2c() const { return runtime_2c_; }
      /// runtime for 2-body 3-center integrals
      const Ref<TwoBodyThreeCenterMOIntsRuntime>& runtime_3c() const { return runtime_3c_; }
      /// runtime for 2-body 4-center integrals
      const Ref<TwoBodyFourCenterMOIntsRuntime>& runtime_4c() const { return runtime_4c_; }

    private:
      static ClassDesc class_desc_;

      Ref<MOIntsTransformFactory> factory_;
      Ref<DensityFittingParams> dfparams_;
      Ref<DensityFittingInfo> dfinfo_;
      Ref<DensityFittingRuntime> runtime_df_;
      Ref<TwoBodyTwoCenterMOIntsRuntime> runtime_2c_;
      Ref<TwoBodyThreeCenterMOIntsRuntime> runtime_3c_;
      Ref<TwoBodyFourCenterMOIntsRuntime> runtime_4c_;
  };

  //////////////////////

} // end of namespace sc

#endif /* end of header guard */

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
