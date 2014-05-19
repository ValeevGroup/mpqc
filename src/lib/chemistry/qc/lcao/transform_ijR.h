//
// transform_ijR.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_transform_ijR_h
#define _mpqc_src_lib_chemistry_qc_lcao_transform_ijR_h

#include <chemistry/qc/basis/intdescr.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/lcao/transform_factory.h>

namespace sc {

  /** TwoBodyThreeCenterMOIntsTransform_ijR computes (ij|R) integrals,
      where R are atomic orbitals, using parallel integral-direct AO->MO transformation.
      */
  class TwoBodyThreeCenterMOIntsTransform_ijR: public TwoBodyThreeCenterMOIntsTransform {
    public:
      TwoBodyThreeCenterMOIntsTransform_ijR(StateIn&);
      TwoBodyThreeCenterMOIntsTransform_ijR(const std::string& name,
                                 const Ref<MOIntsTransformFactory>& factory,
                                 const Ref<TwoBodyThreeCenterIntDescr>& tbintdescr,
                                 const Ref<OrbitalSpace>& space1,
                                 const Ref<OrbitalSpace>& space2,
                                 const Ref<OrbitalSpace>& space3);
      ~TwoBodyThreeCenterMOIntsTransform_ijR();
      void save_data_state(StateOut&);

      // implements TwoBodyThreeCenterMOIntsTransform::type();
      std::string type() const {
        return "ijR";
      }

      /// Computes transformed integrals
      void compute();

    private:

      static ClassDesc class_desc_;

      int compute_transform_batchsize(size_t mem_static, int rank_R);
      // implements TwoBodyThreeCenterMOIntsTransform::compute_transform_dynamic_memory()
      distsize_t compute_transform_dynamic_memory(int batchsize = -1) const;
      // implements TwoBodyThreeCenterMOIntsTransform::init_acc()
      void init_acc();
      // implements TwoBodyThreeCenterMOIntsTransform::extra_memory_report()
      void extra_memory_report(std::ostream& os = ExEnv::out0()) const;

      /// specialization of compute() for the case where space1 and space2 is AO
      void compute_pqR();
      /// specialization of compute() for the case where space1, but not space2, is AO
      void compute_pjR();
      /// specialization of compute() for the case where neither space1 nor space2 is AO
      void compute_ijR();

  };

  /** TwoBodyThreeCenterMOIntsTransform_ijR computes (ij|R) integrals,
      where R are atomic orbitals, using (iq|R) integrals
      */
  class TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR: public TwoBodyThreeCenterMOIntsTransform_ijR {
    public:
      TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR(StateIn&);
      TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR(const std::string& name,
                                                      const Ref<TwoBodyThreeCenterMOIntsTransform_ijR>& iqR_tform,
                                                      const Ref<OrbitalSpace>& space2);
      ~TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR();
      void save_data_state(StateOut&);

      // implements TwoBodyThreeCenterMOIntsTransform::type();
      std::string type() const {
        return "ijR";
      }

      /// Computes transformed integrals
      void compute();

    private:

      static ClassDesc class_desc_;

      Ref<TwoBodyThreeCenterMOIntsTransform_ijR> iqR_tform_;

      int compute_transform_batchsize(size_t mem_static, int rank_R);
      // reimplements TwoBodyThreeCenterMOIntsTransform::compute_transform_dynamic_memory()
      distsize_t compute_transform_dynamic_memory(int batchsize = -1) const;
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
