//
// transform_iRjS.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_transformiRjS_h
#define _mpqc_src_lib_chemistry_qc_lcao_transformiRjS_h

#include <chemistry/qc/lcao/transform_tbint.h>

namespace sc {

  /** TwoBodyMOIntsTransform_iRjS computes (iR|jS), or <ij|RS> integrals,
      where R and S are atomic orbitals,
      using parallel integral-direct AO->MO transformation.
      */
  class TwoBodyMOIntsTransform_iRjS: public TwoBodyMOIntsTransform {

      // Initialize the target integrals accumulator
      void init_acc();
      // Compute required dynamic memory for a given batch size
      distsize_t compute_transform_dynamic_memory_(int ni) const;

    public:

      TwoBodyMOIntsTransform_iRjS(StateIn&);
      TwoBodyMOIntsTransform_iRjS(const std::string& name,
                                  const Ref<MOIntsTransformFactory>& factory,
                                  const Ref<TwoBodyIntDescr>& tbintdescr,
                                  const Ref<OrbitalSpace>& space1,
                                  const Ref<OrbitalSpace>& space2,
                                  const Ref<OrbitalSpace>& space3,
                                  const Ref<OrbitalSpace>& space4);
      ~TwoBodyMOIntsTransform_iRjS();

      void save_data_state(StateOut&);

      /// Implementation of TwoBodyMOIntsTransform::type()
      std::string type() const {
        return "iRjS";
      }

      /** Returns the number of bytes allocated for each ij-block of integrals of one type
       in MemoryGrp */
      size_t memgrp_blksize() const;

      /// Computes transformed integrals
      void compute();
      /// Check symmetry of transformed integrals
      void check_int_symm(double threshold =
          TwoBodyMOIntsTransform::zero_integral);

    private:
      Ref<GaussianBasisSet> basis2_;
      Ref<GaussianBasisSet> basis4_;
  };

} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
