//
// transform_ixjy.h
//
// Copyright (C) 2004 Edward Valeev
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

#ifndef _chemistry_qc_lcao_transformixjy_h
#define _chemistry_qc_lcao_transformixjy_h

#include <string>
#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/lcao/transform_tbint.h>

namespace sc {

  /** TwoBodyMOIntsTransform_ixjy computes (ix|jy) integrals
      using parallel integrals-direct AO->MO transformation. */

class TwoBodyMOIntsTransform_ixjy : public TwoBodyMOIntsTransform {

  // Initialize the MO integrals accumulator
  void init_acc();
  // Compute required dynamic memory for a given batch size
  distsize_t compute_transform_dynamic_memory_(int ni) const;

  Ref<DistArray4> partially_tformed_ints_;  // partially transformed integrals

public:

  TwoBodyMOIntsTransform_ixjy(StateIn&);
  TwoBodyMOIntsTransform_ixjy(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                              const Ref<TwoBodyIntDescr>& tbintdescr,
                              const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                              const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4);
  ~TwoBodyMOIntsTransform_ixjy();

  void save_data_state(StateOut&);

  /// Implementation of TwoBodyMOIntsTransform::type()
  std::string type() const { return "ixjy"; }

  /** Returns the number of bytes allocated for each ij-block of integrals of one type
    in MemoryGrp */
  size_t memgrp_blksize() const;

  /// Computes transformed integrals
  void compute();
  /// Check symmetry of transformed integrals
  void check_int_symm(double threshold = TwoBodyMOIntsTransform::zero_integral);

  void partially_transformed_ints(const Ref<DistArray4>&);

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


