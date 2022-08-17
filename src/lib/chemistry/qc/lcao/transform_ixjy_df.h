//
// transform_ixjy_df.h
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

#ifndef _chemistry_qc_lcao_transformixjydf_h
#define _chemistry_qc_lcao_transformixjydf_h

#include <string>
#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/lcao/df_runtime.h>
#include <chemistry/qc/lcao/transform_tbint.h>

namespace sc {

  class DensityFittingInfo;

  /** TwoBodyMOIntsTransform_ixjy_df computes (ix|jy) integrals
      using parallel integral-direct density-fitting. */

class TwoBodyMOIntsTransform_ixjy_df : public TwoBodyMOIntsTransform {

  // DensityFitting objects are managed by DensityFittingRuntime
  Ref<DensityFittingRuntime> runtime_;
  Ref<GaussianBasisSet> dfbasis12_;
  Ref<GaussianBasisSet> dfbasis34_;

  // Initialize the MO integrals accumulator
  void init_acc();
  // Compute required dynamic memory for a given batch size
  distsize_t compute_transform_dynamic_memory_(int ni) const;

public:

  TwoBodyMOIntsTransform_ixjy_df(StateIn&);
  TwoBodyMOIntsTransform_ixjy_df(const std::string& name, const DensityFittingInfo* df_info,
                                 const Ref<TwoBodyIntDescr>& tbintdescr,
                                 const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                 const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4);
  ~TwoBodyMOIntsTransform_ixjy_df();

  void save_data_state(StateOut&);

  /// Implementation of TwoBodyMOIntsTransform::type()
  std::string type() const { return "ixjy_df"; }

  /** Returns the number of bytes allocated for each ij-block of integrals of one type
    in MemoryGrp */
  size_t memgrp_blksize() const;

  /// Computes transformed integrals
  void compute();
  /// Check symmetry of transformed integrals
  void check_int_symm(double threshold = TwoBodyMOIntsTransform::zero_integral);

  /// the runtime used to compute this object
  const Ref<DensityFittingRuntime>& runtime() const { return runtime_; }
  /// returns the basis set used to fit the product of space1 and space2
  const Ref<GaussianBasisSet>& dfbasis12() const { return dfbasis12_; }
  /// returns the basis set used to fit the product of space3 and space4
  const Ref<GaussianBasisSet>& dfbasis34() const { return dfbasis34_; }

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


