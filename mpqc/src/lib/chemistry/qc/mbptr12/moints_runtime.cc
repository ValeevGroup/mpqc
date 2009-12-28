//
// moints_runtime.cc
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
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/moints_runtime.h>

using namespace sc;

ClassDesc
MOIntsRuntime::class_desc_(typeid(this_type),
                           "MOIntsRuntime",
                           1,
                           "virtual public SavableState",
                           0,
                           0,
                           create<this_type> );

MOIntsRuntime::MOIntsRuntime(const Ref<MOIntsTransformFactory>& factory,
                             const Ref<GaussianBasisSet>& dfbasis) :
  factory_(factory), dfbasis_(dfbasis),
  runtime_2c_(new TwoBodyTwoCenterMOIntsRuntime(factory_)),
  runtime_3c_(new TwoBodyThreeCenterMOIntsRuntime(factory_)),
  runtime_4c_(new TwoBodyFourCenterMOIntsRuntime(factory_))
{
  if (dfbasis_.nonnull()) {
    runtime_df_ = new DensityFittingRuntime(new DensityFitting::MOIntsRuntime(factory_,runtime_2c_,runtime_3c_));
    typedef TwoBodyFourCenterMOIntsRuntime::Params DFParams;
    const DFParams* params = new DFParams(dfbasis_, dfbasis_, runtime_df_);
    runtime_4c_->params(params);
  }
}

MOIntsRuntime::~MOIntsRuntime() {
}

MOIntsRuntime::MOIntsRuntime(StateIn& si) {
  factory_ << SavableState::restore_state(si);
  dfbasis_ << SavableState::restore_state(si);
  runtime_2c_ << SavableState::restore_state(si);
  runtime_3c_ << SavableState::restore_state(si);
  runtime_4c_ << SavableState::restore_state(si);
}

void
MOIntsRuntime::save_data_state(StateOut& so) {
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(dfbasis_.pointer(),so);
  SavableState::save_state(runtime_2c_.pointer(),so);
  SavableState::save_state(runtime_3c_.pointer(),so);
  SavableState::save_state(runtime_4c_.pointer(),so);
}

void
MOIntsRuntime::obsolete() {
  runtime_2c_->obsolete();
  runtime_3c_->obsolete();
  runtime_4c_->obsolete();
  if (runtime_df_.nonnull())
    runtime_df_->obsolete();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
