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

#include <chemistry/qc/lcao/moints_runtime.h>

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
                             const Ref<DensityFittingParams>& dfparams) :
  factory_(factory), dfparams_(dfparams), dfinfo_(0),
  runtime_2c_(new TwoBodyTwoCenterMOIntsRuntime(factory_)),
  runtime_3c_(new TwoBodyThreeCenterMOIntsRuntime(factory_)),
  runtime_4c_(new TwoBodyFourCenterMOIntsRuntime(factory_))
{
  if (dfparams_) {
    runtime_df_ = new DensityFittingRuntime(new DensityFitting::MOIntsRuntime(factory_,runtime_2c_,runtime_3c_),
                                            dfparams_.pointer());
    typedef TwoBodyFourCenterMOIntsRuntime::Params DFInfo;
    dfinfo_ = new DFInfo(dfparams_, runtime_df_);
    runtime_4c_->params(dfinfo_);
  }
}

MOIntsRuntime::~MOIntsRuntime() {
}

MOIntsRuntime::MOIntsRuntime(StateIn& si) {
  factory_ << SavableState::restore_state(si);
  dfparams_ << SavableState::restore_state(si);
  dfinfo_ << SavableState::restore_state(si);
  runtime_2c_ << SavableState::restore_state(si);
  runtime_3c_ << SavableState::restore_state(si);
  runtime_4c_ << SavableState::restore_state(si);
}

void
MOIntsRuntime::save_data_state(StateOut& so) {
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(dfparams_.pointer(),so);
  SavableState::save_state(dfinfo_.pointer(),so);
  SavableState::save_state(runtime_2c_.pointer(),so);
  SavableState::save_state(runtime_3c_.pointer(),so);
  SavableState::save_state(runtime_4c_.pointer(),so);
}

void
MOIntsRuntime::obsolete() {
  runtime_2c_->obsolete();
  runtime_3c_->obsolete();
  runtime_4c_->obsolete();
  if (runtime_df_)
    runtime_df_->obsolete();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
