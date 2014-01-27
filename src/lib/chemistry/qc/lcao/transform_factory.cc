//
// transform_factory.cc
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

#include <stdexcept>
#include <cassert>

#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/lcao/transform_factory.h>
#include <chemistry/qc/lcao/transform_ikjy.h>
#include <chemistry/qc/lcao/transform_ijxy.h>
#include <chemistry/qc/lcao/transform_ixjy.h>
#include <chemistry/qc/lcao/transform_iRjS.h>
#include <chemistry/qc/lcao/transform_ijR.h>
#include <chemistry/qc/lcao/transform_ixjy_df.h>
#include <chemistry/qc/lcao/transform_factory.timpl.h>

// Set to 1 if want to use ixjy transforms only
#define USE_IXJY_ALWAYS 0

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  MOIntsTransformFactory
 -----------*/
static ClassDesc MOIntsTransformFactory_cd(
  typeid(MOIntsTransformFactory),"MOIntsTransformFactory",1,"virtual public SavableState",
  0, 0, create<MOIntsTransformFactory>);

MOIntsTransformFactory::MOIntsTransformFactory(const Ref<Integral>& integral) :
  integral_(integral), tbintdescr_(new DefaultTwoBodyIntDescr(integral)),
  space1_(0), space2_(0), space3_(0), space4_(0), df_info_(0)
{
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  oreg_ = OrbitalSpaceRegistry::instance();
  aoreg_ = AOSpaceRegistry::instance();

  // Default values
  debug_ = 0;
  dynamic_ = false;
  print_percent_ = 10.0;
  ints_method_ = MOIntsTransform::StoreMethod::mem_posix;
  file_prefix_ = "/tmp/moints";
  log2_precision_ = -50.0; // 2^{-50} \sim 10^{-15}
}

MOIntsTransformFactory::MOIntsTransformFactory(StateIn& si) : SavableState(si)
{
  integral_ << SavableState::restore_state(si);
  //tbintdescr_ << SavableState::restore_state(si);
  oreg_ = OrbitalSpaceRegistry::restore_instance(si);
  aoreg_ = AOSpaceRegistry::restore_instance(si);
  space1_ << SavableState::restore_state(si);
  space2_ << SavableState::restore_state(si);
  space3_ << SavableState::restore_state(si);
  space4_ << SavableState::restore_state(si);

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  si.get(debug_);
  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;
  si.get(print_percent_);
  int ints_method; si.get(ints_method);
  ints_method_ = static_cast<MOIntsTransform::StoreMethod::type>(ints_method);
  si.get(file_prefix_);
  si.get(log2_precision_);
}

MOIntsTransformFactory::~MOIntsTransformFactory()
{
}

void
MOIntsTransformFactory::save_data_state(StateOut& so)
{
  SavableState::save_state(integral_.pointer(),so);
  //SavableState::save_state(tbintdescr_.pointer(),so);
  OrbitalSpaceRegistry::save_instance(oreg_, so);
  AOSpaceRegistry::save_instance(aoreg_, so);
  SavableState::save_state(space1_.pointer(),so);
  SavableState::save_state(space2_.pointer(),so);
  SavableState::save_state(space3_.pointer(),so);
  SavableState::save_state(space4_.pointer(),so);

  so.put(debug_);
  so.put((int)dynamic_);
  so.put(print_percent_);
  so.put((int)ints_method_);
  so.put(file_prefix_);
  so.put(log2_precision_);
}

void
MOIntsTransformFactory::obsolete() {
  oreg_->clear();
  aoreg_->clear();
  space1_ = space2_ = space3_ = space4_ = 0;
}

void
MOIntsTransformFactory::set_spaces(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                   const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4)
{
  space1_ = space1;
  if (space1_.null())
    throw ProgrammingError("MOIntsTransformFactory::set_spaces() -- space1 cannot be null",__FILE__,__LINE__);
  if (space2.null())
    space2_ = space1_;
  else
    space2_ = space2;
  if (space3.null())
    space3_ = space2_;
  else
    space3_ = space3;
  if (space4.null())
    space4_ = space3_;
  else
    space4_ = space4;
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform_13(const std::string& name,
                                             const Ref<TwoBodyIntDescr>& descrarg)
{
  Ref<TwoBodyMOIntsTransform> result;
  const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);

  if (df_info_ == 0) { // no density-fitting
    if (space2_->rank() <= space2_->basis()->nbasis()) {
#if USE_IXJY_ALWAYS
      result = new TwoBodyMOIntsTransform_ixjy(name,this,descr,space1_,space2_,space3_,space4_);
#else
      result = new TwoBodyMOIntsTransform_ikjy(name,this,descr,space1_,space2_,space3_,space4_);
#endif
    }
    else {
      result = new TwoBodyMOIntsTransform_ixjy(name,this,descr,space1_,space2_,space3_,space4_);
    }
  }
  else {
    result = new TwoBodyMOIntsTransform_ixjy_df(name,df_info_,descr,space1_,space2_,space3_,space4_);
  }

  if (top_mole_)
    result->set_top_mole(top_mole_);

  return result;
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform_12(const std::string& name,
                                             const Ref<TwoBodyIntDescr>& descrarg)
{
  Ref<TwoBodyMOIntsTransform> result;
  const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);
  result = new TwoBodyMOIntsTransform_ijxy(name,this,descr,space1_,space2_,space3_,space4_);

  if (top_mole_)
    result->set_top_mole(top_mole_);

  return result;
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform(StorageType storage,
                                          const std::string& name,
                                          const Ref<TwoBodyIntDescr>& descrarg)
{
  const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);

  switch (storage) {
    case MOIntsTransform::StorageType_12:
    return twobody_transform_12(name,descr);

    case MOIntsTransform::StorageType_13:
    return twobody_transform_13(name,descr);

    default:
    throw ProgrammingError("MOIntsTransformFactory::twobody_transform() -- unknown storage type requested",__FILE__,__LINE__);
  }
}

#if 0
Ref<TwoBodyThreeCenterMOIntsTransform>
MOIntsTransformFactory::twobody_transform(const std::string& name,
                                          const Ref<TwoBodyThreeCenterIntDescr>& descr)
{
  Ref<TwoBodyThreeCenterMOIntsTransform> result = new TwoBodyThreeCenterMOIntsTransform_ijR(name,this,descr,space1_,space2_,space3_);

#if 0
  if (top_mole_)
    result->set_top_mole(top_mole_);
#endif
  return result;
}
#endif

#if 1
Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform(MOIntsTransform::TwoBodyTransformType T,
                                          const std::string& name,
                                          const Ref<TwoBodyIntDescr>& descrarg)
{
  switch (T) {
    case MOIntsTransform::TwoBodyTransformType_ixjy:
      return twobody_transform<TwoBodyMOIntsTransform_ixjy>(name,descrarg);
    case MOIntsTransform::TwoBodyTransformType_ixjy_df:
      return twobody_transform<TwoBodyMOIntsTransform_ixjy_df>(name,descrarg);
    case MOIntsTransform::TwoBodyTransformType_ikjy:
      return twobody_transform<TwoBodyMOIntsTransform_ikjy>(name,descrarg);
    case MOIntsTransform::TwoBodyTransformType_ijxy:
      return twobody_transform<TwoBodyMOIntsTransform_ijxy>(name,descrarg);
    case MOIntsTransform::TwoBodyTransformType_iRjS:
      return twobody_transform<TwoBodyMOIntsTransform_iRjS>(name,descrarg);
    case MOIntsTransform::TwoBodyTransformType_ijR:
      MPQC_ASSERT(false);
  }
  MPQC_ASSERT(false); // should be unreachable
  return Ref<TwoBodyMOIntsTransform>(); // dummy return statement to pacify picky compilers
}
#endif

Ref<TwoBodyThreeCenterMOIntsTransform>
MOIntsTransformFactory::twobody_transform(MOIntsTransform::TwoBodyTransformType T,
                                          const std::string& name,
                                          const Ref<TwoBodyThreeCenterIntDescr>& descrarg)
{
  switch (T) {
    case MOIntsTransform::TwoBodyTransformType_ijR:
      return twobody_transform<TwoBodyThreeCenterMOIntsTransform_ijR>(name,descrarg);
    default:
      ;
  }
  MPQC_ASSERT(false);  // should be unreachable
  return Ref<TwoBodyThreeCenterMOIntsTransform>(); // dummy return statement to pacify picky compilers
}

void
MOIntsTransformFactory::set_debug(int d) {
  debug_ = d;
  TwoBodyMOIntsTransform::set_debug(debug_);
  TwoBodyThreeCenterMOIntsTransform::set_debug(debug_);
}

void
MOIntsTransformFactory::set_print_percent(double pp) {
  print_percent_ = pp;
  TwoBodyMOIntsTransform::set_print_percent(print_percent_);
  TwoBodyThreeCenterMOIntsTransform::set_print_percent(print_percent_);
}

/////////////////////////////////////////////////////////////////////////////

CreateTransformHints::CreateTransformHints() :
  data_persistent_(false)
{
}

CreateTransformHints::CreateTransformHints(const CreateTransformHints& other) :
  data_persistent_(other.data_persistent_)
{
}

CreateTransformHints&
CreateTransformHints::operator=(const CreateTransformHints& other)
{
  data_persistent_ = other.data_persistent_;

  return *this;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
