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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/transform_ikjy.h>
#include <chemistry/qc/mbptr12/transform_ijxy.h>
#include <chemistry/qc/mbptr12/transform_ixjy.h>
#include <chemistry/qc/mbptr12/transform_iRjS.h>
#include <chemistry/qc/mbptr12/transform_factory.timpl.h>

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

MOIntsTransformFactory::MOIntsTransformFactory(const Ref<Integral>& integral,
                                               const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                               const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4) :
  integral_(integral), tbintdescr_(new DefaultTwoBodyIntDescr(integral)),
  space1_(space1), space2_(space2), space3_(space3), space4_(space4)
{
  if (space2.null())
    space2_ = space1_;
  if (space3.null())
    space3_ = space2_;
  if (space4.null())
    space4_ = space3_;

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  // Default values
  memory_ = DEFAULT_SC_MEMORY;
  debug_ = 0;
  dynamic_ = false;
  print_percent_ = 10.0;
  ints_method_ = StoreMethod::mem_posix;
  file_prefix_ = "/tmp/moints";

  // allocate all memory at once
  mem_->set_localsize(memory_);
}

MOIntsTransformFactory::MOIntsTransformFactory(StateIn& si) : SavableState(si)
{
  integral_ << SavableState::restore_state(si);
  //tbintdescr_ << SavableState::restore_state(si);
  space1_ << SavableState::restore_state(si);
  space2_ << SavableState::restore_state(si);
  space3_ << SavableState::restore_state(si);
  space4_ << SavableState::restore_state(si);

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  double memory; si.get(memory); memory_ = (size_t) memory;
  si.get(debug_);
  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;
  si.get(print_percent_);
  int ints_method; si.get(ints_method); ints_method_ = static_cast<StoreMethod::type>(ints_method);
  si.get(file_prefix_);
}

MOIntsTransformFactory::~MOIntsTransformFactory()
{
}

void
MOIntsTransformFactory::save_data_state(StateOut& so)
{
  SavableState::save_state(integral_.pointer(),so);
  //SavableState::save_state(tbintdescr_.pointer(),so);
  SavableState::save_state(space1_.pointer(),so);
  SavableState::save_state(space2_.pointer(),so);
  SavableState::save_state(space3_.pointer(),so);
  SavableState::save_state(space4_.pointer(),so);

  so.put((double)memory_);
  so.put(debug_);
  so.put((int)dynamic_);
  so.put(print_percent_);
  so.put((int)ints_method_);
  so.put(file_prefix_);
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

  if (top_mole_.nonnull())
    result->set_top_mole(top_mole_);
  result->set_debug(debug());
  reserve_memory(result->memory());

  return result;
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform_12(const std::string& name,
                                             const Ref<TwoBodyIntDescr>& descrarg)
{
  Ref<TwoBodyMOIntsTransform> result;
  const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);
  result = new TwoBodyMOIntsTransform_ijxy(name,this,descr,space1_,space2_,space3_,space4_);

  if (top_mole_.nonnull())
    result->set_top_mole(top_mole_);
  result->set_debug(debug());
  reserve_memory(result->memory());

  return result;
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform(StorageType storage,
                                          const std::string& name,
                                          const Ref<TwoBodyIntDescr>& descrarg)
{
  const Ref<TwoBodyIntDescr> descr = (descrarg.null() ? tbintdescr() : descrarg);

  switch (storage) {
    case StorageType_12:
    return twobody_transform_12(name,descr);

    case StorageType_13:
    return twobody_transform_13(name,descr);

    default:
    throw ProgrammingError("MOIntsTransformFactory::twobody_transform() -- unknown storage type requested",__FILE__,__LINE__);
  }
}

Ref<TwoBodyMOIntsTransform>
MOIntsTransformFactory::twobody_transform(TwoBodyTransformType T,
                                          const std::string& name,
                                          const Ref<TwoBodyIntDescr>& descrarg)
{
  switch (T) {
    case TwoBodyTransformType_ixjy:
      return twobody_transform<TwoBodyMOIntsTransform_ixjy>(name,descrarg);
    case TwoBodyTransformType_ikjy:
      return twobody_transform<TwoBodyMOIntsTransform_ikjy>(name,descrarg);
    case TwoBodyTransformType_ijxy:
      return twobody_transform<TwoBodyMOIntsTransform_ijxy>(name,descrarg);
    case TwoBodyTransformType_iRjS:
      return twobody_transform<TwoBodyMOIntsTransform_iRjS>(name,descrarg);
  }
}

void
MOIntsTransformFactory::reserve_memory(size_t bytes) {
  memory_ -= bytes;
}

void
MOIntsTransformFactory::release_memory(size_t bytes) {
  memory_ += bytes;
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
