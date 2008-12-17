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

#include <cstdlib>
#include <sstream>
#include <cassert>
#include <chemistry/qc/mbptr12/moints_runtime.h>


using namespace sc;

namespace {

  // pop off str from beginning up to token.
  std::string
  pop_till_token(std::string& str,
                 char token) {
    const size_t next_token_pos = str.find_first_of(token);
    std::string result;
    if (next_token_pos != std::string::npos) {
      result = str.substr(0,next_token_pos);
      str.erase(0,next_token_pos+1);
    }
    else {
      result = str;
      str.clear();
    }
    return result;
  }

}

ParsedTwoBodyIntKey::ParsedTwoBodyIntKey(const std::string& key) :
  key_(key)
{
  typedef std::string::const_iterator citer;
  std::string keycopy(key);

  // pop off the leading '<'
  assert(keycopy[0] == '<');
  keycopy.erase(keycopy.begin());
  // get bra1
  bra1_ = pop_till_token(keycopy,' ');
  // get bra2
  bra2_ = pop_till_token(keycopy,'|');
  // get oper (+ params)
  const std::string oper_plus_params = pop_till_token(keycopy,'|');
  // get ket1
  ket1_ = pop_till_token(keycopy,' ');
  // get ket2
  ket2_ = pop_till_token(keycopy,'>');

  // parse oper_plus_params
  // find oper first
  std::string oper_plus_params_copy(oper_plus_params);
  oper_ = pop_till_token(oper_plus_params_copy,'[');
  // then check if params is part of the key
  const size_t params_pos = oper_plus_params.find_first_of('[');
  if (params_pos != std::string::npos) { // found params
    params_ = oper_plus_params.substr(params_pos);
  }

  // figure out the desired layout
  if (keycopy == std::string("(b1 b2|k1 k2)"))
    layout_ = MOIntsRuntime::Layout_b1b2_k1k2;
  else if (keycopy == std::string("(b1 k1|b2 k2)"))
    layout_ = MOIntsRuntime::Layout_b1k1_b2k2;
  else
    throw ProgrammingError("ParsedTwoBodyIntKey::ParsedTwoBodyIntKey() -- layout not recognized",__FILE__,__LINE__);

#if 1
  ExEnv::out0() << indent << "ParsedTwoBodyIntKey::ParsedTwoBodyIntKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "bra1 = " << bra1_ << std::endl;
  ExEnv::out0() << indent << "bra2 = " << bra2_ << std::endl;
  ExEnv::out0() << indent << "ket1 = " << ket1_ << std::endl;
  ExEnv::out0() << indent << "ket2 = " << ket2_ << std::endl;
  ExEnv::out0() << indent << "oper = " << oper_ << std::endl;
  ExEnv::out0() << indent << "params = " << params_ << std::endl;
  ExEnv::out0() << indent << "layout = " << layout_ << std::endl << decindent;
#endif
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc MOIntsRuntime_cd(
  typeid(MOIntsRuntime),"MOIntsRuntime",1,"virtual public SavableState",
  0, 0, create<MOIntsRuntime>);

MOIntsRuntime::MOIntsRuntime(const Ref<MOIntsTransformFactory>& f) : factory_(f)
{
  // initialize the random device so that I can use random() to generate keys
  ::srandomdev();

  // associate empty params key with IntParamsVoid
  Ref<IntParams> voidparams = new IntParamsVoid;
  params_map_[std::string("")] = voidparams;
}

MOIntsRuntime::MOIntsRuntime(StateIn& si)
{
  factory_ << SavableState::restore_state(si);
  // TODO implement restoring maps
  assert(false);
}

void
MOIntsRuntime::save_data_state(StateOut& so)
{
  SavableState::save_state(factory_.pointer(),so);
  // TODO implement dumping maps
  assert(false);
}

Ref<MOIntsRuntime::TwoBodyIntsAcc>
MOIntsRuntime::get(const std::string& key)
{
  TformMap::const_iterator tform_iter = tform_map_.find(key);
  if (tform_iter == tform_map_.end()) {  // if not found, create
    const Ref<TwoBodyMOIntsTransform>& tform = create_tform(key);
    tform->compute();
    return tform->ints_acc();
  }
  else {
    const Ref<TwoBodyMOIntsTransform>& tform = tform_iter->second;
    return tform->ints_acc();
  }
  assert(false); // unreachable
}


const Ref<TwoBodyMOIntsTransform>&
MOIntsRuntime::create_tform(const std::string& key)
{
  // parse the key
  ParsedTwoBodyIntKey pkey(key);
  const std::string& bra1_str = pkey.bra1();
  const std::string& bra2_str = pkey.bra2();
  const std::string& ket1_str = pkey.ket1();
  const std::string& ket2_str = pkey.ket2();
  const std::string& oper_str = pkey.oper();
  const std::string& params_str = pkey.params();

  // get the spaces and construct the descriptor
  Ref<MOIndexSpaceRegistry> idxreg = MOIndexSpaceRegistry::instance();
  Ref<MOIndexSpace> bra1 = idxreg->find(bra1_str);
  Ref<MOIndexSpace> bra2 = idxreg->find(bra2_str);
  Ref<MOIndexSpace> ket1 = idxreg->find(ket1_str);
  Ref<MOIndexSpace> ket2 = idxreg->find(ket2_str);
  factory()->set_spaces(bra1,ket1,bra2,ket2);    // factory assumes chemists' convention
  Ref<TwoBodyIntDescr> descr = create_descr(oper_str, params_str);

  // create the transform
  Ref<TwoBodyMOIntsTransform> tform;
  const MOIntsRuntime::Layout layout = pkey.layout();
  switch(layout) {
    case MOIntsRuntime::Layout_b1b2_k1k2:
      tform = factory()->twobody_transform_13(key,descr);   // factory assumes chemists' convention
      break;
    case MOIntsRuntime::Layout_b1k1_b2k2:
      tform = factory()->twobody_transform_12(key,descr);   // factory assumes chemists' convention
      break;
  }

  // add to the map
  tform_map_[key] = tform;

  // return
  return tform_map_[key];
}

std::string
MOIntsRuntime::key(const Ref<MOIndexSpace>& space1_bra,
                   const Ref<MOIndexSpace>& space2_bra,
                   const Ref<MOIndexSpace>& space1_ket,
                   const Ref<MOIndexSpace>& space2_ket,
                   const std::string& operator_key)
{
  std::ostringstream oss;
  // physicists' notation
  oss << "<" << space1_bra->id() << " " << space2_bra->id() << "|" << operator_key << "|"
      << space1_ket->id() << " " << space2_ket->id() << " >";
  // for case-insensitive file systems append spincase
  oss << "_" << detail::id(detail::spincase2(space1_bra,space2_bra));
  return oss.str();
}

std::string
MOIntsRuntime::key_mulliken(const Ref<MOIndexSpace>& space1_bra,
                            const Ref<MOIndexSpace>& space1_ket,
                            const Ref<MOIndexSpace>& space2_bra,
                            const Ref<MOIndexSpace>& space2_ket,
                            const std::string& operator_key)
{
  std::ostringstream oss;
  // physicists' notation
  oss << "<" << space1_bra->id() << " " << space2_bra->id() << "|" << operator_key << "|"
      << space1_ket->id() << " " << space2_ket->id() << ">";
  // for case-insensitive file systems append spincase
  oss << "_" << detail::id(detail::spincase2(space1_bra,space2_bra));
  return oss.str();
}

std::string
MOIntsRuntime::params_key(const Ref<IntParams>& params) const
{
  // search map for params
  typedef ParamsMap::const_iterator citer;
  for(citer v=params_map_.begin();
      v != params_map_.end();
      ++v) {
    if (v->second == params) {  // if found, return
      return v->first;
    }
  }
  // if not found, register
  return register_params(params);
}

Ref<IntParams>
MOIntsRuntime::params(const std::string& key) const
{
  // search for key
  typedef ParamsMap::const_iterator citer;
  citer v = params_map_.find(key);
  if (v == params_map_.end())
    return 0;
  else
    return v->second;
}

std::string
MOIntsRuntime::register_params(const Ref<IntParams>& params) const
{
  std::string key;
  bool unique;
  do {
    // use a random string computed using ::random()
    std::ostringstream oss;
    oss << "[p" << ::random() << "]";
    key = oss.str();
    // make sure it's unique by searching the map for it
    typedef ParamsMap::const_iterator citer;
    citer v = params_map_.find(key);
    if (v == params_map_.end()) { // fine, it's unique, do the work
      params_map_[key] = params;
      return key;
    }
    else { // not unique, try again
      unique = false;
    }
  } while (!unique);
  // unreachable
  abort();
}

void
MOIntsRuntime::register_params(const std::string& key, const Ref<IntParams>& params) const
{
  // make sure it's unique by searching the map for it
  typedef ParamsMap::const_iterator citer;
  citer v = params_map_.find(key);
  if (v == params_map_.end()) { // fine, it's unique, do the work
    params_map_[key] = params;
    return;
  }
  throw ProgrammingError("MOIntsRuntime::register_params() -- this key already exists",__FILE__,__LINE__);
}

std::string
MOIntsRuntime::descr_key(const Ref<TwoBodyIntDescr>& descr)
{
  Ref<TwoBodyIntDescrERI> eridescr; eridescr << descr;
  Ref<TwoBodyIntDescrR12> r12descr; r12descr << descr;
  Ref<TwoBodyIntDescrG12> g12descr; g12descr << descr;
  Ref<TwoBodyIntDescrG12NC> g12ncdescr; g12ncdescr << descr;
  Ref<TwoBodyIntDescrG12DKH> g12dkhdescr; g12dkhdescr << descr;

  std::string result;
  if (eridescr.nonnull()) {
    result = std::string("ERI");
  }
  if (r12descr.nonnull()) {
    result = std::string("R12");
  }
  if (g12descr.nonnull()) {
    result = std::string("G12'");
  }
  if (g12ncdescr.nonnull()) {
    result = std::string("G12'");
  }
  if (g12dkhdescr.nonnull()) {
    result = std::string("G12DKH");
  }

  Ref<IntParams> p = descr->params();
  result += params_key(p);

  return result;
}

Ref<TwoBodyIntDescr>
MOIntsRuntime::create_descr(const std::string& oper_key,
                            const std::string& params_key)
{
  Ref<IntParams> p = params(params_key);
  const Ref<Integral>& integral = factory()->integral();

  if (oper_key == std::string("ERI")) {
    return new TwoBodyIntDescrERI(integral);
  }
  if (oper_key == std::string("R12")) {
    return new TwoBodyIntDescrR12(integral);
  }
  if (oper_key == std::string("G12")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("MOIntsRuntime::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12(integral,params_cast);
  }
  if (oper_key == std::string("G12'")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("MOIntsRuntime::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12NC(integral,params_cast);
  }
  if (oper_key == std::string("G12DKH")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("MOIntsRuntime::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12DKH(integral,params_cast);
  }
  if (oper_key == std::string("GenG12")) {
    Ref<IntParamsGenG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("MOIntsRuntime::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrGenG12(integral,params_cast);
  }
  throw ProgrammingError("MOIntsRuntime::create_descr() -- unknown oper",
                         __FILE__,__LINE__);
}

namespace sc{   namespace detail {

  /// Convert 2 spaces to SpinCase2
  SpinCase2
  spincase2(const Ref<MOIndexSpace>& space1,
            const Ref<MOIndexSpace>& space2)
  {
    char id1 = space1->id()[0];
    char id2 = space2->id()[0];
    if (id1 < 'a' && id2 < 'a')
      return AlphaAlpha;
    if (id1 < 'a' && id2 >= 'a')
      return AlphaBeta;
    if (id1 >= 'a' && id2 >= 'a')
      return BetaBeta;
    throw ProgrammingError("spincase2(space1,space2) -- BetaAlpha spaces are not allowed",
                           __FILE__,__LINE__);
  }
  std::string
  id(SpinCase2 S) {
    switch(S) {
      case AlphaBeta:  return "ab";
      case AlphaAlpha:  return "aa";
      case BetaBeta:  return "bb";
    }
  }

}}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
