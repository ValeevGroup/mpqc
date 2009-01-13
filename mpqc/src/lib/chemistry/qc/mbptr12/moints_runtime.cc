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
#include <chemistry/qc/mbptr12/registry.h>
#include <chemistry/qc/mbptr12/registry.timpl.h>

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
  // layout is what's left
  layout_ = keycopy;

  // parse oper_plus_params
  // find oper first
  std::string oper_plus_params_copy(oper_plus_params);
  oper_ = pop_till_token(oper_plus_params_copy,'[');
  // then check if params is part of the key
  const size_t params_pos = oper_plus_params.find_first_of('[');
  if (params_pos != std::string::npos) { // found params
    params_ = oper_plus_params.substr(params_pos);
  }

#if 0
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

std::string
ParsedTwoBodyIntKey::key(const std::string& bra1,
                         const std::string& bra2,
                         const std::string& ket1,
                         const std::string& ket2,
                         const std::string& oper,
                         const std::string& params,
                         const std::string& layout)
{
  const std::string descr_key(oper + params);
  return key(bra1,bra2,ket1,ket2,descr_key,layout);
}

std::string
ParsedTwoBodyIntKey::key(const std::string& bra1,
                         const std::string& bra2,
                         const std::string& ket1,
                         const std::string& ket2,
                         const std::string& descr,
                         const std::string& layout)
{
  std::ostringstream oss;
  oss << "<" << bra1 << " " << bra2 << "|" << descr << "|" << ket1 << " " << ket2 << ">" << layout;
  return oss.str();
}


std::string
ParsedTwoBodyIntKey::key(const Ref<TwoBodyIntDescr>& descr)
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
    result = std::string("G12");
  }
  if (g12ncdescr.nonnull()) {
    result = std::string("G12'");
  }
  if (g12dkhdescr.nonnull()) {
    result = std::string("G12DKH");
  }

  return result;
}

Ref<TwoBodyIntDescr>
ParsedTwoBodyIntKey::create_descr(const std::string& oper_key,
                                  const Ref<IntParams>& p,
                                  const Ref<Integral>& integral)
{
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
  throw ProgrammingError("ParsedTwoBodyIntKey::create_descr() -- unknown oper",
                         __FILE__,__LINE__);
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc MOIntsRuntime_cd(
  typeid(MOIntsRuntime),"MOIntsRuntime",1,"virtual public SavableState",
  0, 0, create<MOIntsRuntime>);

MOIntsRuntime::MOIntsRuntime(const Ref<MOIntsTransformFactory>& f) : factory_(f),
  params_(ParamRegistry::instance()), tforms_(TformRegistry::instance())
{
  // associate empty params key with IntParamsVoid
  Ref<IntParams> voidparams = new IntParamsVoid;
  register_params(std::string(""),voidparams);
}

MOIntsRuntime::MOIntsRuntime(StateIn& si)
{
  factory_ << SavableState::restore_state(si);
  params_ = ParamRegistry::restore_instance(si);
  tforms_ = TformRegistry::restore_instance(si);
}

void
MOIntsRuntime::save_data_state(StateOut& so)
{
  SavableState::save_state(factory_.pointer(),so);
  ParamRegistry::save_instance(params_,so);
  TformRegistry::save_instance(tforms_,so);
}

bool
MOIntsRuntime::exists(const std::string& key) const
{
  return tforms_->key_exists(key);
}

Ref<MOIntsRuntime::TwoBodyIntsTransform>
MOIntsRuntime::get(const std::string& key)
{
  if (tforms_->key_exists(key)) {
    return tforms_->value(key);
  }
  else {  // if not found
    try { ParsedTwoBodyIntKey parsedkey(key); }
    catch (...) {
      std::ostringstream oss;
      oss << "MOIntsRuntime::get() -- key " << key << " does not match the format";
      throw ProgrammingError(oss.str().c_str(),__FILE__,__LINE__);
    }
    // then create tform
    const Ref<TwoBodyMOIntsTransform>& tform = create_tform(key);
    tform->compute();
    return tform;
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
  Ref<MOIndexSpace> bra1 = idxreg->value(bra1_str);
  Ref<MOIndexSpace> bra2 = idxreg->value(bra2_str);
  Ref<MOIndexSpace> ket1 = idxreg->value(ket1_str);
  Ref<MOIndexSpace> ket2 = idxreg->value(ket2_str);
  factory()->set_spaces(bra1,ket1,bra2,ket2);    // factory assumes chemists' convention
  Ref<TwoBodyIntDescr> descr = create_descr(oper_str, params_str);

  // create the transform
  Ref<TwoBodyMOIntsTransform> tform;
  const Layout layout(pkey.layout());
  if(layout == Layout_b1b2_k1k2) {
    Ref<AOIndexSpaceRegistry> aoidxreg = AOIndexSpaceRegistry::instance();
    Ref<MOIndexSpaceRegistry> idxreg = MOIndexSpaceRegistry::instance();


    // is this a partial transform?
    if (aoidxreg->value_exists(ket1) &&
        aoidxreg->value_exists(ket2))
      tform = factory()->twobody_transform(MOIntsTransformFactory::TwoBodyTransformType_iRjS,key,descr);
    else { // if not, look for a partial transform

      Ref<MOIndexSpace> aoket1 = aoidxreg->value(ket1->basis());
      Ref<MOIndexSpace> aoket2 = aoidxreg->value(ket2->basis());
      const std::string half_tform_key =
        ParsedTwoBodyIntKey::key(idxreg->key(bra1),
                                 idxreg->key(bra2),
                                 idxreg->key(aoket1),
                                 idxreg->key(aoket2),
                                 oper_str,
                                 params_str,
                                 pkey.layout());

#define ALWAYS_USE_PARTIAL_TRANSFORMS 1
#define ALWAYS_USE_IXJY 0

      if (tforms_->key_exists(half_tform_key)) { // partially tformed integrals exist, use them
        tform = factory()->twobody_transform(MOIntsTransformFactory::TwoBodyTransformType_ixjy,key,descr);
        tform->partially_transformed_ints( tforms_->value(half_tform_key)->ints_acc() );
      }
#if ALWAYS_USE_PARTIAL_TRANSFORMS
      else { // create partial transform and use it
        const Ref<TwoBodyMOIntsTransform>& half_tform = create_tform(half_tform_key);
        half_tform->compute();
        return create_tform(key);
      }
#else
      else { // decide the best algorithm
#if !ALWAYS_USE_IXJY
        if (ket1->rank() <= ket1->basis()->nbasis()) {
          tform = factory()->twobody_transform(MOIntsTransformFactory::TwoBodyTransformType_ikjy,key,descr);
        }
        else
#endif
        {
          tform = factory()->twobody_transform(MOIntsTransformFactory::TwoBodyTransformType_ixjy,key,descr);
        }
      }
#endif
    }
  }
  else if (layout == Layout_b1k1_b2k2)
    tform = factory()->twobody_transform(MOIntsTransformFactory::TwoBodyTransformType_ijxy,key,descr);

  // add to the map
  tforms_->add(key,tform);

  return tforms_->value(key);
}

std::string
MOIntsRuntime::params_key(const Ref<IntParams>& params) const
{
  if (params_->value_exists(params))
    return params_->key(params);
  else  // if not found, register
    return register_params(params);
}

Ref<IntParams>
MOIntsRuntime::params(const std::string& key) const
{
  if (params_->key_exists(key))
    return params_->value(key);
  else
    return 0;
}

std::string
MOIntsRuntime::register_params(const Ref<IntParams>& params) const
{
  std::string key;
  bool unique;
  do {
    // use a random string
    std::ostringstream oss;
#if HAVE_LRAND48
    const int random_number = ::lrand48();
#else
    const int random_number = ::rand();
#endif
    oss << "[p" << random_number << "]";
    key = oss.str();
    // make sure it's unique by searching the map for it
    unique = params_->key_exists(key) ? false : true;
    if (unique) {
      params_->add(key,params);
      return key;
    }
  } while (!unique); // if not unique -- try again
  // unreachable
  abort();
}

void
MOIntsRuntime::register_params(const std::string& key, const Ref<IntParams>& params) const
{
  assert(! params_->key_exists(key));
  assert(! params_->value_exists(params));
  params_->add(key,params);
}

std::string
MOIntsRuntime::descr_key(const Ref<TwoBodyIntDescr>& descr)
{
  std::string result = ParsedTwoBodyIntKey::key(descr);
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
  return ParsedTwoBodyIntKey::create_descr(oper_key,p,integral);
}

/////////////////////////////////////////////////////////////////////////////

MOIntsRuntime::Layout::Layout(const std::string& str)
{
  if (str == std::string("(b1 b2|k1 k2)")) {
    type_ = b1b2_k1k2;
  }
  else if (str == std::string("(b1 k1|b2 k2)")) {
    type_ = b1k1_b2k2;
  }
  else
    throw ProgrammingError("MOIntsRuntime::Layout::Layout() -- unknown initializer string",__FILE__,__LINE__);
}

MOIntsRuntime::Layout::Layout(const Layout& other) :
  type_(other.type_)
{
}

MOIntsRuntime::Layout::operator std::string() {
  switch (type_) {
    case b1b2_k1k2:
      return std::string("(b1 b2|k1 k2)");
    case b1k1_b2k2:
      return std::string("(b1 k1|b2 k2)");
  }
}

MOIntsRuntime::Layout&
MOIntsRuntime::Layout::operator=(const Layout& other)
{
  type_ = other.type_;
  return *this;
}

bool
MOIntsRuntime::Layout::operator==(const Layout& other) const
{
  return type_ == other.type_;
}

MOIntsRuntime::Layout MOIntsRuntime::Layout_b1b2_k1k2(std::string("(b1 b2|k1 k2)"));
MOIntsRuntime::Layout MOIntsRuntime::Layout_b1k1_b2k2(std::string("(b1 k1|b2 k2)"));

/////////////////////////////////////////////////////////////////////////////

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
