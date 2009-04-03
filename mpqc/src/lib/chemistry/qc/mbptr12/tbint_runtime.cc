//
// tbint_runtime.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <cstdlib>
#include <sstream>
#include <cassert>
#include <chemistry/qc/mbptr12/tbint_runtime.h>
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

ParsedTwoBodyFourCenterIntKey::ParsedTwoBodyFourCenterIntKey(const std::string& key) :
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
  ExEnv::out0() << indent << "ParsedTwoBodyFourCenterIntKey::ParsedTwoBodyFourCenterIntKey():" << std::endl << incindent;
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
ParsedTwoBodyFourCenterIntKey::key(const std::string& bra1,
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
ParsedTwoBodyFourCenterIntKey::key(const std::string& bra1,
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
ParsedTwoBodyFourCenterIntKey::key(const Ref<TwoBodyIntDescr>& descr)
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
ParsedTwoBodyFourCenterIntKey::create_descr(const std::string& oper_key,
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
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyFourCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12(integral,params_cast);
  }
  if (oper_key == std::string("G12'")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyFourCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12NC(integral,params_cast);
  }
  if (oper_key == std::string("G12DKH")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyFourCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrG12DKH(integral,params_cast);
  }
  if (oper_key == std::string("GenG12")) {
    Ref<IntParamsGenG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyFourCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyIntDescrGenG12(integral,params_cast);
  }
  throw ProgrammingError("ParsedTwoBodyFourCenterIntKey::create_descr() -- unknown oper",
                         __FILE__,__LINE__);
}

////

ParsedTwoBodyThreeCenterIntKey::ParsedTwoBodyThreeCenterIntKey(const std::string& key) :
  pkey_(key)
{
}

std::string
ParsedTwoBodyThreeCenterIntKey::key(const std::string& bra1,
                       const std::string& bra2,
                       const std::string& ket1,
                       const std::string& oper,
                       const std::string& params) {
  return ParsedTwoBodyFourCenterIntKey::key(bra1,bra2,ket1,"",oper,params,"");
}
std::string
ParsedTwoBodyThreeCenterIntKey::key(const std::string& bra1,
                                    const std::string& bra2,
                                    const std::string& ket1,
                                    const std::string& descr) {
  return ParsedTwoBodyFourCenterIntKey::key(bra1,bra2,ket1,"",descr,"");
}
std::string
ParsedTwoBodyThreeCenterIntKey::key(const Ref<TwoBodyThreeCenterIntDescr>& descr) {
    Ref<TwoBodyThreeCenterIntDescrERI> eridescr; eridescr << descr;
    Ref<TwoBodyThreeCenterIntDescrR12> r12descr; r12descr << descr;
    Ref<TwoBodyThreeCenterIntDescrG12> g12descr; g12descr << descr;
    Ref<TwoBodyThreeCenterIntDescrG12NC> g12ncdescr; g12ncdescr << descr;
    Ref<TwoBodyThreeCenterIntDescrG12DKH> g12dkhdescr; g12dkhdescr << descr;

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

Ref<TwoBodyThreeCenterIntDescr>
ParsedTwoBodyThreeCenterIntKey::create_descr(const std::string& oper_key,
                                             const Ref<IntParams>& p,
                                             const Ref<Integral>& integral) {
  if (oper_key == std::string("ERI")) {
    return new TwoBodyThreeCenterIntDescrERI(integral);
  }
  if (oper_key == std::string("R12")) {
    return new TwoBodyThreeCenterIntDescrR12(integral);
  }
  if (oper_key == std::string("G12")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyThreeCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyThreeCenterIntDescrG12(integral,params_cast);
  }
  if (oper_key == std::string("G12'")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyThreeCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyThreeCenterIntDescrG12NC(integral,params_cast);
  }
  if (oper_key == std::string("G12DKH")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyThreeCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyThreeCenterIntDescrG12DKH(integral,params_cast);
  }
  if (oper_key == std::string("GenG12")) {
    Ref<IntParamsGenG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyThreeCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyThreeCenterIntDescrGenG12(integral,params_cast);
  }
  throw ProgrammingError("ParsedTwoBodyThreeCenterIntKey::create_descr() -- unknown oper",
                         __FILE__,__LINE__);
}

////

ParsedTwoBodyTwoCenterIntKey::ParsedTwoBodyTwoCenterIntKey(const std::string& key) :
  pkey_(key)
{
}

std::string
ParsedTwoBodyTwoCenterIntKey::key(const std::string& bra1,
                       const std::string& bra2,
                       const std::string& oper,
                       const std::string& params) {
  return ParsedTwoBodyFourCenterIntKey::key(bra1,bra2,"","",oper,params,"");
}
std::string
ParsedTwoBodyTwoCenterIntKey::key(const std::string& bra1,
                                    const std::string& bra2,
                                    const std::string& descr) {
  return ParsedTwoBodyFourCenterIntKey::key(bra1,bra2,"","",descr,"");
}
std::string
ParsedTwoBodyTwoCenterIntKey::key(const Ref<TwoBodyTwoCenterIntDescr>& descr) {
    Ref<TwoBodyTwoCenterIntDescrERI> eridescr; eridescr << descr;
    Ref<TwoBodyTwoCenterIntDescrR12> r12descr; r12descr << descr;
    Ref<TwoBodyTwoCenterIntDescrG12> g12descr; g12descr << descr;
    Ref<TwoBodyTwoCenterIntDescrG12NC> g12ncdescr; g12ncdescr << descr;
    Ref<TwoBodyTwoCenterIntDescrG12DKH> g12dkhdescr; g12dkhdescr << descr;

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

Ref<TwoBodyTwoCenterIntDescr>
ParsedTwoBodyTwoCenterIntKey::create_descr(const std::string& oper_key,
                                             const Ref<IntParams>& p,
                                             const Ref<Integral>& integral) {
  if (oper_key == std::string("ERI")) {
    return new TwoBodyTwoCenterIntDescrERI(integral);
  }
  if (oper_key == std::string("R12")) {
    return new TwoBodyTwoCenterIntDescrR12(integral);
  }
  if (oper_key == std::string("G12")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyTwoCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyTwoCenterIntDescrG12(integral,params_cast);
  }
  if (oper_key == std::string("G12'")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyTwoCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyTwoCenterIntDescrG12NC(integral,params_cast);
  }
  if (oper_key == std::string("G12DKH")) {
    Ref<IntParamsG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyTwoCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyTwoCenterIntDescrG12DKH(integral,params_cast);
  }
  if (oper_key == std::string("GenG12")) {
    Ref<IntParamsGenG12> params_cast; params_cast << p;
    if (params_cast.null()) throw ProgrammingError("ParsedTwoBodyTwoCenterIntKey::create_descr() -- mismatch between oper and param",__FILE__,__LINE__);
    return new TwoBodyTwoCenterIntDescrGenG12(integral,params_cast);
  }
  throw ProgrammingError("ParsedTwoBodyTwoCenterIntKey::create_descr() -- unknown oper",
                         __FILE__,__LINE__);
}

/////////////////////////////////////////////////////////////////////////////

Ref<ParamsRegistry> ParamsRegistry::instance_ = new ParamsRegistry;

const Ref<ParamsRegistry>&
ParamsRegistry::instance() { return instance_; }

ParamsRegistry::ParamsRegistry() : params_(RegistryType::instance()){
  // associate empty params key with IntParamsVoid
  Ref<IntParams> voidparams = new IntParamsVoid;
  const std::string key("");
  params_->add(key,voidparams);
}

std::string
ParamsRegistry::key(const Ref<IntParams>& params) const
{
  if (params_->value_exists(params))
    return params_->key(params);
  else  // if not found, register
    return this->add(params);
}

Ref<IntParams>
ParamsRegistry::value(const std::string& key) const
{
  if (params_->key_exists(key))
    return params_->value(key);
  else
    return 0;
}

std::string
ParamsRegistry::add(const Ref<IntParams>& params) const
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
ParamsRegistry::add(const std::string& key, const Ref<IntParams>& params) const
{
  assert(! params_->key_exists(key));
  assert(! params_->value_exists(params));
  params_->add(key,params);
}

/////////////////////////////////////////////////////////////////////////////

TwoBodyIntLayout::TwoBodyIntLayout(const std::string& str)
{
  if (str == std::string("(b1 b2|k1 k2)")) {
    type_ = _b1b2_k1k2;
  }
  else if (str == std::string("(b1 k1|b2 k2)")) {
    type_ = _b1k1_b2k2;
  }
  else
    throw ProgrammingError("TwoBodyIntLayout::Layout() -- unknown initializer string",__FILE__,__LINE__);
}

TwoBodyIntLayout::TwoBodyIntLayout(const TwoBodyIntLayout& other) :
  type_(other.type_)
{
}

TwoBodyIntLayout::operator std::string() {
  switch (type_) {
    case _b1b2_k1k2:
      return std::string("(b1 b2|k1 k2)");
    case _b1k1_b2k2:
      return std::string("(b1 k1|b2 k2)");
  }
}

TwoBodyIntLayout&
TwoBodyIntLayout::operator=(const TwoBodyIntLayout& other)
{
  type_ = other.type_;
  return *this;
}

bool
TwoBodyIntLayout::operator==(const TwoBodyIntLayout& other) const
{
  return type_ == other.type_;
}

TwoBodyIntLayout TwoBodyIntLayout::b1b2_k1k2(std::string("(b1 b2|k1 k2)"));
TwoBodyIntLayout TwoBodyIntLayout::b1k1_b2k2(std::string("(b1 k1|b2 k2)"));

/////////////////////////////////////////////////////////////////////////////

namespace sc{   namespace detail {

  /// Convert 2 spaces to SpinCase2
  SpinCase2
  spincase2(const Ref<OrbitalSpace>& space1,
            const Ref<OrbitalSpace>& space2)
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

namespace sc {

template <>
const TwoBodyMOIntsRuntime<4>::TwoBodyIntEvalRef&
TwoBodyMOIntsRuntime<4>::create_eval(const std::string& key)
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
  Ref<OrbitalSpaceRegistry> idxreg = OrbitalSpaceRegistry::instance();
  Ref<OrbitalSpace> bra1 = idxreg->value(bra1_str);
  Ref<OrbitalSpace> bra2 = idxreg->value(bra2_str);
  Ref<OrbitalSpace> ket1 = idxreg->value(ket1_str);
  Ref<OrbitalSpace> ket2 = idxreg->value(ket2_str);
  Ref<TwoBodyIntDescr> descr = create_descr(oper_str, params_str);
  factory()->set_spaces(bra1,ket1,bra2,ket2);    // factory assumes chemists' convention

#define USE_DENSITY_FITTING 1
#if USE_DENSITY_FITTING
  factory()->df_info(this->params());
#endif

  // create the transform
  Ref<TwoBodyIntEval> tform;
  const TwoBodyIntLayout layout(pkey.layout());
  if(layout == TwoBodyIntLayout::b1b2_k1k2) {
    Ref<AOSpaceRegistry> aoidxreg = AOSpaceRegistry::instance();

    // is this a partial transform?
    if (aoidxreg->value_exists(ket1) &&
        aoidxreg->value_exists(ket2))
      tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_iRjS,key,descr);
    else { // if not, look for a partial transform

      Ref<OrbitalSpace> aoket1 = aoidxreg->value(ket1->basis());
      Ref<OrbitalSpace> aoket2 = aoidxreg->value(ket2->basis());
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

      if (evals_->key_exists(half_tform_key)) { // partially tformed integrals exist, use them
        tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ixjy,key,descr);
        tform->partially_transformed_ints( evals_->value(half_tform_key)->ints_acc() );
      }
#if ALWAYS_USE_PARTIAL_TRANSFORMS
      else { // create partial transform and use it
        const Ref<TwoBodyIntEval>& half_tform = create_eval(half_tform_key);
        half_tform->compute();
        return create_eval(key);
      }
#else
      else { // decide the best algorithm
#if !ALWAYS_USE_IXJY
        if (ket1->rank() <= ket1->basis()->nbasis()) {
          tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ikjy,key,descr);
        }
        else
#endif
        {
          tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ixjy,key,descr);
        }
      }
#endif
    }
  }
  else if (layout == TwoBodyIntLayout::b1k1_b2k2)
    tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ijxy,key,descr);

  // add to the map
  evals_->add(key,tform);

  return evals_->value(key);
}

template <>
const TwoBodyMOIntsRuntime<3>::TwoBodyIntEvalRef&
TwoBodyMOIntsRuntime<3>::create_eval(const std::string& key)
{
  // parse the key
  ParsedTwoBodyIntKey pkey(key);
  const std::string& bra1_str = pkey.bra1();
  const std::string& bra2_str = pkey.bra2();
  const std::string& ket1_str = pkey.ket1();
  const std::string& oper_str = pkey.oper();
  const std::string& params_str = pkey.params();

  // get the spaces and construct the descriptor
  Ref<OrbitalSpaceRegistry> idxreg = OrbitalSpaceRegistry::instance();
  Ref<OrbitalSpace> bra1 = idxreg->value(bra1_str);
  Ref<OrbitalSpace> bra2 = idxreg->value(bra2_str);
  Ref<OrbitalSpace> ket1 = idxreg->value(ket1_str);
  factory()->set_spaces(bra1,ket1,bra2,0);    // factory assumes chemists' convention
  Ref<TwoBodyIntDescr> descr = create_descr(oper_str, params_str);

  // create the transform
  Ref<TwoBodyIntEval> tform =
    factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ijR,key,descr);

  // add to the map
  evals_->add(key,tform);

  return evals_->value(key);
}

template <>
const TwoBodyMOIntsRuntime<2>::TwoBodyIntEvalRef&
TwoBodyMOIntsRuntime<2>::create_eval(const std::string& key)
{
  // parse the key
  ParsedTwoBodyIntKey pkey(key);
  const std::string& bra1_str = pkey.bra1();
  const std::string& bra2_str = pkey.bra2();
  const std::string& oper_str = pkey.oper();
  const std::string& params_str = pkey.params();

  // currently can only handle Coulomb integrals
  assert(oper_str == "ERI");

  // get the spaces and construct the descriptor
  Ref<OrbitalSpaceRegistry> idxreg = OrbitalSpaceRegistry::instance();
  Ref<OrbitalSpace> bra1 = idxreg->value(bra1_str);
  Ref<OrbitalSpace> bra2 = idxreg->value(bra2_str);
  Ref<TwoBodyIntDescr> descr = create_descr(oper_str, params_str);

  // compute the matrix
  RefSCMatrix result;
  {
    Ref<GaussianBasisSet> brabas = bra1->basis();
    Ref<GaussianBasisSet> ketbas = bra2->basis();
    Ref<Integral> localints = factory()->integral()->clone();
    localints->set_basis(brabas,ketbas);

    // form 2-center Coulomb in AO basis
    RefSCMatrix ao(brabas->basisdim(), ketbas->basisdim(), brabas->matrixkit());
    ao.assign(0.0);
    Ref<SCElementOp> sc =
        new TwoBodyTwoCenterIntOp(localints->electron_repulsion2());
    ao.element_op(sc);
    sc = 0;

    result = ao;
  }

  // add to the map
  evals_->add(key,result);

  return evals_->value(key);
}

}; // end of namespace sc

/////////////////////////////////////////////////////////////////////////////

ClassDesc
TwoBodyMOIntsRuntimeUnion23::class_desc_(typeid(this_type),
                           "TwoBodyMOIntsRuntimeUnion23",
                           1,
                           "virtual public SavableState",
                           0,
                           0,
                           create<this_type> );

TwoBodyMOIntsRuntimeUnion23::TwoBodyMOIntsRuntimeUnion23(const Ref<MOIntsTransformFactory>& factory,
                                                         const Ref<TwoBodyTwoCenterMOIntsRuntime>& r2c,
                                                         const Ref<TwoBodyThreeCenterMOIntsRuntime>& r3c) :
  factory_(factory),
  runtime_2c_(r2c.null() ? Ref<TwoBodyTwoCenterMOIntsRuntime>(new TwoBodyTwoCenterMOIntsRuntime(factory_)) : r2c),
  runtime_3c_(r3c.null() ? Ref<TwoBodyThreeCenterMOIntsRuntime>(new TwoBodyThreeCenterMOIntsRuntime(factory_)) : r3c) {
}

TwoBodyMOIntsRuntimeUnion23::~TwoBodyMOIntsRuntimeUnion23() {
}

TwoBodyMOIntsRuntimeUnion23::TwoBodyMOIntsRuntimeUnion23(StateIn& si) {
  factory_ << SavableState::restore_state(si);
  runtime_2c_ << SavableState::restore_state(si);
  runtime_3c_ << SavableState::restore_state(si);
}

void
TwoBodyMOIntsRuntimeUnion23::save_data_state(StateOut& so) {
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(runtime_2c_.pointer(),so);
  SavableState::save_state(runtime_3c_.pointer(),so);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
