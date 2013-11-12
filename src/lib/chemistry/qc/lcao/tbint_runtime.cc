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

#include <cstdlib>
#include <sstream>
#include <cassert>
#include <chemistry/qc/lcao/tbint_runtime.h>
#include <util/misc/registry.h>
#include <util/misc/registry.timpl.h>
#include <chemistry/qc/lcao/transform_ijR.h>
#include <chemistry/qc/lcao/df_runtime.h>

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

////

ParsedTwoBodyOperKey::ParsedTwoBodyOperKey() {}

ParsedTwoBodyOperKey::ParsedTwoBodyOperKey(const std::string& key) :
  key_(key)
{
  typedef std::string::const_iterator citer;
  std::string keycopy(key);

  // parse oper_plus_params
  // find oper first
  oper_ = pop_till_token(keycopy,'[');
  if (not keycopy.empty()) { // found params
    params_ = '[' + pop_till_token(keycopy,']') + ']';
  }

#if 0
  ExEnv::out0() << indent << "ParsedTwoBodyOperKey::ParsedTwoBodyOperKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "oper = " << oper_ << std::endl;
  ExEnv::out0() << indent << "params = " << params_ << std::endl << decindent;
#endif
}

std::string
ParsedTwoBodyOperKey::key(const std::string& oper,
                          const std::string& params)
{
  const std::string descr_key(oper + params);
  return descr_key;
}

////

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
  oper_pkey_ = ParsedTwoBodyOperKey(oper_plus_params);
  // get ket1
  ket1_ = pop_till_token(keycopy,' ');
  // get ket2
  ket2_ = pop_till_token(keycopy,'>');
  // layout is what's left
  layout_ = keycopy;


#if 0
  ExEnv::out0() << indent << "ParsedTwoBodyFourCenterIntKey::ParsedTwoBodyFourCenterIntKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "bra1 = " << bra1_ << std::endl;
  ExEnv::out0() << indent << "bra2 = " << bra2_ << std::endl;
  ExEnv::out0() << indent << "ket1 = " << ket1_ << std::endl;
  ExEnv::out0() << indent << "ket2 = " << ket2_ << std::endl;
  ExEnv::out0() << indent << "oper = " << oper() << std::endl;
  ExEnv::out0() << indent << "params = " << params() << std::endl;
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
  return key(bra1,bra2,ket1,ket2,ParsedTwoBodyOperKey::key(oper,params),layout);
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

/////////////////////////////////////////////////////////////////////////////

namespace sc {

template <int NumCenters> struct ParsedTwoBodyNCenterIntKeyInvolvesSpace {
    ParsedTwoBodyNCenterIntKeyInvolvesSpace(const std::string& skey) : space_key(skey) {}
    bool operator()(const std::pair<std::string, typename detail::TwoBodyIntEval<NumCenters>::refvalue >& i) const;
    std::string space_key;
};

template <> bool
ParsedTwoBodyNCenterIntKeyInvolvesSpace<4>::operator()(const std::pair<std::string, detail::TwoBodyIntEval<4>::refvalue>& i) const
{
  const ParsedTwoBodyFourCenterIntKey pkey(i.first);
  return pkey.bra1() == space_key || pkey.bra2() == space_key || pkey.ket1() == space_key || pkey.ket2() == space_key;
}

template <> bool
ParsedTwoBodyNCenterIntKeyInvolvesSpace<3>::operator()(const std::pair<std::string, detail::TwoBodyIntEval<3>::refvalue>& i) const
{
  const ParsedTwoBodyThreeCenterIntKey pkey(i.first);
  return pkey.bra1() == space_key || pkey.bra2() == space_key || pkey.ket1() == space_key;
}

template <> bool
ParsedTwoBodyNCenterIntKeyInvolvesSpace<2>::operator()(const std::pair<std::string, detail::TwoBodyIntEval<2>::refvalue>& i) const
{
  const ParsedTwoBodyTwoCenterIntKey pkey(i.first);
  return pkey.bra1() == space_key || pkey.bra2() == space_key;
}

template <>
void
TwoBodyMOIntsRuntime<4>::remove_if(const std::string& space_key)
{
  ParsedTwoBodyNCenterIntKeyInvolvesSpace<4> pred(space_key);
  evals_->remove_if(pred);
}

template <>
void
TwoBodyMOIntsRuntime<3>::remove_if(const std::string& space_key)
{
  ParsedTwoBodyNCenterIntKeyInvolvesSpace<3> pred(space_key);
  evals_->remove_if(pred);
}

template <>
void
TwoBodyMOIntsRuntime<2>::remove_if(const std::string& space_key)
{
  ParsedTwoBodyNCenterIntKeyInvolvesSpace<2> pred(space_key);
  evals_->remove_if(pred);
}

};

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

void
ParamsRegistry::clear() {
  params_->clear();
}

bool
ParamsRegistry::key_exists(const std::string& key) const {
  return params_->key_exists(key);
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
  if (str == std::string("")) {
    type_ = _b1b2_k1k2;
  }
  else if (str == std::string("(11|22)")) {
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
      return std::string("");
    case _b1k1_b2k2:
      return std::string("(11|22)");
    default:
      assert(false);
      return std::string();
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

TwoBodyIntLayout TwoBodyIntLayout::b1b2_k1k2(std::string(""));
TwoBodyIntLayout TwoBodyIntLayout::b1k1_b2k2(std::string("(11|22)"));

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
      default: break;
    }
    assert(false);
    return std::string(); // dummy return statement to pacify picky compilers
  }

}}

/////////////////////////////////////////////////////////////////////////////

namespace sc {

template <>
TwoBodyMOIntsRuntime<4>::TwoBodyMOIntsRuntime(StateIn& si)
{
  factory_ << SavableState::restore_state(si);
  evals_ << EvalRegistry::restore_instance(si);
  params_ = require_dynamic_cast<Params*>(SavableState::restore_state(si),"restored pointer of unexpected type");
}

template <>
TwoBodyMOIntsRuntime<3>::TwoBodyMOIntsRuntime(StateIn& si)
{
  factory_ << SavableState::restore_state(si);
  evals_ << EvalRegistry::restore_instance(si);
  params_ = require_dynamic_cast<Params*>(SavableState::restore_state(si),"restored pointer of unexpected type");
}

template <>
TwoBodyMOIntsRuntime<2>::TwoBodyMOIntsRuntime(StateIn& si)
{
  factory_ << SavableState::restore_state(si);
  evals_ << EvalRegistry::restore_instance(si);
  params_ = require_dynamic_cast<Params*>(SavableState::restore_state(si),"restored pointer of unexpected type");
}

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
  Ref<OrbitalSpaceRegistry> idxreg = this->factory()->orbital_registry();
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
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

// if not using density-fitting, use partial transforms
#define ALWAYS_USE_PARTIAL_TRANSFORMS 1
#define ALWAYS_USE_IXJY 1

      if (evals_->key_exists(half_tform_key)) { // partially tformed integrals exist, use them
        tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ixjy,key,descr);
        tform->partially_transformed_ints( evals_->value(half_tform_key)->ints_distarray4() );
      }
#if ALWAYS_USE_PARTIAL_TRANSFORMS
      else if (factory()->df_info() == 0) { // if not doing density-fitting create partial transform and use it
        const Ref<TwoBodyIntEval>& half_tform = create_eval(half_tform_key);
        half_tform->compute();
        return create_eval(key);
      }
#endif
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
    }
  }
  else if (layout == TwoBodyIntLayout::b1k1_b2k2)
    tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ijxy,key,descr);

  // add to the map
  evals_->add(key,tform);

  return evals_->value(key);
}

#define USE_IQR_TFORM 1
#define ALWAYS_CREATE_IQR_TFORM 1

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
  Ref<OrbitalSpaceRegistry> idxreg = this->factory()->orbital_registry();
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
  Ref<OrbitalSpace> bra1 = idxreg->value(bra1_str);
  Ref<OrbitalSpace> bra2 = idxreg->value(bra2_str);
  Ref<OrbitalSpace> ket1 = idxreg->value(ket1_str);
  Ref<TwoBodyThreeCenterIntDescr> descr = create_descr(oper_str, params_str);

  // create the transform
  Ref<TwoBodyIntEval> result;
#if USE_IQR_TFORM
  const bool ket1_is_ao = aoidxreg->value_exists(ket1);
  Ref<OrbitalSpace> ket1_ao = aoidxreg->value(ket1->basis());

  const std::string key_ao = ParsedTwoBodyIntKey::key(bra1->id(), bra2->id(),
                                                      ket1_ao->id(), oper_str, params_str);
  if (!ket1_is_ao && this->exists(key_ao)) { // partially transformed iqR transform exists (or this is a partially transformed tform)
    Ref<TwoBodyThreeCenterMOIntsTransform_ijR> iqR_tform;
    iqR_tform << this->get(key_ao);
    factory()->set_spaces(bra1,ket1,bra2,0);    // factory assumes chemists' convention
    result = new TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR(key,iqR_tform,ket1);
    // add to the map
    evals_->add(key,result);
  } else
#endif
  {
#if USE_IQR_TFORM && ALWAYS_CREATE_IQR_TFORM
    if (!ket1_is_ao) {
      // create partial transform, then try again
      Ref<OrbitalSpace> ket1_ao = aoidxreg->value(ket1->basis());
      const std::string key_ao = ParsedTwoBodyIntKey::key(bra1->id(), bra2->id(),
                                                          ket1_ao->id(), oper_str, params_str);
      factory()->set_spaces(bra1,ket1_ao,bra2,0);    // factory assumes chemists' convention
      Ref<TwoBodyIntEval> iqR_tform = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ijR,key_ao,descr);
      evals_->add(key_ao,iqR_tform);
      return this->create_eval(key);
    } else
#endif
    {
      factory()->set_spaces(bra1,ket1,bra2,0);    // factory assumes chemists' convention
      result = factory()->twobody_transform(MOIntsTransform::TwoBodyTransformType_ijR,key,descr);
      // add to the map
      evals_->add(key,result);
    }
  }

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

  Ref<OrbitalSpaceRegistry> idxreg = this->factory()->orbital_registry();
  Ref<OrbitalSpace> bra1 = idxreg->value(bra1_str);
  Ref<OrbitalSpace> bra2 = idxreg->value(bra2_str);

  // compute the matrix
  RefSCMatrix result;
  {
    Ref<GaussianBasisSet> brabas = bra1->basis();
    Ref<GaussianBasisSet> ketbas = bra2->basis();
    Ref<Integral> localints = factory()->integral()->clone();
    localints->set_basis(brabas,ketbas);
    Ref<TwoBodyIntDescr> descr = ParsedTwoBodyOperKey::create_descr<2>(oper_str,
                                                                       ParamsRegistry::instance()->value(params_str),
                                                                       localints);

    // form 2-center kernel in AO basis
    RefSCMatrix ao(brabas->basisdim(), ketbas->basisdim(), brabas->matrixkit());
    ao.assign(0.0);
    Ref<SCElementOp> sc;
    TwoBodyOperSet::type operset = TwoBodyOperSet::to_type(oper_str);
    sc = new TwoBodyTwoCenterIntOp(descr->inteval(), TwoBodyOperSetDescr::instance(operset)->opertype(0));
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
  runtime_3c_(r3c.null() ? Ref<TwoBodyThreeCenterMOIntsRuntime>(new TwoBodyThreeCenterMOIntsRuntime(factory_)) : r3c),
  runtime_2c_inv_(KernelInverseRegistry::instance())
{
}

TwoBodyMOIntsRuntimeUnion23::~TwoBodyMOIntsRuntimeUnion23() {
}

TwoBodyMOIntsRuntimeUnion23::TwoBodyMOIntsRuntimeUnion23(StateIn& si) {
  factory_ << SavableState::restore_state(si);
  runtime_2c_ << SavableState::restore_state(si);
  runtime_3c_ << SavableState::restore_state(si);
  runtime_2c_inv_ = KernelInverseRegistry::restore_instance(si);
}

void
TwoBodyMOIntsRuntimeUnion23::save_data_state(StateOut& so) {
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(runtime_2c_.pointer(),so);
  SavableState::save_state(runtime_3c_.pointer(),so);
  KernelInverseRegistry::save_instance(runtime_2c_inv_, so);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
