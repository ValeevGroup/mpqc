//
// df_runtime.cc
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

#include <cassert>
#include <chemistry/qc/lcao/df_runtime.h>
#include <math/optimize/gaussianfit.h>
#include <math/optimize/gaussianfit.timpl.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////////

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

ParsedDensityFittingKey::ParsedDensityFittingKey(const std::string& key) :
  key_(key)
{
  std::string keycopy(key);

  // pop off the leading '('
  MPQC_ASSERT(keycopy[0] == '(');
  keycopy.erase(keycopy.begin());
  // get space1
  space1_ = pop_till_token(keycopy,' ');
  // get space2
  space2_ = pop_till_token(keycopy,'|');
  // get rid of the "DF" token
  std::string crap = pop_till_token(keycopy,'(');
  // get kernel
  std::string kernel = pop_till_token(keycopy,')');
  kernel_pkey_ = ParsedTwoBodyOperKey(kernel);
  // get rid of |
  crap = pop_till_token(keycopy,'|');
  // get fspace
  fspace_ = pop_till_token(keycopy,')');

#if 0
  ExEnv::out0() << indent << "ParsedDensityFittingKey::ParsedDensityFittingKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "space1 = " << space1_ << std::endl;
  ExEnv::out0() << indent << "space2 = " << space2_ << std::endl;
  ExEnv::out0() << indent << "fspace = " << fspace_ << std::endl;
  ExEnv::out0() << indent << "kernel = " << kernel_pkey_.key() << std::endl;
#endif
}

std::string
ParsedDensityFittingKey::key(const std::string& space1,
                             const std::string& space2,
                             const std::string& fspace,
                             const std::string& kernel)
{
  std::ostringstream oss;
  oss << "(" << space1 << " " << space2 << "|DF(" << kernel << ")|" << fspace << ")";
  return oss.str();
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingRuntime::class_desc_(typeid(this_type),
                                   "DensityFittingRuntime",
                                   1,
                                   "virtual public SavableState", 0, 0,
                                   create<this_type> );

DensityFittingRuntime::DensityFittingRuntime(const Ref<MOIntsRuntime>& r,
                                             const DensityFittingParams* dfp) :
  moints_runtime_(r),
  dfparams_(dfp),
  results_(ResultRegistry::instance())
{
}

DensityFittingRuntime::DensityFittingRuntime(StateIn& si)
{
  moints_runtime_ << SavableState::restore_state(si);
  Ref<DensityFittingParams> dfp; dfp << SavableState::restore_state(si);
  dfp->reference();
  dfparams_ = dfp.pointer();

  results_ = ResultRegistry::restore_instance(si);
}

void
DensityFittingRuntime::save_data_state(StateOut& so)
{
  SavableState::save_state(moints_runtime_.pointer(),so);
  SavableState::save_state(const_cast<DensityFittingParams*>(dfparams_),so);
  ResultRegistry::save_instance(results_,so);
}

void
DensityFittingRuntime::obsolete() {
  results_->clear();
}

bool
DensityFittingRuntime::exists(const std::string& key) const
{
  return results_->key_exists(key);
}

struct ParsedDensityFittingKeyInvolvesSpace {
    ParsedDensityFittingKeyInvolvesSpace(const std::string& skey) : space_key(skey) {}
    bool operator()(const std::pair<std::string, sc::Ref<sc::DistArray4> >& i) const {
      const ParsedDensityFittingKey pkey(i.first);
      return pkey.space1() == space_key || pkey.space2() == space_key || pkey.fspace() == space_key;
    }
    std::string space_key;
};

void
DensityFittingRuntime::remove_if(const std::string& space_key)
{
  ParsedDensityFittingKeyInvolvesSpace pred(space_key);
  results_->remove_if(pred);
}

DensityFittingRuntime::ResultRef
DensityFittingRuntime::get(const std::string& key)
{
  if (results_->key_exists(key)) {
    return results_->value(key);
  }
  else {  // if not found
    try { ParsedResultKey parsedkey(key); }
    catch (...) {
      std::ostringstream oss;
      oss << "DensityFittingRuntime::get() -- key " << key << " does not match the format";
      throw ProgrammingError(oss.str().c_str(),__FILE__,__LINE__);
    }
    // then create evaluated tform
    const ResultRef& result = create_result(key);
    return result;
  }
  MPQC_ASSERT(false); // unreachable
}

#define USE_TRANSFORMED_DF 1
#define ALWAYS_MAKE_AO2_DF 1

const DensityFittingRuntime::ResultRef&
DensityFittingRuntime::create_result(const std::string& key)
{
  // parse the key
  ParsedResultKey pkey(key);
  const std::string& space1_str = pkey.space1();
  const std::string& space2_str = pkey.space2();
  const std::string& fspace_str = pkey.fspace();
  std::string dfkernel_key = pkey.kernel();
  if (dfkernel_key.empty())
    dfkernel_key = TwoBodyOper::to_string(TwoBodyOper::eri);

  // get the spaces and construct the descriptor
  Ref<OrbitalSpaceRegistry> idxreg = this->moints_runtime()->factory()->orbital_registry();
  Ref<OrbitalSpace> space1 = idxreg->value(space1_str);
  Ref<OrbitalSpace> space2 = idxreg->value(space2_str);
  Ref<OrbitalSpace> fspace = idxreg->value(fspace_str);

  //
  // create the result assuming that fits have been already created to allow getting to the target with minimal effort
  // by minimal effort I mean: permutation or transforming one of spaces from AO basis. Thus the procedure is:
  // 1) look for (space2 space1|
  // 2) look for (space1 AO(space2)|, (space2 AO(space1)|, (AO(space1) space2|
  // 3) else construct from scratch
  //   a) make (space_i space_j| where i and j are determined such that space_i is an AO space, if possible
  //      (to make 3-center integrals easier); if space1 and space2 and both AO choose space_i such that
  //      rank(space_i) = min(rank(space1),rank(space2))
  //   b) call itself

  // 1) look for (space2 space1|
  {
    const std::string bkey = ParsedResultKey::key(space2->id(), space1->id(), fspace->id(),
                                                  dfkernel_key);
    if (this->exists(bkey)) {
      Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                          space1, space2, fspace->basis(),
                                                          this->get(bkey));
      df->compute();
      ResultRef result = df->C();
      results_->add(key,result);
      return results_->value(key);
    }
  }

  Ref<AOSpaceRegistry> aoreg = this->moints_runtime()->factory()->ao_registry();

  // 2) look for existing AO-basis fittings
#if USE_TRANSFORMED_DF
  {
    Ref<OrbitalSpace> space1_ao = aoreg->value( space1->basis() );
    Ref<OrbitalSpace> space2_ao = aoreg->value( space2->basis() );
    const bool space1_is_ao = aoreg->value_exists( space1 );
    const bool space2_is_ao = aoreg->value_exists( space2 );

    // look for (space1 AO(space2)| -> compute (space1 space2| and return
    if (!space2_is_ao) {
      const std::string bkey = ParsedResultKey::key(space1->id(), space2_ao->id(), fspace->id(),
                                                    dfkernel_key);
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new TransformedDensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                               space1, space2, fspace->basis(),
                                                               this->get(bkey));
        df->compute();
        ResultRef result = df->C();
        results_->add(key,result);
        return results_->value(key);
      }
    }

    // look for (space2 AO(space1)| -> compute (space2 space1| and call itself
    if (!space1_is_ao) {
      const std::string bkey = ParsedResultKey::key(space2->id(), space1_ao->id(), fspace->id(),
                                                    dfkernel_key);
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new TransformedDensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                               space2, space1, fspace->basis(),
                                                               this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space2->id(), space1->id(), fspace->id(),
                                                      dfkernel_key);
        df->compute();
        ResultRef result = df->C();
        results_->add(tkey,result);
        return this->create_result(key);
      }
    }

    // look for (AO(space1) space2| -> compute (space2 AO(space1)| and call itself
    if (!space1_is_ao) {
      const std::string bkey = ParsedResultKey::key(space1_ao->id(), space2->id(), fspace->id(),
                                                    dfkernel_key);
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                            space2, space1_ao, fspace->basis(),
                                                            this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space2->id(), space1_ao->id(), fspace->id(),
                                                      dfkernel_key);
        df->compute();
        ResultRef result = df->C();
        results_->add(tkey,result);
        return this->create_result(key);
      }
    }

    // look for (AO(space2) space1| -> compute (space1 AO(space2)| and call itself
    if (!space2_is_ao) {
      const std::string bkey = ParsedResultKey::key(space2_ao->id(), space1->id(), fspace->id(),
                                                    dfkernel_key);
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                            space1, space2_ao, fspace->basis(),
                                                            this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space1->id(), space2_ao->id(), fspace->id(),
                                                      dfkernel_key);
        df->compute();
        ResultRef result = df->C();
        results_->add(tkey,result);
        return this->create_result(key);
      }
    }

  }
#endif

  // 3) construct from scratch
  {
#if ALWAYS_MAKE_AO2_DF && USE_TRANSFORMED_DF
    Ref<OrbitalSpace> space1_ao = aoreg->value( space1->basis() );
    Ref<OrbitalSpace> space2_ao = aoreg->value( space2->basis() );
    const bool space1_is_ao = aoreg->value_exists( space1 );
    const bool space2_is_ao = aoreg->value_exists( space2 );
    // if both spaces are MO, compute (AO(space_i) space_j|, (where rank is increased the least) call itself
    if (!space1_is_ao && !space2_is_ao) {
      Ref<OrbitalSpace> mospace, aospace;
      if (space1->rank() * space2_ao->rank() >= space1_ao->rank() * space2->rank()) {
        aospace = space1_ao;
        mospace = space2;
      }
      else {
        aospace = space2_ao;
        mospace = space1;
      }
      const std::string bkey = ParsedResultKey::key(aospace->id(), mospace->id(), fspace->id(),
                                                    dfkernel_key);
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                  aospace, mospace, fspace->basis());
      df->compute();
      results_->add(bkey, df->C());
      return this->create_result(key);
    }
    // if space1 is MO and space2 is AO, compute (space2 space_1|, call itself
    if (!space1_is_ao && space2_is_ao) {
      const std::string bkey = ParsedResultKey::key(space2->id(), space1_ao->id(), fspace->id(),
                                                    dfkernel_key);
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                  space2, space1_ao, fspace->basis());
      df->compute();
      results_->add(bkey, df->C());
      return this->create_result(key);
    }
#endif
    // otherwise just compute (space1 space2|
    {
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfkernel_key, dfparams_->solver(),
                                                  space1, space2, fspace->basis());
      df->compute();
      ResultRef result = df->C();
      results_->add(key,result);
      return results_->value(key);
    }
  }

  MPQC_ASSERT(false);  // unreachable
}

std::string
DensityFittingRuntime::default_dfbs_name(const std::string& obs_name, int incX, bool force_aug) {
  int X = 0;
  bool aug = force_aug;
  if (obs_name == "cc-pVDZ-F12")
    X = 3;
  else if (obs_name == "cc-pVTZ-F12")
    X = 4;
  else if (obs_name == "cc-pVQZ-F12")
    X = 5;
  else if (obs_name == "aug-cc-pVDZ")
  { X = 3; aug = true; }
  else if (obs_name == "aug-cc-pVTZ")
  { X = 4; aug = true; }
  else if (obs_name == "aug-cc-pVQZ")
  { X = 5; aug = true; }
  else if (obs_name == "aug-cc-pV5Z")
  { X = 6; aug = true; }

  std::string dfbs_name;
  switch (X) {
    case 2: dfbs_name = "cc-pVDZ-RI"; break;
    case 3: dfbs_name = "cc-pVTZ-RI"; break;
    case 4: dfbs_name = "cc-pVQZ-RI"; break;
    case 5: dfbs_name = "cc-pV5Z-RI"; break;
    case 6: dfbs_name = "cc-pV6Z-RI"; break;
  }
  if (aug && not dfbs_name.empty())
    dfbs_name = std::string("aug-") + dfbs_name;

  return dfbs_name;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingParams::class_desc_(typeid(DensityFittingParams),
                     "DensityFittingParams",
                     1,               // version
                     "virtual public SavableState", // must match parent
                     0, 0, create<DensityFittingParams>
                     );

DensityFittingParams::DensityFittingParams(const Ref<GaussianBasisSet>& basis,
                                           const std::string& kernel,
                                           const std::string& solver) :
                                           basis_(basis),
                                           kernel_(kernel)
{
  if (solver == "cholesky_inv")
    solver_ = DensityFitting::SolveMethod_InverseCholesky;
  else if (solver == "cholesky")
    solver_ = DensityFitting::SolveMethod_Cholesky;
  else if (solver == "cholesky_refine")
    solver_ = DensityFitting::SolveMethod_RefinedCholesky;
  else if (solver == "bunchkaufman_inv")
    solver_ = DensityFitting::SolveMethod_InverseBunchKaufman;
  else if (solver == "bunchkaufman")
      solver_ = DensityFitting::SolveMethod_BunchKaufman;
  else if (solver == "bunchkaufman_refine")
      solver_ = DensityFitting::SolveMethod_RefinedBunchKaufman;
  else
    throw ProgrammingError("invalid solver", __FILE__, __LINE__, class_desc());

  if (not kernel_.empty()) { // throw if not valid kernel
    ParsedTwoBodyOperKey kernel_pkey(kernel_);
    TwoBodyOperSet::type kernel_oper = TwoBodyOperSet::to_type(kernel_pkey.oper());
    Ref<IntParams> kernel_params = ParamsRegistry::instance()->value(kernel_pkey.params());
  }
}

DensityFittingParams::DensityFittingParams(StateIn& si) : SavableState(si) {
  basis_ << SavableState::restore_state(si);
  si.get(kernel_);
  int s; si.get(s); solver_ = static_cast<DensityFitting::SolveMethod>(s);
}

DensityFittingParams::~DensityFittingParams() {
}

void
DensityFittingParams::save_data_state(StateOut& so) {
  SavableState::save_state(basis_.pointer(), so);
  so.put(kernel_);
  so.put((int)solver_);
}

void
DensityFittingParams::print(std::ostream& o) const {
  o << indent << "Density-Fitting Parameters:" << std::endl;
  o << incindent;
    o << indent << "basis set:" << std::endl;
    o << incindent;
      basis_->print(o);
    o << decindent;
    o << indent << "kernel = " << kernel_ << std::endl;
    o << indent << "solver = ";
    switch(solver_) {
      case DensityFitting::SolveMethod_InverseBunchKaufman:  o << "Bunch-Kaufman (inverse)"; break;
      case DensityFitting::SolveMethod_BunchKaufman:         o << "Bunch-Kaufman"; break;
      case DensityFitting::SolveMethod_RefinedBunchKaufman:  o << "Bunch-Kaufman (refine)"; break;
      case DensityFitting::SolveMethod_InverseCholesky:      o << "Cholesky (inverse)"; break;
      case DensityFitting::SolveMethod_Cholesky:             o << "Cholesky"; break;
      case DensityFitting::SolveMethod_RefinedCholesky:      o << "Cholesky (refine)"; break;
      default: MPQC_ASSERT(false); // unreachable
    }
    o << std::endl;
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingInfo::class_desc_(typeid(this_type),
                                "DensityFittingInfo",
                                1,
                                "virtual public SavableState", 0, 0,
                                create<this_type> );

void
DensityFittingInfo::obsolete() {
  runtime()->obsolete();
  runtime()->moints_runtime()->runtime_2c()->obsolete();
  runtime()->moints_runtime()->runtime_3c()->obsolete();
  runtime()->moints_runtime()->runtime_2c_inv()->clear();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
