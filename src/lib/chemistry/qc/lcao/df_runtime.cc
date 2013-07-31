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
  assert(keycopy[0] == '(');
  keycopy.erase(keycopy.begin());
  // get space1
  space1_ = pop_till_token(keycopy,' ');
  // get space2
  space2_ = pop_till_token(keycopy,'|');
  // get rid of the "DF" token
  std::string crap = pop_till_token(keycopy,'[');
  // get kernel
  kernel_ = pop_till_token(keycopy,']');
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
  ExEnv::out0() << indent << "kernel = " << kernel_ << std::endl;
#endif
}

std::string
ParsedDensityFittingKey::key(const std::string& space1,
                             const std::string& space2,
                             const std::string& fspace,
                             const std::string& kernel)
{
  std::ostringstream oss;
  oss << "(" << space1 << " " << space2 << "|DF[" << kernel << "]|" << fspace << ")";
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
  assert(false); // unreachable
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
                                                  dfparams()->kernel_key());
    if (this->exists(bkey)) {
      Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
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
                                                    dfparams()->kernel_key());
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new TransformedDensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
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
                                                    dfparams()->kernel_key());
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new TransformedDensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                               space2, space1, fspace->basis(),
                                                               this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space2->id(), space1->id(), fspace->id(),
                                                      dfparams()->kernel_key());
        df->compute();
        ResultRef result = df->C();
        results_->add(tkey,result);
        return this->create_result(key);
      }
    }

    // look for (AO(space1) space2| -> compute (space2 AO(space1)| and call itself
    if (!space1_is_ao) {
      const std::string bkey = ParsedResultKey::key(space1_ao->id(), space2->id(), fspace->id(),
                                                    dfparams()->kernel_key());
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                            space2, space1_ao, fspace->basis(),
                                                            this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space2->id(), space1_ao->id(), fspace->id(),
                                                      dfparams()->kernel_key());
        df->compute();
        ResultRef result = df->C();
        results_->add(tkey,result);
        return this->create_result(key);
      }
    }

    // look for (AO(space2) space1| -> compute (space1 AO(space2)| and call itself
    if (!space2_is_ao) {
      const std::string bkey = ParsedResultKey::key(space2_ao->id(), space1->id(), fspace->id(),
                                                    dfparams()->kernel_key());
      if (this->exists(bkey)) {
        Ref<DensityFitting> df = new PermutedDensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                            space1, space2_ao, fspace->basis(),
                                                            this->get(bkey));
        const std::string tkey = ParsedResultKey::key(space1->id(), space2_ao->id(), fspace->id(),
                                                      dfparams()->kernel_key());
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
                                                    dfparams()->kernel_key());
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                  aospace, mospace, fspace->basis());
      df->compute();
      results_->add(bkey, df->C());
      return this->create_result(key);
    }
#endif
    // if space1 is MO and space2 is AO, compute (space2 space_1|, call itself
    if (!space1_is_ao && space2_is_ao) {
      const std::string bkey = ParsedResultKey::key(space2->id(), space1_ao->id(), fspace->id(),
                                                    dfparams()->kernel_key());
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                  space2, space1_ao, fspace->basis());
      df->compute();
      results_->add(bkey, df->C());
      return this->create_result(key);
    }
    // otherwise just compute (space1 space2|
    {
      Ref<DensityFitting> df = new DensityFitting(moints_runtime_, dfparams_->kernel_key(), dfparams_->solver(),
                                                  space1, space2, fspace->basis());
      df->compute();
      ResultRef result = df->C();
      results_->add(key,result);
      return results_->value(key);
    }
  }

  assert(false);  // unreachable
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



/////////////////////////////////////////////////////////////////////////////

ClassDesc
DensityFittingParams::class_desc_(typeid(DensityFittingParams),
                     "DensityFittingParams",
                     2,               // version
                     "virtual public SavableState", // must match parent
                     0, 0, create<DensityFittingParams>
                     );

DensityFittingParams::DensityFittingParams(const Ref<GaussianBasisSet>& basis,
                                           const std::string& kernel,
                                           const std::string& solver,
                                           bool local_coulomb,
                                           bool local_exchange
                                           ) :
                                           basis_(basis),
                                           kernel_(kernel),
                                           kernel_intparams_key_("default"),
                                           local_coulomb_(local_coulomb),
                                           local_exchange_(local_exchange)
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

  if (not valid_kernel(kernel_))
    throw ProgrammingError("invalid kernel", __FILE__, __LINE__, class_desc());
}

DensityFittingParams::DensityFittingParams(StateIn& si) : SavableState(si) {
  basis_ << SavableState::restore_state(si);
  si.get(kernel_);
  si.get(kernel_intparams_key_);
  int s; si.get(s); solver_ = static_cast<DensityFitting::SolveMethod>(s);
  si.get(local_coulomb_);
  si.get(local_exchange_);
}

DensityFittingParams::~DensityFittingParams() {
}

void
DensityFittingParams::save_data_state(StateOut& so) {
  SavableState::save_state(basis_.pointer(), so);
  so.put(kernel_);
  so.put(kernel_intparams_key_);
  so.put((int)solver_);
  so.put(local_coulomb_);
  so.put(local_exchange_);
}

void
DensityFittingParams::print(std::ostream& o) const {
  o << indent << "Density-Fitting Parameters:" << std::endl;
  o << incindent;
    o << indent << "basis set:" << std::endl;
    o << incindent;
      basis_->print(o);
    o << decindent;
    if(local_coulomb_) o << "strongly local fitting active in coulomb operator" << std::endl;
    if(local_exchange_) o << "strongly local fitting active in exchange operator" << std::endl;
    o << indent << "kernel = " << kernel_ << std::endl;
    o << indent << "solver = ";
    switch(solver_) {
      case DensityFitting::SolveMethod_InverseBunchKaufman:  o << "Bunch-Kaufman (inverse)"; break;
      case DensityFitting::SolveMethod_BunchKaufman:         o << "Bunch-Kaufman"; break;
      case DensityFitting::SolveMethod_RefinedBunchKaufman:  o << "Bunch-Kaufman (refine)"; break;
      case DensityFitting::SolveMethod_InverseCholesky:      o << "Cholesky (inverse)"; break;
      case DensityFitting::SolveMethod_Cholesky:             o << "Cholesky"; break;
      case DensityFitting::SolveMethod_RefinedCholesky:      o << "Cholesky (refine)"; break;
      default: assert(false); // unreachable
    }
    o << std::endl;
  o << decindent;
}

TwoBodyOper::type
DensityFittingParams::kernel_otype() const {
  if (kernel_ == "coulomb") return TwoBodyOper::eri;
  if (kernel_ == "delta") return TwoBodyOper::delta;
  if (kernel_.find("exp") != std::string::npos) return TwoBodyOper::r12_0_g12;
  assert(false); // unreachable
}

std::string
DensityFittingParams::intparams_key() const {
  if (kernel_intparams_key_ != "default")
    return kernel_intparams_key_;

  if (kernel_ == "coulomb" || kernel_ == "delta")
    return "";
  else if (kernel_.find("exp") != std::string::npos) {
    std::string kernel_param = kernel_params(kernel_);
    std::istringstream iss(kernel_param);
    double lengthscale;
    iss >> lengthscale;
    assert(lengthscale > 0.0);
    const double gamma = 1.0/lengthscale;

    // for now, fit to 6 geminals
    const int ngtg = 6;
    typedef IntParamsG12::ContractedGeminal CorrParams;
    CorrParams params;
    using namespace sc::math;
    // use exp(-gamma*r_{12}) as the weight also
    PowerExponential1D* w = new PowerExponential1D(gamma,1,0);
    typedef GaussianFit<Slater1D,PowerExponential1D> GTGFit;
    // fit on [0,2*lengthscale]
    GTGFit gtgfit(ngtg, *w, 0.0, 2*lengthscale, 1001);
    delete w;

    // fit exp(-gamma*r_{12})
    Slater1D stg(gamma);
    typedef GTGFit::Gaussians Gaussians;
    Gaussians gtgs = gtgfit(stg);

    // feed to the constructor of CorrFactor
    typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
    typedef IntParamsG12::ContractedGeminal ContractedGeminal;
    ContractedGeminal geminal;
    typedef Gaussians::const_iterator citer;
    typedef Gaussians::iterator iter;
    for (iter g = gtgs.begin(); g != gtgs.end(); ++g) {
      geminal.push_back(*g);
    }
    Ref<IntParams> intparams = new IntParamsG12(geminal);

    std::cout << "Fit exp(-" << gamma <<"*r) to " << ngtg << " Gaussians" << std::endl;
    for(int g=0; g<ngtg; ++g) {
      std::cout << "  " << geminal[g].first << " " << geminal[g].second << std::endl;
    }

    kernel_intparams_key_ = ParamsRegistry::instance()->add(intparams);
  }
  else
    assert(false); // should be unreachable
  return kernel_intparams_key_;
}

bool
DensityFittingParams::valid_kernel(const std::string& kernel) {
  if (kernel == "coulomb")
    return true;
  if (kernel == "delta")
    return true;

  std::string::size_type s = kernel.find("exp");
  if (s != std::string::npos) {
    std::string exp_params = DensityFittingParams::kernel_params(kernel);
    if (not exp_params.empty())
      return true;
    // exponential kernel_key must have positive lengthscale
    std::istringstream iss(exp_params);
    double param;  iss >> param;
    if (param <= 0.0)
      return false;
  }

  return false;
}

std::string
DensityFittingParams::kernel_params(std::string kernel) {
  pop_till_token(kernel, '(');
  if (not kernel.empty() && kernel.find(')') != std::string::npos) {
    std::string result = pop_till_token(kernel, ')');
    if (kernel.empty())
      return result;
  }
  return "";
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
