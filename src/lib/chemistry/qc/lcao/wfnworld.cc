//
// wfnworld.cc
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

// includes go here
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/basis/petite.h>
#include <util/misc/scexception.h>
#include <util/misc/consumableresources.h>
#include <util/misc/xmlwriter.h>
#include <chemistry/qc/lcao/df_runtime.h>
#include <Eigen/Dense>
#include <math/optimize/gaussianfit.h>
#include <math/optimize/gaussianfit.timpl.h>

using namespace sc;
using boost::property_tree::ptree;

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

/*---------------
  WavefunctionWorld
 ---------------*/
static ClassDesc WavefunctionWorld_cd(
  typeid(WavefunctionWorld),"WavefunctionWorld",12,"virtual public SavableState",
  0, create<WavefunctionWorld>, create<WavefunctionWorld>);

WavefunctionWorld::WavefunctionWorld(const Ref<KeyVal>& keyval)
{
  debug_ = 0;
  print_percent_ = 10.0;

  wfn_ = dynamic_cast<Wavefunction*>(keyval->describedclassvalue("wfn").pointer());

  df_ = keyval->booleanvalue("df", KeyValValueboolean(false));

  bs_df_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("df_basis").pointer(),
      "WavefunctionWorld::WavefunctionWorld\n"
      );
  if (bs_df_.nonnull())
    df_ = true;
  if (df_ == true) {
    std::string df_kernel = keyval->stringvalue("df_kernel", KeyValValuestring("coulomb"));
    //std::string df_kernel = keyval->stringvalue("df_kernel");
    if (not df_kernel.empty()) {
      std::pair<TwoBodyOperSet::type, Ref<IntParams> > kernel = init_df_kernel(df_kernel);
      df_kernel_opertype_ = kernel.first;
      df_kernel_params_ = kernel.second;
      df_kernel_ = ParsedTwoBodyOperKey::key(TwoBodyOperSet::to_string(df_kernel_opertype_),
                                             ParamsRegistry::instance()->key(df_kernel_params_));
    }

    df_solver_ = keyval->stringvalue("df_solver", KeyValValuestring("cholesky_inv"));
    if (df_solver_ != "bunchkaufman" &&
        df_solver_ != "bunchkaufman_inv" &&
        df_solver_ != "bunchkaufman_refine" &&
        df_solver_ != "householder" &&
        df_solver_ != "householder_colpiv" &&
        df_solver_ != "householder_fullpiv" &&
        df_solver_ != "cholesky" &&
        df_solver_ != "cholesky_inv" &&
        df_solver_ != "cholesky_refine"
       )
      throw InputError("invalid value",
                       __FILE__,
                       __LINE__,
                       "df_solver",
                       df_solver_.c_str(),
                       class_desc());

    df_local_coulomb_ = keyval->booleanvalue("df_local_coulomb", KeyValValueboolean(false));
    df_local_exchange_ = keyval->booleanvalue("df_local_exchange", KeyValValueboolean(false));
    exact_diag_J_ = keyval->booleanvalue("exact_diag_J", KeyValValueboolean(false));
    exact_diag_K_ = keyval->booleanvalue("exact_diag_K", KeyValValueboolean(false));

    int nout_data = keyval->count("xml_data");
    for(int i = 0; i < nout_data; ++i){
      std::string data_name = keyval->stringvalue("xml_data", i);
      if(data_name == "df_basis"){
        out_data_.push_back(DFBasis);
      }
      else if(data_name == "df_coefficients"){
        out_data_.push_back(DFCoefficients);
      }
      else if(data_name == "df_integrals_eri"){
        out_data_.push_back(DFIntegralsERI);
      }
      else if(data_name == "exact_integrals_eri"){
        out_data_.push_back(ExactIntegralsERI);
      }
      else{
        throw InputError("invalid value",
                         __FILE__,
                         __LINE__,
                         "xml_data",
                         data_name.c_str(),
                         class_desc());

      }

    }
  }

  // Determine how to store MO integrals
  std::string ints_str = keyval->stringvalue("store_ints",KeyValValuestring("posix"));
  if (ints_str == std::string("mem")) {
    ints_method_ = StoreMethod::mem_only;
  }
  else if (ints_str == std::string("posix")) {
    ints_method_ = StoreMethod::posix;
  }
  else if (ints_str == std::string("mem-posix")) {
    ints_method_ = StoreMethod::mem_posix;
  }
  else if (ints_str == std::string("mpi")) {
#ifdef HAVE_MPIIO
    ints_method_ = StoreMethod::mpi;
#else
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- store_ints=mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else if (ints_str == std::string("mem-mpi")) {
#ifdef HAVE_MPIIO
    ints_method_ = StoreMethod::mem_mpi;
#else
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- store_ints=mem_mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else {
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- invalid value for keyword r12ints");
  }

  // Get the filename prefix to store the integrals
  const std::string default_basename_prefix = SCFormIO::fileext_to_filename(".moints");
  const std::string default_full_prefix = ConsumableResources::get_default_instance()->disk_location() +
                                          default_basename_prefix;
  ints_file_ = keyval->stringvalue("ints_file",KeyValValuestring(default_full_prefix));
  // if the last character of ints_file is '/' then append the default basename
  if (*(ints_file_.rbegin()) == '/')
    ints_file_ += default_basename_prefix;
  ExEnv::out0() << indent << "ints_file = " << ints_file_ << std::endl;

  // dynamic load balancing?
  dynamic_ = static_cast<bool>(keyval->booleanvalue("dynamic",KeyValValueboolean(false)));

  // ints precision
  // the default will indicate that ints_precision will be determined heuristically
  ints_precision_ = keyval->doublevalue("ints_precision", KeyValValuedouble(DBL_MAX));

  // world should have their own communicators, thus currently only one world will work
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  if (wfn_ != 0) {
    integral()->set_basis(wfn_->basis());
    initialize();
  }

  // allocate MemoryGrp storage if will use it for integrals
  if (ints_method_ == StoreMethod::mem_only ||
      ints_method_ == StoreMethod::mem_posix ||
      ints_method_ == StoreMethod::mem_mpi) {
    const size_t memgrp_size = 2 * ConsumableResources::get_default_instance()->memory()/3;    
    mem_->set_localsize(memgrp_size);
  }
}

WavefunctionWorld::WavefunctionWorld(StateIn& si) : SavableState(si)
{
  Ref<Wavefunction> wfn; wfn << SavableState::restore_state(si);
  wfn_ = wfn.pointer();

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  int ints_method; si.get(ints_method);
  ints_method_ = static_cast<WavefunctionWorld::StoreMethod::type>(ints_method);
  si.get(ints_file_);

  si.get(debug_);
  si.get(dynamic_);
  si.get(print_percent_);
  si.get(ints_precision_);
  si.get(df_local_coulomb_);
  si.get(exact_diag_J_);
  si.get(df_local_exchange_);
  si.get(exact_diag_K_);

  tfactory_ << SavableState::restore_state(si);
  moints_runtime_ << SavableState::restore_state(si);
  fockbuild_runtime_ << SavableState::restore_state(si);

  if (si.version(::class_desc<WavefunctionWorld>()) < 12)
    throw FeatureNotImplemented("Cannot restore WavefunctionWorld from versions prior to 12",
                                __FILE__, __LINE__);
  else
    si.get(df_);

  bs_df_ << SavableState::restore_state(si);

  // allocate MemoryGrp storage if will use it for integrals
  if (ints_method_ == StoreMethod::mem_only ||
      ints_method_ == StoreMethod::mem_posix ||
      ints_method_ == StoreMethod::mem_mpi) {
    const size_t memgrp_size = 2 * ConsumableResources::get_default_instance()->memory()/3;
    mem_->set_localsize(memgrp_size);
  }
}

WavefunctionWorld::~WavefunctionWorld()
{
    //mem_->set_localsize(0);
}

void WavefunctionWorld::save_data_state(StateOut& so)
{
  SavableState::save_state(wfn_,so);

  so.put((int)ints_method_);
  so.put(ints_file_);
  so.put(debug_);
  so.put(dynamic_);
  so.put(print_percent_);
  so.put(ints_precision_);
  so.put(df_local_coulomb_);
  so.put(exact_diag_J_);
  so.put(df_local_exchange_);
  so.put(exact_diag_K_);

  SavableState::save_state(tfactory_.pointer(),so);
  SavableState::save_state(moints_runtime_.pointer(),so);
  SavableState::save_state(fockbuild_runtime_.pointer(),so);
  so.put(df_);
  SavableState::save_state(bs_df_.pointer(),so);
}

void
WavefunctionWorld::set_wfn(Wavefunction* w) {
  if (w != wfn_) {
    wfn_ = w;
    if (this->tfactory_.null()) { // if it had not been initialized previosuly, initialize
      integral()->set_basis(wfn_->basis());
      initialize();
    }
    else { // otherwise complain
      throw ProgrammingError("should not call reset wfn of a WavefunctionWorld",
                             __FILE__, __LINE__, this->class_desc());
    }
  }
}

void
WavefunctionWorld::initialize()
{
  tfactory_ = new MOIntsTransformFactory(integral());
  tfactory_->set_dynamic(dynamic_);
  tfactory_->set_ints_method(ints_method_);
  tfactory_->set_file_prefix(ints_file_);
  double tfactory_ints_precision;
  if (ints_precision_ != DBL_MAX) { // precision provided in the constructor -- override wfn
    tfactory_ints_precision = ints_precision_;
  }
  else { // determine the precision using wfn accuracy
    // but in practice I have no idea how to do this reliably
    const double heuristic_precision = 1e-15;
    tfactory_ints_precision = heuristic_precision;
  }
  const double tfactory_ints_log2_precision = log(tfactory_ints_precision) / log(2.0);
  tfactory_->set_log2_precision( tfactory_ints_log2_precision );

  {
//    // also create AO spaces
//    Ref<OrbitalSpaceRegistry> idxreg = tfactory_->orbital_registry();
//    Ref<AOSpaceRegistry> aoidxreg = tfactory_->ao_registry();
//    Ref<Integral> localints = integral()->clone();
//    // OBS
//    Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu", "OBS(AO)", wfn()->basis(), localints);
//    idxreg->add(make_keyspace_pair(mu));
//    aoidxreg->add(mu->basis(),mu);

    // create MO integrals runtime
    Ref<DensityFittingParams> dfparams;
    if (bs_df_.nonnull()) {
      dfparams = new DensityFittingParams(
          bs_df_,
          df_kernel_,
          df_solver_,
          df_local_coulomb_,
          df_local_exchange_,
          exact_diag_J_,
          exact_diag_K_
      );
    }
    moints_runtime_ = new MOIntsRuntime(tfactory_, dfparams);
    tfactory_->df_info( const_cast<DensityFittingInfo*>(moints_runtime_->runtime_4c()->params()) );

    // to boot Fock build runtime we need densities
    // make zero densities to start with
    Ref<GaussianBasisSet> bs = wfn()->basis();
    Ref<PetiteList> plist = integral()->petite_list();
    RefSymmSCMatrix Pa = bs->so_matrixkit()->symmmatrix(plist->AO_basisdim()); Pa.assign(0.0);
    RefSymmSCMatrix Pb = Pa;
    // use Integral used by reference wfn!
    fockbuild_runtime_ = new FockBuildRuntime(this->tfactory()->orbital_registry(),
                                              this->tfactory()->ao_registry(),
                                              wfn()->basis(), Pa, Pb, integral(),
                                              wfn()->electric_field(),
                                              msg(),
                                              thr());
    fockbuild_runtime_->dfinfo( const_cast<DensityFittingInfo*>(moints_runtime_->runtime_4c()->params()) );
    const double fock_log2_precision = tfactory_ints_log2_precision;
    fockbuild_runtime_->set_log2_precision( tfactory_ints_log2_precision );

//    if (bs_df_) { // DF-BS
//      Ref<Integral> integral = moints_runtime_->factory()->integral();
//      // TODO how to generate unique labels
//      Ref<OrbitalSpace> fbs_space = new AtomicOrbitalSpace("Mu", "AO(FBS)", bs_df_, integral);
//      idxreg->add(make_keyspace_pair(fbs_space));
//      aoidxreg->add(bs_df_, fbs_space);
//    }
    initialize_ao_spaces();
  }
}

void
WavefunctionWorld::obsolete() {
  tfactory_->obsolete();
  moints_runtime_->obsolete();
  fockbuild_runtime_->obsolete();
}

void
WavefunctionWorld::initialize_ao_spaces()
{
  // reinitialize AO spaces
  Ref<OrbitalSpaceRegistry> idxreg = tfactory_->orbital_registry();
  Ref<AOSpaceRegistry> aoidxreg = tfactory_->ao_registry();
  // OBS
  if (aoidxreg->key_exists(wfn()->basis()) == false &&
      idxreg->key_exists("mu") == false) {
    Ref<Integral> localints = integral()->clone();
    Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu", "OBS(AO)", wfn()->basis(), localints);
    idxreg->add(make_keyspace_pair(mu));
    aoidxreg->add(mu->basis(),mu);
  }
  if (bs_df_) { // DF-BS
    if (aoidxreg->key_exists(bs_df_) == false &&
        idxreg->key_exists("Mu") == false) {
      Ref<Integral> localints = integral()->clone();
      // TODO how to generate unique labels
      Ref<OrbitalSpace> fbs_space = new AtomicOrbitalSpace("Mu", "AO(FBS)", bs_df_, localints);
      idxreg->add(make_keyspace_pair(fbs_space));
      aoidxreg->add(bs_df_, fbs_space);
    }
  }
}

const std::string& WavefunctionWorld::ints_file() const
{
  return ints_file_;
}

void
WavefunctionWorld::print(std::ostream& o) const {


  o << indent << "WavefunctionWorld:" << std::endl;
  o << incindent;

  if (df_) {
    Ref<DensityFittingParams> dfparams = new DensityFittingParams(bs_df_, df_kernel_, df_solver_);
    dfparams->print(o);
  }

  std::string ints_str;
  switch (ints_method_) {
  case WavefunctionWorld::StoreMethod::mem_only:
    ints_str = std::string("mem"); break;
  case WavefunctionWorld::StoreMethod::mem_posix:
    ints_str = std::string("mem-posix"); break;
  case WavefunctionWorld::StoreMethod::posix:
    ints_str = std::string("posix"); break;
#ifdef HAVE_MPIIO
  case WavefunctionWorld::StoreMethod::mem_mpi:
    ints_str = std::string("mem-mpi"); break;
  case WavefunctionWorld::StoreMethod::mpi:
    ints_str = std::string("mpi"); break;
#endif
  default:
    throw std::runtime_error("WavefunctionWorld::print -- invalid value of ints_method_");
  }
  o << indent << "How to Store Transformed Integrals: " << ints_str << std::endl;
  o << indent << "Transformed Integrals file suffix: " << ints_file_ << std::endl;
  o << decindent << std::endl;
}

ptree&
WavefunctionWorld::write_xml(
    ptree& parent,
    const XMLWriter& writer
)
{
  //----------------------------------------------------------------------------//
  // Only run this on master (should only be doing that anyway)
  ptree& child = get_my_ptree(parent);
  //----------------------------------------------------------------------------//
  for(int idata = 0; idata < out_data_.size(); ++idata){
    switch(out_data_[idata]){
    case DFBasis:
      writer.insert_child(child, bs_df_, "df_basis");
      break;
    case DFCoefficients:
    {
      /*=======================================================================================*/
      /* DFCoefficients (local and non-local)                 		                        {{{1 */ #if 1 // begin fold
      const Ref<DensityFittingInfo>& df_info = fockbuild_runtime()->dfinfo();
      std::string metric_key = df_info->params()->kernel_key();
      //----------------------------------------------------------------------------//
      bool did_local = false;
      if(df_local_coulomb_ || df_local_exchange_){
        did_local = true;
        ptree& df_tree = child.add_child("df_coefficients", ptree());
        df_tree.put("<xmlattr>.local", true);
        xml_data_local(false, df_tree, writer);
        df_tree.put("metric", metric_key);
      }
      //----------------------------------------------------------------------------//
      if(!df_local_coulomb_ || !df_local_exchange_){
        ptree& df_tree = child.add_child("df_coefficients", ptree());
        df_tree.put("<xmlattr>.local", false);
        xml_data_nonlocal(false, df_tree, writer);
        df_tree.put("metric", metric_key);
      }
      //----------------------------------------------------------------------------//
      break;
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/

    }
    case DFIntegralsERI:
    {
      /*=======================================================================================*/
      /* DFIntegralsERI (local and non-local)                  		                        {{{1 */ #if 1 // begin fold
      const Ref<DensityFittingInfo>& df_info = fockbuild_runtime()->dfinfo();
      std::string metric_key = df_info->params()->kernel_key();
      //----------------------------------------------------------------------------//
      bool did_local = false;
      if(df_local_coulomb_ || df_local_exchange_){
        did_local = true;
        ptree& df_tree = child.add_child("df_integrals_eri", ptree());
        df_tree.put("<xmlattr>.local", true);
        xml_data_local(true, df_tree, writer);
        df_tree.put("metric", metric_key);
      }
      //----------------------------------------------------------------------------//
      if(!df_local_coulomb_ || !df_local_exchange_){
        ptree& df_tree = child.add_child("df_integrals_eri", ptree());
        df_tree.put("<xmlattr>.local", false);
        xml_data_nonlocal(true, df_tree, writer);
        df_tree.put("metric", metric_key);
      }
      //----------------------------------------------------------------------------//
      break;
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
    }
    case ExactIntegralsERI:
    {
      /*=======================================================================================*/
      /* ExactIntegralsERI                                    		                        {{{1 */ #if 1 // begin fold
      ptree& eri_tree = child.add_child("eri_integrals", ptree());
      //----------------------------------------//
      const Ref<DensityFittingInfo>& df_info = fockbuild_runtime()->dfinfo();
      std::string metric_key = df_info->params()->kernel_key();
      const Ref<GaussianBasisSet> obs = wfn()->basis();
      // convenient alias for comparing with code in coulomb_df_local() and exchange_df_local() methods
      const Ref<GaussianBasisSet> dfbs = bs_df_;
      //----------------------------------------//
      const Ref<AOSpaceRegistry> ao_registry = moints_runtime()->factory()->ao_registry();
      const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
      //----------------------------------------//
      const int obsnbf = obs->nbasis();
      //----------------------------------------//
      const std::string eri_key = ParsedTwoBodyFourCenterIntKey::key(
          obs_space->id(),
          obs_space->id(),
          obs_space->id(),
          obs_space->id(),
          "ERI",
          TwoBodyIntLayout::b1k1_b2k2 // Chemists' notation
      );
      //----------------------------------------------------------------------------//
      const Ref<TwoBodyMOIntsTransform>& eri_tform = moints_runtime4()->get(eri_key);
      eri_tform->compute();
      Ref<DistArray4> g = eri_tform->ints_distarray4();  g->activate();
      unsigned int g_type_idx = TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI)->opertype(TwoBodyOper::eri);
      //----------------------------------------//
      for(int mu = 0; mu < obsnbf; ++mu){
        for(int nu = 0; nu <= mu; ++nu){
          const double* g_mu_nu_ptr = g->retrieve_pair_block(mu, nu, g_type_idx);
          Eigen::MatrixXd g_mu_nu(obsnbf, obsnbf);
          g_mu_nu = Eigen::MatrixXd::Zero(obsnbf, obsnbf);
          for(int la = 0; la < obsnbf; ++la){
            for(int si = 0; si < obsnbf; ++si){
              g_mu_nu(la, si) = g_mu_nu_ptr[la*obsnbf + si];
            }
          }
          ptree& int_parent = eri_tree.add_child("integral_block", ptree());
          int_parent.put("<xmlattr>.ao_index1", mu);
          int_parent.put("<xmlattr>.ao_index2", nu);
          writer.write_to_xml(g_mu_nu, int_parent);
          g->release_pair_block(mu, nu, g_type_idx);
        }
      }
      g->deactivate();
      break;
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
    }
    };
  }
  return child;
}

//----------------------------------------------------------------------------//
// WavefunctionWorld::write_xml() helper functions

void
WavefunctionWorld::xml_data_local(bool do_integrals, ptree& pt, const XMLWriter& writer)
{
  /*=======================================================================================*/
  /* Setup                                                		                        {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//
  const Ref<DensityFittingInfo>& df_info = fockbuild_runtime()->dfinfo();
  std::string metric_key = df_info->params()->kernel_key();
  const Ref<GaussianBasisSet> obs = wfn()->basis();
  const Ref<GaussianBasisSet> dfbs = bs_df_;
  //----------------------------------------//
  const Ref<DensityFittingRuntime> df_rtime = df_info->runtime();
  const Ref<TwoBodyThreeCenterMOIntsRuntime> int3c_rtime = df_rtime->moints_runtime()->runtime_3c();
  const Ref<TwoBodyTwoCenterMOIntsRuntime> int2c_rtime = df_rtime->moints_runtime()->runtime_2c();
  //----------------------------------------//
  const Ref<AOSpaceRegistry> ao_registry = moints_runtime()->factory()->ao_registry();
  const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
  const Ref<OrbitalSpace>& df_space = ao_registry->value(bs_df_);
  //----------------------------------------//
  const bool noncoulomb_kernel = (metric_key.find("exp") != std::string::npos);
  std::string params_key = df_info->params()->intparams_key();
  std::string operset_key = noncoulomb_kernel ? "G12'" : "ERI";
  TwoBodyOper::type metric_oper =
      noncoulomb_kernel ? TwoBodyOper::r12_0_g12 : TwoBodyOper::eri;
  unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(
      noncoulomb_kernel ? TwoBodyOperSet::G12NC : TwoBodyOperSet::ERI
  )->opertype(metric_oper);
  //----------------------------------------------------------------------------//
  std::string dfkey = ParsedDensityFittingKey::key(
      obs_space->id(),
      obs_space->id(),
      df_space->id(),
      metric_key
  );
  //----------------------------------------------------------------------------//
  const int obsnbf = obs->nbasis();
  const int dfnbf = dfbs->nbasis();
  //----------------------------------------------------------------------------//
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Transform ket if we're doing integrals               		                        {{{1 */ #if 1 // begin fold
  std::vector<std::shared_ptr<Eigen::MatrixXd> > g_tilde(obsnbf);
  if(do_integrals){
    //----------------------------------------//
    // Get the (X|Y) integrals
    const std::string coulomb2c_key = ParsedTwoBodyTwoCenterIntKey::key(
        df_space->id(),
        df_space->id(),
        "ERI", ""
    );
    RefSCMatrix coulomb_2c_ints = int2c_rtime->get(coulomb2c_key);
    //----------------------------------------//
    double* coulomb_2c_ints_ptr = allocate<double>(dfnbf*dfnbf);
    coulomb_2c_ints.convert(coulomb_2c_ints_ptr);
    Eigen::Map<Eigen::MatrixXd> X_g_Y(coulomb_2c_ints_ptr, dfnbf, dfnbf);
    //----------------------------------------//
    // Contract C_mu_nu^X with (X|Y)
    //----------------------------------------//
    for(int mu = 0; mu < obsnbf; ++mu){
      const int ishA = obs->function_to_shell(mu);
      const int atomA = obs->shell_to_center(ishA);
      const int dfnbfA = dfbs->nbasis_on_center(atomA);
      const int dfshoffA = dfbs->shell_on_center(atomA, 0);
      const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
      //----------------------------------------//
      std::shared_ptr<Eigen::MatrixXd> gtilde_mu(new Eigen::MatrixXd(obsnbf, dfnbf));
      *gtilde_mu = Eigen::MatrixXd::Zero(obsnbf, dfnbf);
      //----------------------------------------//
      // note: only loop over nu <= mu
      for(int nu = 0; nu < obsnbf; ++nu){
        const int jshB = obs->function_to_shell(nu);
        const int atomB = obs->shell_to_center(jshB);
        const int dfnbfB = dfbs->nbasis_on_center(atomB);
        const int dfshoffB = dfbs->shell_on_center(atomB, 0);
        const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
        //----------------------------------------//
        Eigen::VectorXd C(dfnbf);
        C = Eigen::VectorXd::Zero(dfnbf);
        std::shared_ptr<Eigen::VectorXd> Cpart = df_rtime->get(dfkey, mu, nu);
        C.segment(dfbfoffA, dfnbfA) += Cpart->head(dfnbfA);
        if(atomA != atomB){
          C.segment(dfbfoffB, dfnbfB) += Cpart->tail(dfnbfB);
        }
        //----------------------------------------//
        gtilde_mu->row(nu) += X_g_Y * C;
      }
      g_tilde[mu] = gtilde_mu;
    }
    deallocate(coulomb_2c_ints_ptr);
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Get the coefficients (and transform the bra if we're doing integrals)            {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//
  for(int mu = 0; mu < obsnbf; ++mu){
    const int ishA = obs->function_to_shell(mu);
    const int atomA = obs->shell_to_center(ishA);
    const int dfnbfA = dfbs->nbasis_on_center(atomA);
    const int dfshoffA = dfbs->shell_on_center(atomA, 0);
    const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
    //----------------------------------------//
    // note: only loop over nu <= mu
    for(int nu = 0; nu <= mu; ++nu){
      const int jshB = obs->function_to_shell(nu);
      const int atomB = obs->shell_to_center(jshB);
      const int dfnbfB = dfbs->nbasis_on_center(atomB);
      const int dfshoffB = dfbs->shell_on_center(atomB, 0);
      const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
      //----------------------------------------//
      Eigen::VectorXd C(dfnbf);
      C = Eigen::VectorXd::Zero(dfnbf);
      std::shared_ptr<Eigen::VectorXd> Cpart = df_rtime->get(dfkey, mu, nu);
      C.segment(dfbfoffA, dfnbfA) += Cpart->head(dfnbfA);
      if(atomA != atomB){
        C.segment(dfbfoffB, dfnbfB) += Cpart->tail(dfnbfB);
      }
      if(not do_integrals){
        ptree& coeff_parent = pt.add_child("coefficient_vector", ptree());
        coeff_parent.put("<xmlattr>.ao_index1", mu);
        coeff_parent.put("<xmlattr>.ao_index2", nu);
        writer.write_to_xml(C, coeff_parent);
      }
      else{
        Eigen::MatrixXd g_mu_nu(obsnbf, obsnbf);
        g_mu_nu = Eigen::MatrixXd::Zero(obsnbf, obsnbf);
        for(int la = 0; la < obsnbf; ++la){
          for(int si = 0; si < obsnbf; ++si){
            g_mu_nu(la, si) += g_tilde[la]->row(si) * C;
          }
        }
        //----------------------------------------//
        ptree& int_parent = pt.add_child("integral_block", ptree());
        int_parent.put("<xmlattr>.ao_index1", mu);
        int_parent.put("<xmlattr>.ao_index2", nu);
        writer.write_to_xml(g_mu_nu, int_parent);
      }
    } // end loop over nu
  } // end loop over mu
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

void
WavefunctionWorld::xml_data_nonlocal(
    bool do_integrals, ptree& pt, const XMLWriter& writer
)
{
  /*=======================================================================================*/
  /* Setup                                                		                        {{{1 */ #if 1 // begin fold
  const Ref<DensityFittingInfo>& df_info = fockbuild_runtime()->dfinfo();
  std::string metric_key = df_info->params()->kernel_key();
  const Ref<GaussianBasisSet> obs = wfn()->basis();
  // convenient alias for comparing with code in coulomb_df_local() and exchange_df_local() methods
  const Ref<GaussianBasisSet> dfbs = bs_df_;
  //----------------------------------------//
  const Ref<DensityFittingRuntime> df_rtime = df_info->runtime();
  const Ref<TwoBodyTwoCenterMOIntsRuntime> int2c_rtime = df_rtime->moints_runtime()->runtime_2c();
  //----------------------------------------//
  const Ref<AOSpaceRegistry> ao_registry = moints_runtime()->factory()->ao_registry();
  const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
  const Ref<OrbitalSpace>& df_space = ao_registry->value(bs_df_);
  //----------------------------------------//
  const bool noncoulomb_kernel = (metric_key.find("exp") != std::string::npos);
  std::string params_key = df_info->params()->intparams_key();
  std::string operset_key = noncoulomb_kernel ? "G12'" : "ERI";
  TwoBodyOper::type metric_oper =
      noncoulomb_kernel ? TwoBodyOper::r12_0_g12 : TwoBodyOper::eri;
  unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(
      noncoulomb_kernel ? TwoBodyOperSet::G12NC : TwoBodyOperSet::ERI
  )->opertype(metric_oper);
  //----------------------------------------------------------------------------//
  const int obsnbf = obs->nbasis();
  const int dfnbf = dfbs->nbasis();
  //----------------------------------------------------------------------------//
  const std::string C_key = ParsedDensityFittingKey::key(
      obs_space->id(),
      obs_space->id(),
      df_space->id(),
      metric_key
  );
  //----------------------------------------------------------------------------//
  Ref<DistArray4> C = df_rtime->get(C_key); C->activate();
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Transform ket if we're doing integrals               		                        {{{1 */ #if 1 // begin fold
  std::vector<std::shared_ptr<Eigen::MatrixXd> > g_tilde(obsnbf);
  if(do_integrals){
    //----------------------------------------//
    // Get the (X|Y) integrals
    const std::string coulomb2c_key = ParsedTwoBodyTwoCenterIntKey::key(
        df_space->id(),
        df_space->id(),
        "ERI", ""
    );
    RefSCMatrix coulomb_2c_ints = int2c_rtime->get(coulomb2c_key);
    //----------------------------------------//
    double* coulomb_2c_ints_ptr = allocate<double>(dfnbf*dfnbf);
    coulomb_2c_ints.convert(coulomb_2c_ints_ptr);
    Eigen::Map<Eigen::MatrixXd> X_g_Y(coulomb_2c_ints_ptr, dfnbf, dfnbf);
    //----------------------------------------//
    // Contract C_mu_nu^X with (X|Y)
    for(int mu = 0; mu < obsnbf; ++mu){
      const double* C_mu = C->retrieve_pair_block(0, mu, ints_type_idx);
      //----------------------------------------//
      std::shared_ptr<Eigen::MatrixXd> gtilde_mu(new Eigen::MatrixXd(obsnbf, dfnbf));
      *gtilde_mu = Eigen::MatrixXd::Zero(obsnbf, dfnbf);
      //----------------------------------------//
      for(int nu = 0; nu < obsnbf; ++nu){
        const double* C_mu_nu = C_mu + nu*dfnbf;
        //----------------------------------------//
        // transfer one-by-one for now since the Eigen map doesn't seem to
        //   work like I think it does
        Eigen::VectorXd Cvect(dfnbf);
        for(int X = 0; X < dfnbf; ++X){
          Cvect(X) = C_mu_nu[X];
        }
        //----------------------------------------//
        gtilde_mu->row(nu) += X_g_Y * Cvect;
      }
      g_tilde[mu] = gtilde_mu;
      C->release_pair_block(0, mu, ints_type_idx);
    }
    //----------------------------------------//
    deallocate(coulomb_2c_ints_ptr);
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Get the coefficients (and transform the bra if we're doing integrals)            {{{1 */ #if 1 // begin fold
  //----------------------------------------------------------------------------//
  for(int mu = 0; mu < obsnbf; ++mu){
    const double* C_mu = C->retrieve_pair_block(0, mu, ints_type_idx);
    for(int nu = 0; nu < obsnbf; ++nu){
      const double* C_mu_nu = C_mu + nu*dfnbf;
      //----------------------------------------//
      std::shared_ptr<Eigen::MatrixXd> gtilde_mu(new Eigen::MatrixXd(obsnbf, dfnbf));
      *gtilde_mu = Eigen::MatrixXd::Zero(obsnbf, dfnbf);
      //----------------------------------------//
      // transfer one-by-one for now since the Eigen map doesn't seem to
      //   work like I think it does
      Eigen::VectorXd Cvect(dfnbf);
      Cvect = Eigen::VectorXd::Zero(dfnbf);
      for(int X = 0; X < dfnbf; ++X){
        Cvect(X) = C_mu_nu[X];
      }
      if(not do_integrals){
        ptree& coeff_parent = pt.add_child("coefficient_vector", ptree());
        coeff_parent.put("<xmlattr>.ao_index1", mu);
        coeff_parent.put("<xmlattr>.ao_index2", nu);
        writer.write_to_xml(Cvect, coeff_parent);
      }
      else{
        Eigen::MatrixXd g_mu_nu(obsnbf, obsnbf);
        g_mu_nu = Eigen::MatrixXd::Zero(obsnbf, obsnbf);
        for(int la = 0; la < obsnbf; ++la){
          for(int si = 0; si < obsnbf; ++si){
            g_mu_nu(la, si) += g_tilde[la]->row(si) * Cvect;
          }
        }
        //----------------------------------------//
        ptree& int_parent = pt.add_child("integral_block", ptree());
        int_parent.put("<xmlattr>.ao_index1", mu);
        int_parent.put("<xmlattr>.ao_index2", nu);
        writer.write_to_xml(g_mu_nu, int_parent);
      }
    }
    C->release_pair_block(0, mu, ints_type_idx);
  }
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
  /* Clean up                                             		                        {{{1 */ #if 1 // begin fold
  C->deactivate();
  /*****************************************************************************************/ #endif //1}}}
  /*=======================================================================================*/
}

std::pair<TwoBodyOperSet::type, Ref<IntParams> >
WavefunctionWorld::init_df_kernel(std::string kernel_key) {

  if (not (kernel_key == "coulomb" ||
           kernel_key == "delta" ||
           kernel_key.find("exp") != std::string::npos)
     )
    throw InputError("invalid value",
                           __FILE__,
                           __LINE__,
                           "df_kernel",
                           kernel_key.c_str());

  if (kernel_key == "coulomb") {
    return std::make_pair(TwoBodyOperSet::ERI, ParamsRegistry::instance()->value("") );
  }
  if (kernel_key == "delta") {
    return std::make_pair(TwoBodyOperSet::DeltaFunction, ParamsRegistry::instance()->value("") );
  }

  if (kernel_key.find("exp") != std::string::npos) {

    std::string::size_type s = kernel_key.find("exp");
    std::string exp_params;
    std::string kernel = kernel_key;
    pop_till_token(kernel, '(');
    if (not kernel.empty() && kernel.find(')') != std::string::npos) {
      std::string result = pop_till_token(kernel, ')');
      if (kernel.empty())
        exp_params = result;
    }

    if (exp_params.empty())
      throw InputError("improperly formatted exponential df kernel",
                       __FILE__,
                       __LINE__,
                       "df_kernel",
                       kernel_key.c_str());
    // exponential kernel_key must have positive lengthscale
    std::istringstream iss(exp_params);
    double lengthscale;
    iss >> lengthscale;
    double param;  iss >> param;
    if (lengthscale <= 0.0)
      throw InputError("exponential df kernel must have positive range",
                       __FILE__,
                       __LINE__,
                       "df_kernel",
                       kernel_key.c_str());
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

    const std::string params_key = ParamsRegistry::instance()->add(intparams);
    return std::make_pair(TwoBodyOperSet::R12_0_G12, ParamsRegistry::instance()->value(params_key));
  }

  // unreachable
  assert(false);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
