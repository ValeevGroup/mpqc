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

#ifdef __GNUG__
#pragma implementation
#endif

// includes go here
#include <chemistry/qc/mbptr12/wfnworld.h>
#include <chemistry/qc/basis/petite.h>
#include <util/class/scexception.h>
#include <chemistry/qc/df/df_runtime.h>

using namespace sc;

/*---------------
  WavefunctionWorld
 ---------------*/
static ClassDesc WavefunctionWorld_cd(
  typeid(WavefunctionWorld),"WavefunctionWorld",11,"virtual public SavableState",
  0, 0, create<WavefunctionWorld>);

WavefunctionWorld::WavefunctionWorld(
    const Ref<KeyVal>& keyval,
    Wavefunction* wfn) :
    wfn_(wfn)
{
  debug_ = 0;
  print_percent_ = 10.0;

  bs_df_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("df_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  if (bs_df_.nonnull()) {
    df_kernel_ = keyval->stringvalue("df_kernel", KeyValValuestring("coulomb"));
    if (df_kernel_ != "coulomb")
      throw InputError("invalid value",
                             __FILE__,
                             __LINE__,
                             "df_kernel",
                             df_kernel_.c_str(),
                             class_desc());

    df_solver_ = keyval->stringvalue("df_solver", KeyValValuestring("bunchkaufman_refine"));
    if (df_solver_ != "bunchkaufman" &&
        df_solver_ != "bunchkaufman_inv" &&
        df_solver_ != "bunchkaufman_refine" &&
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
#if HAVE_MPIIO
    ints_method_ = StoreMethod::mpi;
#else
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- store_ints=mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else if (ints_str == std::string("mem-mpi")) {
#if HAVE_MPIIO
    ints_method_ = StoreMethod::mem_mpi;
#else
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- store_ints=mem_mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else {
    throw std::runtime_error("WavefunctionWorld::WavefunctionWorld -- invalid value for keyword r12ints");
  }

  // Get the filename prefix to store the integrals
  const std::string ints_file_default = SCFormIO::fileext_to_filename(".moints");
  ints_file_ = keyval->stringvalue("ints_file",KeyValValuestring(ints_file_default));
  // if the last character of ints_file is '/' then append the default basename
  if (*(ints_file_.rbegin()) == '/')
    ints_file_ += ints_file_default;
  ExEnv::out0() << indent << "ints_file = " << ints_file_ << std::endl;

  // dynamic load balancing?
  dynamic_ = static_cast<bool>(keyval->booleanvalue("dynamic",KeyValValueboolean(false)));

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();
  integral()->set_basis(wfn_->basis());

  initialize();
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

  tfactory_ << SavableState::restore_state(si);
  moints_runtime_ << SavableState::restore_state(si);
  fockbuild_runtime_ << SavableState::restore_state(si);
  bs_df_ << SavableState::restore_state(si);
}

WavefunctionWorld::~WavefunctionWorld()
{
}

void WavefunctionWorld::save_data_state(StateOut& so)
{
  SavableState::save_state(wfn_,so);

  so.put((int)ints_method_);
  so.put(ints_file_);
  so.put(debug_);
  so.put(dynamic_);
  so.put(print_percent_);

  SavableState::save_state(tfactory_.pointer(),so);
  SavableState::save_state(moints_runtime_.pointer(),so);
  SavableState::save_state(fockbuild_runtime_.pointer(),so);
  SavableState::save_state(bs_df_.pointer(),so);
}

void
WavefunctionWorld::initialize()
{
  tfactory_ = new MOIntsTransformFactory(integral());
  tfactory_->set_dynamic(dynamic_);
  tfactory_->set_ints_method(ints_method_);
  tfactory_->set_file_prefix(ints_file_);

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
      dfparams = new DensityFittingParams(bs_df_, df_kernel_, df_solver_);
    }
    moints_runtime_ = new MOIntsRuntime(tfactory_, dfparams);

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

  if (bs_df_.nonnull()) {
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
#if HAVE_MPIIO
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
