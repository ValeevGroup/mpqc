//
// vxb_eval_info.cc
//
// Copyright (C) 2003 Edward Valeev
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
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/singlerefinfo.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  R12IntEvalInfo
 ---------------*/
static ClassDesc R12IntEvalInfo_cd(
  typeid(R12IntEvalInfo),"R12IntEvalInfo",9,"virtual public SavableState",
  0, 0, create<R12IntEvalInfo>);

R12IntEvalInfo::R12IntEvalInfo(
    const Ref<KeyVal>& keyval,
    Wavefunction* wfn,
    const Ref<SCF>& ref,
    unsigned int nfzc,
    unsigned int nfzv,
    bool spinadapted,
    bool delayed_initialization
    ) :
    wfn_(wfn), spinadapted_(spinadapted),
    initialized_(false),
    refinfo_(new SingleRefInfo(ref,nfzc,nfzv,delayed_initialization))
{
  // Default values
  memory_ = DEFAULT_SC_MEMORY;
  debug_ = 0;
  print_percent_ = 10.0;

  bs_aux_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("aux_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  if (bs_aux_.pointer() == NULL)
      bs_aux_ = ref->basis();

  bs_vir_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("vir_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  if (bs_vir_.pointer() == NULL)
      bs_vir_ = ref->basis();

  // Determine how to store MO integrals
  std::string ints_str = keyval->stringvalue("store_ints",KeyValValuestring("posix"));
  if (ints_str == std::string("mem")) {
    ints_method_ = MOIntsTransformFactory::StoreMethod::mem_only;
  }
  else if (ints_str == std::string("posix")) {
    ints_method_ = MOIntsTransformFactory::StoreMethod::posix;
  }
  else if (ints_str == std::string("mem-posix")) {
    ints_method_ = MOIntsTransformFactory::StoreMethod::mem_posix;
  }
  else if (ints_str == std::string("mpi")) {
#if HAVE_MPIIO
    ints_method_ = MOIntsTransformFactory::StoreMethod::mpi;
#else
    throw std::runtime_error("R12IntEvalInfo::R12IntEvalInfo -- store_ints=mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else if (ints_str == std::string("mem-mpi")) {
#if HAVE_MPIIO
    ints_method_ = MOIntsTransformFactory::StoreMethod::mem_mpi;
#else
    throw std::runtime_error("R12IntEvalInfo::R12IntEvalInfo -- store_ints=mem_mpi is not valid in this environment (no MPI-I/O detected)");
#endif
  }
  else {
    throw std::runtime_error("R12IntEvalInfo::R12IntEvalInfo -- invalid value for keyword r12ints");
  }

  // Get the filename prefix to store the integrals
  const std::string ints_file_default = SCFormIO::fileext_to_filename(".moints");
  ints_file_ = keyval->stringvalue("ints_file",KeyValValuestring(ints_file_default));
  // if the last character of ints_file is '/' then append the default basename
  if (*(ints_file_.rbegin()) == '/')
    ints_file_ += ints_file_default;
  ExEnv::out0() << indent << "ints_file = " << ints_file_ << endl;

  // dynamic load balancing?
  dynamic_ = static_cast<bool>(keyval->booleanvalue("dynamic",KeyValValueboolean(0)));

  r12tech_ = new R12Technology(keyval,ref->basis(),bs_vir_,bs_aux_);
  // Make sure can use the integral factory for R12 calcs
  r12tech_->check_integral_factory(integral());

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  integral()->set_basis(ref->basis());
  Ref<PetiteList> plist = integral()->petite_list();
  RefSCDimension oso_dim = plist->SO_basisdim();

  if (!delayed_initialization)
      initialize();
}

R12IntEvalInfo::R12IntEvalInfo(StateIn& si) : SavableState(si)
{
  Ref<Wavefunction> wfn; wfn << SavableState::restore_state(si);
  wfn_ = wfn.pointer();
  bs_aux_ << SavableState::restore_state(si);
  bs_vir_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  int spinadapted; si.get(spinadapted); spinadapted_ = (bool)spinadapted;
  int ints_method; si.get(ints_method);
  ints_method_ = static_cast<R12IntEvalInfo::StoreMethod::type>(ints_method);
  si.get(ints_file_);

  double memory; si.get(memory); memory_ = (size_t) memory;
  si.get(debug_);
  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 2) {
    si.get(print_percent_);
  }

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 4) {
    abs_space_ << SavableState::restore_state(si);
    ribs_space_ << SavableState::restore_state(si);
    vir_act_ << SavableState::restore_state(si);
    vir_ << SavableState::restore_state(si);
    vir_sb_ << SavableState::restore_state(si);
    tfactory_ << SavableState::restore_state(si);
  }

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 5) {
    refinfo_ << SavableState::restore_state(si);
  }

  int initialized; si.get(initialized); initialized_ = (bool)initialized;
}

R12IntEvalInfo::~R12IntEvalInfo()
{
}

void R12IntEvalInfo::save_data_state(StateOut& so)
{
  SavableState::save_state(wfn_,so);
  SavableState::save_state(bs_aux_.pointer(),so);
  SavableState::save_state(bs_vir_.pointer(),so);
  SavableState::save_state(bs_ri_.pointer(),so);

  so.put((int)spinadapted_);
  so.put((int)ints_method_);
  so.put(ints_file_);
  so.put((double)memory_);
  so.put(debug_);
  so.put((int)dynamic_);
  so.put(print_percent_);

  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
  SavableState::save_state(vir_act_.pointer(),so);
  SavableState::save_state(vir_.pointer(),so);
  SavableState::save_state(vir_sb_.pointer(),so);
  SavableState::save_state(tfactory_.pointer(),so);
  SavableState::save_state(refinfo_.pointer(),so);
  so.put(initialized_);
}

void
R12IntEvalInfo::initialize()
{
  if (!initialized_) {
      refinfo_->initialize();
      construct_ri_basis_(safety_check());
      construct_orthog_vir_();

      tfactory_ = new MOIntsTransformFactory(integral(),refinfo()->orbs(Alpha));
      tfactory_->set_memory(memory_);
      tfactory_->set_dynamic(dynamic_);
      tfactory_->set_ints_method(ints_method_);
      tfactory_->set_file_prefix(ints_file_);
      initialized_ = true;

      // provide hints to the factory about the likely use of transforms
      // 1) if stdapprox is A' or A'' most transforms will never be reused
      // 2) otherwise use persistent transforms
      if (stdapprox() == LinearR12::StdApprox_Ap || stdapprox() == LinearR12::StdApprox_App)
        tfactory_->hints().data_persistent(false);
      else
        tfactory_->hints().data_persistent(true);
  }
}

const Ref<MOIndexSpace>&
R12IntEvalInfo::ribs_space(const SpinCase1& S) const
{
    if (abs_method() == LinearR12::ABS_CABS || abs_method() == LinearR12::ABS_CABSPlus)
	return vir_spaces_[S].ri_;
    else
	throw ProgrammingError("CABS space requested by abs_method set to ABS/ABS+",__FILE__,__LINE__);
}

const Ref<SingleRefInfo>&
R12IntEvalInfo::refinfo() const {
  return refinfo_;
}

const std::string& R12IntEvalInfo::ints_file() const
{
  return ints_file_;
}

void
R12IntEvalInfo::set_memory(const size_t memory)
{
  if (memory > 0)
    memory_ = memory;
  tfactory_->set_memory(memory_);
}

void
R12IntEvalInfo::throw_if_spin_polarized() const
{
  if (refinfo()->spin_polarized())
    throw ProgrammingError("R12IntEvalInfo -- spin-independent space is requested but the reference function is spin-polarized",__FILE__,__LINE__);
}

void
R12IntEvalInfo::SpinSpaces::init(const Ref<SingleRefInfo>& refinfo, const SpinCase1& spincase)
{
  vir_sb_ = refinfo->uocc_sb(spincase);
  vir_ = refinfo->uocc(spincase);
  vir_act_ = refinfo->uocc_act(spincase);
}

void
R12IntEvalInfo::SpinSpaces::init(const Ref<SingleRefInfo>& refinfo,
                                 const SpinCase1& spincase,
                                 const Ref<MOIndexSpace>& vbs)
{
  std::string id, label;
  if (spincase == Alpha)
    id = "E(sym)";
  else
    id = "e(sym)";
  vir_sb_ = R12IntEvalInfo::orthog_comp(refinfo->occ_sb(spincase), vbs, "e(sym)", "VBS", refinfo->ref()->lindep_tol());
  // Design flaw!!! Need to compute Fock matrix right here but can't since Fock is built into R12IntEval
  // Need to move all relevant code outside of MBPT2-F12 code
  vir_ = vir_sb_;
  vir_act_ = vir_;
}

void
R12IntEvalInfo::construct_orthog_vir_()
{
  if (vir_.nonnull())
    return;

  const bool spin_polarized = refinfo()->spin_polarized();

  Ref<GaussianBasisSet> obs = refinfo()->ref()->basis();
  if (obs == bs_vir_) {
    if (!spin_polarized) {
      vir_ = refinfo()->uocc();
      vir_sb_ = refinfo()->uocc_sb();
      if (refinfo()->nfzv() == 0)
        vir_act_ = vir_;
      else
        vir_act_ = refinfo()->uocc_act();
    }
    vir_spaces_[Alpha].init(refinfo(),Alpha);
    vir_spaces_[Beta].init(refinfo(),Beta);
    nlindep_vir_ = 0;
  }
  else {
    if (refinfo()->nfzv() != 0)
      throw std::runtime_error("R12IntEvalInfo::construct_orthog_vir_() -- nfzv_ != 0 is not allowed yet");

    // This is a set of orthonormal functions that span VBS
    Ref<MOIndexSpace> vir_space = orthogonalize("e","VBS", bs_vir_, integral(), refinfo()->ref()->orthog_method(), refinfo()->ref()->lindep_tol(), nlindep_vir_);
    if (!spin_polarized) {
      // Now project out occupied MOs
      vir_sb_ = orthog_comp(refinfo()->docc_sb(), vir_space, "e(sym)", "VBS", refinfo()->ref()->lindep_tol());
      // Design flaw!!! Need to compute Fock matrix right here but can't since Fock is built into R12IntEval
      // Need to move all relevant code outside of MBPT2-F12 code
      if (refinfo()->nfzv() != 0)
        throw std::runtime_error("R12IntEvalInfo::construct_orthog_vir_() -- nfzv_ != 0 is not allowed yet");
      vir_ = vir_sb_;
      vir_act_ = vir_sb_;
    }
    vir_spaces_[Alpha].init(refinfo(),Alpha,vir_space);
    vir_spaces_[Beta].init(refinfo(),Beta,vir_space);
  }
}

void
R12IntEvalInfo::vir(const SpinCase1& S, const Ref<MOIndexSpace>& sp)
{
  if (refinfo()->spin_polarized()) {
    vir_spaces_[S].vir_ = sp;
  }
  else {
    vir_spaces_[Alpha].vir_ = sp;
    vir_spaces_[Beta].vir_ = sp;
    vir_ = sp;
  }
}

void
R12IntEvalInfo::vir_sb(const SpinCase1& S, const Ref<MOIndexSpace>& sp)
{
  if (refinfo()->spin_polarized()) {
    vir_spaces_[S].vir_sb_ = sp;
  }
  else {
    vir_spaces_[Alpha].vir_sb_ = sp;
    vir_spaces_[Beta].vir_sb_ = sp;
    vir_sb_ = sp;
  }
}

void
R12IntEvalInfo::vir_act(const SpinCase1& S, const Ref<MOIndexSpace>& sp)
{
  if (refinfo()->spin_polarized()) {
    vir_spaces_[S].vir_act_ = sp;
  }
  else {
    vir_spaces_[Alpha].vir_act_ = sp;
    vir_spaces_[Beta].vir_act_ = sp;
    vir_act_ = sp;
  }
}

void
R12IntEvalInfo::print(std::ostream& o) const {

  o << indent << "R12IntEvalInfo:" << endl;
  o << incindent;

  if (bs_vir_->equiv(basis())) {
      o << indent << "Virtuals Basis Set (VBS):" << endl;
      o << incindent; bs_vir_->print(o); o << decindent << endl;
  }
  if (!bs_aux_->equiv(basis())) {
      o << indent << "Auxiliary Basis Set (ABS):" << endl;
      o << incindent; bs_aux_->print(o); o << decindent << endl;
  }

  r12tech()->print(o);

  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;

  std::string ints_str;
  switch (ints_method_) {
  case R12IntEvalInfo::StoreMethod::mem_only:
    ints_str = std::string("mem"); break;
  case R12IntEvalInfo::StoreMethod::mem_posix:
    ints_str = std::string("mem-posix"); break;
  case R12IntEvalInfo::StoreMethod::posix:
    ints_str = std::string("posix"); break;
#if HAVE_MPIIO
  case R12IntEvalInfo::StoreMethod::mem_mpi:
    ints_str = std::string("mem-mpi"); break;
  case R12IntEvalInfo::StoreMethod::mpi:
    ints_str = std::string("mpi"); break;
#endif
  default:
    throw std::runtime_error("R12IntEvalInfo::print -- invalid value of ints_method_");
  }
  o << indent << "How to Store Transformed Integrals: " << ints_str << endl << endl;
  o << indent << "Transformed Integrals file suffix: " << ints_file_ << endl << endl;
}

bool
R12IntEvalInfo::bc() const
{
  // 1) false if VBS != OBS
  const bool vbs_neq_obs = !bs_vir_->equiv(basis());
  // 2) or using ROHF reference
  const bool open_shell_but_not_uhf = refinfo()->ref()->spin_polarized() && !refinfo()->ref()->spin_unrestricted();
  return !(vbs_neq_obs || open_shell_but_not_uhf);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
