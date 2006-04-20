//
// vxb_eval_info.cc
//
// Copyright (C) 2003 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
  typeid(R12IntEvalInfo),"R12IntEvalInfo",8,"virtual public SavableState",
  0, 0, create<R12IntEvalInfo>);

R12IntEvalInfo::R12IntEvalInfo(MBPT2_R12* mbptr12)
{
  // Default values
  memory_ = DEFAULT_SC_MEMORY;
  debug_ = 0;
  dynamic_ = false;
  print_percent_ = 10.0;

  wfn_ = mbptr12;
  refinfo_ = new SingleRefInfo(mbptr12->ref(),mbptr12->nfzcore(),mbptr12->nfzvirt());
  bs_aux_ = mbptr12->aux_basis();
  bs_vir_ = mbptr12->vir_basis();

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  integral()->set_basis(basis());
  Ref<PetiteList> plist = integral()->petite_list();
  RefSCDimension oso_dim = plist->SO_basisdim();

  ints_method_ = mbptr12->r12ints_method();
  ints_file_ = mbptr12->r12ints_file();

  corrfactor_ = mbptr12->corrfactor();
  abs_method_ = mbptr12->abs_method();
  stdapprox_ = mbptr12->stdapprox();
  ansatz_ = mbptr12->ansatz();
  ebc_ = mbptr12->ebc();
  gbc_ = mbptr12->gbc();
  ks_ebcfree_ = mbptr12->ks_ebcfree();
  omit_P_ = mbptr12->omit_P();
  maxnabs_ = mbptr12->maxnabs();
  include_mp1_ = mbptr12->include_mp1();
  construct_ri_basis_(false);
  construct_orthog_vir_();

  tfactory_ = new MOIntsTransformFactory(integral(),refinfo()->orbs(Alpha));
  tfactory_->set_memory(memory_);
  tfactory_->set_ints_method(ints_method_);
  tfactory_->set_file_prefix(ints_file_);
}

R12IntEvalInfo::R12IntEvalInfo(StateIn& si) : SavableState(si)
{
  wfn_ = require_dynamic_cast<Wavefunction*>(SavableState::restore_state(si),
					     "R12IntEvalInfo::R12IntEvalInfo");

  bs_aux_ << SavableState::restore_state(si);
  bs_vir_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  int ints_method; si.get(ints_method); ints_method_ = static_cast<StoreMethod::type>(ints_method);
  si.get(ints_file_);

  double memory; si.get(memory); memory_ = (size_t) memory;
  si.get(debug_);
  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 2) {
    si.get(print_percent_);
  }

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 3) {
    int absmethod; si.get(absmethod); abs_method_ = (LinearR12::ABSMethod) absmethod;
  }
  
  if (si.version(::class_desc<R12IntEvalInfo>()) >= 7) {
    int gbc; si.get(gbc); gbc_ = (bool) gbc;
    int ebc; si.get(ebc); ebc_ = (bool) ebc;
    int stdapprox; si.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
    if (si.version(::class_desc<R12IntEvalInfo>()) >= 8) {
      int ansatz; si.get(ansatz); ansatz_ = (LinearR12::Ansatz) ansatz;
    }
    int spinadapted; si.get(spinadapted); spinadapted_ = (bool) spinadapted;
    int includemp1; si.get(includemp1); include_mp1_ = (bool) includemp1;
    int ks_ebcfree; si.get(ks_ebcfree); ks_ebcfree_ = static_cast<bool>(ks_ebcfree);
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
  
  maxnabs_ = 2;
  if (si.version(::class_desc<R12IntEvalInfo>()) >= 6) {
    si.get(maxnabs_);
  }
  
  if (si.version(::class_desc<R12IntEvalInfo>()) >= 7) {
    int ks_ebcfree; si.get(ks_ebcfree); ks_ebcfree_ = static_cast<bool>(ks_ebcfree);
    int omitP; si.get(omitP); omit_P_ = static_cast<bool>(omitP);
  }
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

  so.put((int)ints_method_);
  so.put(ints_file_);

  so.put((double)memory_);
  so.put(debug_);
  so.put((int)dynamic_);
  so.put(print_percent_);
  so.put((int)abs_method_);
  so.put((int)gbc_);
  so.put((int)ebc_);
  so.put((int)stdapprox_);
  so.put((int)ansatz_);
  so.put((int)spinadapted_);
  so.put((int)include_mp1_);
  
  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
  SavableState::save_state(vir_act_.pointer(),so);
  SavableState::save_state(vir_.pointer(),so);
  SavableState::save_state(vir_sb_.pointer(),so);
  SavableState::save_state(tfactory_.pointer(),so);
  SavableState::save_state(refinfo_.pointer(),so);
  
  so.put(maxnabs_);
  so.put((int)ks_ebcfree_);
  so.put((int)omit_P_);
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
R12IntEvalInfo::set_absmethod(LinearR12::ABSMethod abs_method)
{
  if (abs_method != abs_method_) {
    abs_method = abs_method_;
    construct_ri_basis_(false);
  }
}

void
R12IntEvalInfo::throw_if_spin_polarized() const
{
  if (refinfo()->ref()->spin_polarized())
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
  
  const bool spin_polarized = refinfo()->ref()->spin_polarized();

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
  if (refinfo()->ref()->spin_polarized()) {
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
  if (refinfo()->ref()->spin_polarized()) {
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
  if (refinfo()->ref()->spin_polarized()) {
    vir_spaces_[S].vir_act_ = sp;
  }
  else {
    vir_spaces_[Alpha].vir_act_ = sp;
    vir_spaces_[Beta].vir_act_ = sp;
    vir_act_ = sp;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
