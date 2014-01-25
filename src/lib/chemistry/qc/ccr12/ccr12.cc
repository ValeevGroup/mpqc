//
// ccr12.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: TS & EFV
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

#include <util/misc/scexception.h>
#include <chemistry/qc/ccr12/ccr12.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/scf/clscf.h>

using namespace std;
using namespace sc;


/*--------------------------------
  CCR12
 --------------------------------*/

static ClassDesc CCR12_cd(
  typeid(CCR12),"CCR12",1,"public Wavefunction",
  0,create<CCR12>,create<CCR12>);

CCR12::CCR12(StateIn& s): SavableState(s), Wavefunction(s), ccr12_info_(0) {
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructor not yet implemented",__FILE__,__LINE__);
}


CCR12::CCR12(const Ref<KeyVal>& keyval): Wavefunction(keyval), ccr12_info_(0) {

  reference_ << keyval->describedclassvalue("reference");
  if (reference_.null()) {
    ExEnv::err0() << "MBPT2::MBPT2: no reference wavefunction" << endl;
    abort();
  }
  copy_orthog_info(reference_);

  thrgrp_ = ThreadGrp::get_default_threadgrp();
  msggrp_ = MessageGrp::get_default_messagegrp();
  mem_ = MemoryGrp::get_default_memorygrp();
  timer_ = new RegionTimer();

  // if world not given, make this the center of a new World
  Ref<WavefunctionWorld> world; world << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world.null())
    world = new WavefunctionWorld(keyval);
  if (world.null())
    throw InputError("CCR12 requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world->wfn() == 0) world->set_wfn(this);

  const bool spin_restricted = false;   // do not use spin-restricted orbitals -> for ROHF use semicanonical orbitals
  Ref<OrbitalSpace> vbs;

  nfzc_ = keyval->intvalue("nfzc");
  string nfzc_charval = keyval->stringvalue("nfzc");
  if (nfzc_charval == "auto") {
    if (molecule()->max_z() > 30) {
      ExEnv::err0() << "CCR12: cannot use \"nfzc = auto\" for Z > 30" << endl;
      abort();
    }
    nfzc_ = molecule()->n_core_electrons()/2;
    ExEnv::out0() << "  CCR12: auto-freezing " << nfzc_ << " core orbitals" << endl;
  }
  nfzv_ = keyval->intvalue("nfzv");

  Ref<RefWavefunction> refinfo = new SD_RefWavefunction(world, ref(), spin_restricted,
                                                        nfzc_, nfzv_,
                                                        vbs);
  r12world_ = new R12WavefunctionWorld(keyval, refinfo);
  Ref<R12Technology> r12tech = r12world_->r12tech();

#if 0
  // CABS singles is not supported yet in CCR12.
  cabs_singles_ = keyval->booleanvalue("cabs_singles",KeyValValueboolean(true));
#endif

  this->set_desired_value_accuracy(desired_value_accuracy());

  ExEnv::out0() << endl << indent << "-------- CCR12 calculation --------"<< endl;
  if (sizeof(long) < 8)
    ExEnv::out0() << indent << "!!!!!!!!!! \"long int\" less than 8 bytes !!!!!!!!!!" << endl;

  ndiis_=keyval->intvalue("ndiis", KeyValValueint(2));
  diis_start_ = keyval->intvalue("diis_start", KeyValValueint(0));

  CLSCF* clscfref=dynamic_cast<CLSCF*>(ref().pointer());
  rhf_=(clscfref!=0);

  // maxiter
  maxiter_=keyval->intvalue("maxiter",   KeyValValueint(100));
  // cctresh
  ccthresh_=keyval->doublevalue("ccthresh", KeyValValuedouble(1.0e-9));
  // get the memory sizes
  memorysize_ = keyval->longvalue("memory",   KeyValValuelong(200000000));
  ExEnv::out0() << indent << "Memory size per node: " << memorysize_ << endl;
  worksize_ = keyval->longvalue("workmemory",   KeyValValuelong(50000000));
#ifdef DISK_BASED_SMITH
  worksize_ = memorysize_; 
#endif
  ExEnv::out0() << indent << "Work   size per node: " << worksize_   << endl;
  tilesize_forced_ = keyval->intvalue("force_tilesize",KeyValValueint(0));
}


CCR12::~CCR12(){
  mem_->sync();
  if (ccr12_info_ != 0) delete ccr12_info_;
}

void CCR12::compute(){

  reference_->set_desired_value_accuracy(desired_value_accuracy()/ref_to_ccr12_acc());

  r12world()->initialize();

  // CCR12_Info will do integral evaluation, before MemoryGrp is used by Tensors
  if (ccr12_info_ != 0) delete ccr12_info_;
  ccr12_info_=new CCR12_Info(r12world(),mem_,memorysize_,ref(),nfzc_,nfzv_,
                  molecule()->point_group()->char_table().nirrep(),worksize_,memorysize_,mem_->n(),ndiis_,
                  theory_,perturbative_,tilesize_forced_);

}


void CCR12::obsolete() {
  if (reference_) reference_->obsolete();
  Wavefunction::obsolete();
}

/// utilities >>>>>>> form here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void CCR12::print_theory() {
  if (theory_ == "notheory") throw InputError("CCR12::CCR12 -- no theory specified",__FILE__,__LINE__);
  ExEnv::out0() << endl << indent << "Theory:       " << theory_ << endl;
  if (perturbative_ != "")
    ExEnv::out0() << endl << indent << "Perturbative: " << perturbative_ << endl << endl;
}

void CCR12::print_iteration_header(string theory){
  if (mem_->me()==0) {
    ExEnv::out0() << endl << indent << theory << " iterations:" << endl;
    ExEnv::out0() << indent << "iter      corr energy        residual RMS        Wall " << endl;
    ExEnv::out0() << indent << "======================================================" << endl;
  }
}

void CCR12::print_iteration_header_short(string theory){
  if (mem_->me()==0) {
    ExEnv::out0() << endl << indent << theory << " iterations:" << endl;
    ExEnv::out0() << indent << "iter      residual RMS        Wall " << endl;
    ExEnv::out0() << indent << "===================================" << endl;
  }
}

void CCR12::print_iteration_footer(){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << "======================================================" << endl;
  }
}

void CCR12::print_iteration_footer_short(){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << "===================================" << endl;
  }
}

void CCR12::print_iteration(int iter,double energy,double rnorm,double iter_start,double iter_end){
  double iter_duration=iter_end-iter_start;
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << setw(4) << iter
                            << setw(20) << setprecision(13) << energy
                            << setw(20) << setprecision(13) << rnorm
                            << setw(10) << setprecision(2)  << iter_duration << endl;
  }
}

void CCR12::print_iteration_short(int iter,double rnorm,double iter_start,double iter_end){
  double iter_duration=iter_end-iter_start;
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << setw(4) << iter
                            << setw(20) << setprecision(13) << rnorm
                            << setw(10) << setprecision(2)  << iter_duration << endl;
  }
}

void CCR12::print_timing(double timing, std::string tag){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << "Elapsed time [ "<< tag << " ]: "
                            << setw(10) << setprecision(2)  << timing << endl << endl;
  }
}


void CCR12::print_correction(double corr, double base, string theory){
 if (mem_->me()==0) {
  ExEnv::out0() << endl;
  ExEnv::out0() << indent << theory << fixed << setprecision(10) << " correction: " << corr      << endl;
  ExEnv::out0() << indent << theory << fixed << setprecision(10) << " energy    : " << corr+base << endl;
  ExEnv::out0() << endl;
 }
}


void CCR12::print(ostream&o) const {
  o << indent << "CCR12:" << endl;
  o << incindent;
  Wavefunction::print(o);
  info()->print(o);
  o << decindent;
}

