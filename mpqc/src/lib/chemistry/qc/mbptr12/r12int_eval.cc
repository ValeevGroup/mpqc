//
// r12int_eval.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------------
  R12IntEval
 -----------------*/
static ClassDesc R12IntEval_cd(
  typeid(R12IntEval),"R12IntEval",1,"virtual public SavableState",
  0, 0, 0);

R12IntEval::R12IntEval(const Ref<R12IntEvalInfo>& r12info) :
  r12info_(r12info)
{
    int nocc_act = r12info_->nocc_act();
    dim_aa_ = new SCDimension((nocc_act*(nocc_act-1))/2);
    dim_ab_ = new SCDimension(nocc_act*nocc_act);
    dim_s_ = new SCDimension((nocc_act*(nocc_act+1))/2);
    dim_t_ = dim_aa_;

    Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
    Vaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Vab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    Xaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Xab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    Baa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Bab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    Aaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
    Aab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
    emp2pair_aa_ = local_matrix_kit->vector(dim_aa_);
    emp2pair_ab_ = local_matrix_kit->vector(dim_ab_);
    init_intermeds_();
    init_tforms_();
    
    // Default values
    gebc_ = true;
    abs_method_ = LinearR12::ABS_CABSPlus;
    stdapprox_ = LinearR12::StdApprox_Ap;
    spinadapted_ = true;
    evaluated_ = false;
    debug_ = 0;
}

R12IntEval::R12IntEval(StateIn& si) : SavableState(si)
{
  r12info_ << SavableState::restore_state(si);
  dim_aa_ << SavableState::restore_state(si);
  dim_ab_ << SavableState::restore_state(si);
  dim_s_ << SavableState::restore_state(si);
  dim_t_ << SavableState::restore_state(si);

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  Vaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Vab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Xaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Xab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Baa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Bab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  Aaa_ = local_matrix_kit->matrix(dim_aa_,dim_aa_);
  Aab_ = local_matrix_kit->matrix(dim_ab_,dim_ab_);
  emp2pair_aa_ = local_matrix_kit->vector(dim_aa_);
  emp2pair_ab_ = local_matrix_kit->vector(dim_ab_);

  Vaa_.restore(si);
  Vab_.restore(si);
  Xaa_.restore(si);
  Xab_.restore(si);
  Baa_.restore(si);
  Bab_.restore(si);
  emp2pair_aa_.restore(si);
  emp2pair_ab_.restore(si);

  ipjq_tform_ << SavableState::restore_state(si);
  ikjy_tform_ << SavableState::restore_state(si);
  init_tforms_();

  int gebc; si.get(gebc); gebc_ = (bool) gebc;
  int absmethod; si.get(absmethod); abs_method_ = (LinearR12::ABSMethod) absmethod;
  int stdapprox; si.get(stdapprox); stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  int spinadapted; si.get(spinadapted); spinadapted_ = (bool) spinadapted;
  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
  si.get(debug_);
}

R12IntEval::~R12IntEval()
{
  r12info_ = 0;
  dim_aa_ = 0;
  dim_ab_ = 0;
  dim_s_ = 0;
  dim_t_ = 0;
  ipjq_tform_ = 0;
  ikjy_tform_ = 0;
}

void
R12IntEval::save_data_state(StateOut& so)
{
  SavableState::save_state(r12info_.pointer(),so);
  SavableState::save_state(dim_aa_.pointer(),so);
  SavableState::save_state(dim_ab_.pointer(),so);
  SavableState::save_state(dim_s_.pointer(),so);
  SavableState::save_state(dim_t_.pointer(),so);

  Vaa_.save(so);
  Vab_.save(so);
  Xaa_.save(so);
  Xab_.save(so);
  Baa_.save(so);
  Bab_.save(so);
  Aaa_.save(so);
  Aab_.save(so);
  emp2pair_aa_.save(so);
  emp2pair_ab_.save(so);

  SavableState::save_state(ipjq_tform_.pointer(),so);
  SavableState::save_state(ikjy_tform_.pointer(),so);

  so.put((int)gebc_);
  so.put((int)abs_method_);
  so.put((int)stdapprox_);
  so.put((int)spinadapted_);
  so.put((int)evaluated_);
  so.put(debug_);
}

void
R12IntEval::obsolete()
{
  evaluated_ = false;
  ipjq_tform_->obsolete();
  ikjy_tform_->obsolete();
  init_intermeds_();
}

void R12IntEval::set_gebc(bool gebc) { gebc_ = gebc; };
void R12IntEval::set_absmethod(LinearR12::ABSMethod abs_method) { abs_method_ = abs_method; };
void R12IntEval::set_stdapprox(LinearR12::StandardApproximation stdapprox) { stdapprox_ = stdapprox; };
void R12IntEval::set_spinadapted(bool spinadapted) { spinadapted_ = spinadapted; };
void R12IntEval::set_debug(int debug) { if (debug >= 0) { debug_ = debug; r12info_->set_debug_level(debug_); }};
void R12IntEval::set_dynamic(bool dynamic) { r12info_->set_dynamic(dynamic); };
void R12IntEval::set_print_percent(double pp) { r12info_->set_print_percent(pp); };
void R12IntEval::set_memory(size_t nbytes) { r12info_->set_memory(nbytes); };

Ref<R12IntEvalInfo> R12IntEval::r12info() const { return r12info_; };
RefSCDimension R12IntEval::dim_aa() const { return dim_aa_; };
RefSCDimension R12IntEval::dim_ab() const { return dim_ab_; };
RefSCDimension R12IntEval::dim_s() const { return dim_s_; };
RefSCDimension R12IntEval::dim_t() const { return dim_t_; };

RefSCMatrix R12IntEval::V_aa() {
  compute();
  return Vaa_;
}

RefSCMatrix R12IntEval::X_aa() {
  compute();
  return Xaa_;
}

RefSCMatrix R12IntEval::B_aa() {
  compute();
  return Baa_;
}

RefSCMatrix R12IntEval::A_aa() {
  compute();
  return Aaa_;
}

RefSCMatrix R12IntEval::V_ab() {
  compute();
  return Vab_;
}

RefSCMatrix R12IntEval::X_ab() {
  compute();
  return Xab_;
}

RefSCMatrix R12IntEval::B_ab() {
  compute();
  return Bab_;
}

RefSCMatrix R12IntEval::A_ab() {
  compute();
  return Aab_;
}

RefSCVector R12IntEval::emp2_aa() {
  compute();
  return emp2pair_aa_;
}

RefSCVector R12IntEval::emp2_ab() {
  compute();
  return emp2pair_ab_;
}

RefDiagSCMatrix R12IntEval::evals() const { return r12info_->evals(); };

void
R12IntEval::checkpoint_() const
{
  int me = r12info_->msg()->me();
  Wavefunction* wfn = r12info_->wfn();

  if (me == 0 && wfn->if_to_checkpoint()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }
}

void
R12IntEval::init_tforms_()
{
  Ref<MOIntsTransformFactory> tfactory = r12info_->tfactory();

  if (ipjq_tform_.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->obs_space(),
                         r12info_->act_occ_space(),r12info_->obs_space());
    ipjq_tform_ = tfactory->twobody_transform("(ip|jq)");
  }

  if (ikjy_tform_.null()) {
    tfactory->set_spaces(r12info_->act_occ_space(),r12info_->occ_space(),
                         r12info_->act_occ_space(),r12info_->ribs_space());
    ikjy_tform_ = tfactory->twobody_transform("(ik|jy)");
  }
}

void
R12IntEval::init_intermeds_()
{
  Vaa_->unit();
  Vab_->unit();
  Baa_->unit();
  Bab_->unit();
  Aaa_.assign(0.0);
  Aab_.assign(0.0);

  Xaa_.assign(0.0);
  Xab_.assign(0.0);
  r2_contrib_to_X_();

  emp2pair_aa_.assign(0.0);
  emp2pair_ab_.assign(0.0);
}

void
R12IntEval::r2_contrib_to_X_()
{
  /*-----------------------------------------------
    Compute dipole and quadrupole moment integrals
   -----------------------------------------------*/
  RefSymmSCMatrix MX, MY, MZ, MXX, MYY, MZZ;
  r12info_->compute_multipole_ints(MX,MY,MZ,MXX,MYY,MZZ);
  if (debug_)
    ExEnv::out0() << indent << "Computed multipole moment integrals" << endl;
    
  const int nproc = r12info_->msg()->n();
  const int me = r12info_->msg()->me();

  MOPairIter_SD ij_iter(r12info_->act_occ_space());
  MOPairIter_SD kl_iter(r12info_->act_occ_space());

  for(kl_iter.start();int(kl_iter);kl_iter.next()) {

    const int kl = kl_iter.ij();
    int kl_proc = kl%nproc;
    if (kl_proc != me)
      continue;
    const int k = kl_iter.i();
    const int l = kl_iter.j();
    const int kl_aa = kl_iter.ij_aa();
    const int kl_ab = kl_iter.ij_ab();
    const int lk_ab = kl_iter.ji_ab();

    for(ij_iter.start();int(ij_iter);ij_iter.next()) {

      const int i = ij_iter.i();
      const int j = ij_iter.j();
      const int ij_aa = ij_iter.ij_aa();
      const int ij_ab = ij_iter.ij_ab();
      const int ji_ab = ij_iter.ji_ab();

      /*----------------------------------
        Compute (r12)^2 contribution to X
       ----------------------------------*/
      double r1r1_ik = -1.0*(MXX->get_element(i,k) + MYY->get_element(i,k) + MZZ->get_element(i,k));
      double r1r1_il = -1.0*(MXX->get_element(i,l) + MYY->get_element(i,l) + MZZ->get_element(i,l));
      double r1r1_jk = -1.0*(MXX->get_element(j,k) + MYY->get_element(j,k) + MZZ->get_element(j,k));
      double r1r1_jl = -1.0*(MXX->get_element(j,l) + MYY->get_element(j,l) + MZZ->get_element(j,l));
      double r1r2_ijkl = MX->get_element(i,k)*MX->get_element(j,l) +
        MY->get_element(i,k)*MY->get_element(j,l) +
        MZ->get_element(i,k)*MZ->get_element(j,l);
      double r1r2_ijlk = MX->get_element(i,l)*MX->get_element(j,k) +
        MY->get_element(i,l)*MY->get_element(j,k) +
        MZ->get_element(i,l)*MZ->get_element(j,k);
      double delta_ik = (i==k ? 1.0 : 0.0);
      double delta_il = (i==l ? 1.0 : 0.0);
      double delta_jk = (j==k ? 1.0 : 0.0);
      double delta_jl = (j==l ? 1.0 : 0.0);
      
      double Xab_ijkl = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
      Xab_.accumulate_element(ij_ab,kl_ab,Xab_ijkl);

      if (ij_ab != ji_ab) {
        double Xab_jikl = r1r1_jk * delta_il + r1r1_il * delta_jk - 2.0*r1r2_ijlk;
        Xab_.accumulate_element(ji_ab,kl_ab,Xab_jikl);
      }

      if (kl_ab != lk_ab) {
        double Xab_ijlk = r1r1_il * delta_jk + r1r1_jk * delta_il - 2.0*r1r2_ijlk;
        Xab_.accumulate_element(ij_ab,lk_ab,Xab_ijlk);
      }

      if (ij_ab != ji_ab && kl_ab != lk_ab) {
        double Xab_jilk = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
        Xab_.accumulate_element(ji_ab,lk_ab,Xab_jilk);
      }

      if (ij_aa != -1 && kl_aa != -1) {
        double Xaa_ijkl = r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl -
          r1r1_jk * delta_il - r1r1_il * delta_jk + 2.0*r1r2_ijlk;
        Xaa_.accumulate_element(ij_aa,kl_aa,Xaa_ijkl);
      }

    }
  }          
}

const int
R12IntEval::tasks_with_ints_(const Ref<R12IntsAcc> ints_acc, vector<int>& map_to_twi)
{
  int nproc = r12info_->msg()->n();
  
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) nproc_with_ints++;
 
  map_to_twi.resize(nproc);
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (ints_acc->has_access(proc)) {
      map_to_twi[proc] = count;
      count++;
    }
      else
        map_to_twi[proc] = -1;

  ExEnv::out0() << indent << "Computing intermediates on " << nproc_with_ints
    << " processors" << endl;
    
  return nproc_with_ints;
}


void
R12IntEval::compute()
{
  obs_contrib_to_VXB_gebc_();
  if (r12info_->basis() != r12info_->basis_ri())
    abs1_contrib_to_VXB_gebc_();
  evaluated_ = true;
}

void
R12IntEval::globally_sum_intermeds_()
{
  Ref<MessageGrp> msg = r12info_->msg();
  // If there's only one task then there's nothing to do
  if (msg->n() == 1)
    return;

  const int naa = dim_aa_.n();
  const int nab = dim_ab_.n();
  
  double *Vaa_ijkl = new double[naa*naa];
  double *Baa_ijkl = new double[naa*naa];
  double *Xaa_ijkl = new double[naa*naa];
  double *Aaa_ijkl = new double[naa*naa];
  double *Vab_ijkl = new double[nab*nab];
  double *Bab_ijkl = new double[nab*nab];
  double *Xab_ijkl = new double[nab*nab];
  double *Aab_ijkl = new double[nab*nab];
  double *emp2_aa = new double[naa];
  double *emp2_ab = new double[nab];

  Vaa_.convert(Vaa_ijkl);
  Xaa_.convert(Xaa_ijkl);
  Baa_.convert(Baa_ijkl);
  Aaa_.convert(Aaa_ijkl);
  Vab_.convert(Vab_ijkl);
  Xab_.convert(Xab_ijkl);
  Bab_.convert(Bab_ijkl);
  Aab_.convert(Aab_ijkl);
  emp2pair_aa_.convert(emp2_aa);
  emp2pair_ab_.convert(emp2_ab);

  msg->sum(Vaa_ijkl,naa*naa,0,-1);
  msg->sum(Vab_ijkl,nab*nab,0,-1);
  msg->sum(Xaa_ijkl,naa*naa,0,-1);
  msg->sum(Xab_ijkl,nab*nab,0,-1);
  msg->sum(Baa_ijkl,naa*naa,0,-1);
  msg->sum(Bab_ijkl,nab*nab,0,-1);
  msg->sum(Aaa_ijkl,naa*naa,0,-1);
  msg->sum(Aab_ijkl,nab*nab,0,-1);
  msg->sum(emp2_aa,naa,0,-1);
  msg->sum(emp2_aa,nab,0,-1);

  Vaa_.assign(Vaa_ijkl);
  Xaa_.assign(Xaa_ijkl);
  Baa_.assign(Baa_ijkl);
  Aaa_.assign(Aaa_ijkl);
  Vab_.assign(Vab_ijkl);
  Xab_.assign(Xab_ijkl);
  Bab_.assign(Bab_ijkl);
  Aab_.assign(Aab_ijkl);
  emp2pair_aa_.assign(emp2_aa);
  emp2pair_ab_.assign(emp2_ab);

  delete[] Vaa_ijkl;
  delete[] Xaa_ijkl;
  delete[] Baa_ijkl;
  delete[] Aaa_ijkl;
  delete[] Vab_ijkl;
  delete[] Xab_ijkl;
  delete[] Bab_ijkl;
  delete[] Aab_ijkl;
  delete[] emp2_aa;
  delete[] emp2_ab;
  
  if (debug_)
    ExEnv::out0() << indent << "Collected contributions to the intermediates from all tasks" << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
