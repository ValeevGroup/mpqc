//
// mp2r12_energy.cc
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
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-------------
  MP2R12Energy
 -------------*/
static ClassDesc MP2R12Energy_cd(
  typeid(MP2R12Energy),"MP2R12Energy",1,"virtual public SavableState",
  0, 0, create<MP2R12Energy>);

MP2R12Energy::MP2R12Energy(Ref<R12IntEval>& r12eval, LinearR12::StandardApproximation stdapp, int debug)
{
  r12eval_ = r12eval;
  stdapprox_ = stdapp;
  if (debug >= 0)
    debug_ = debug;
  else
    debug_ = 0;
  evaluated_ = false;
}

MP2R12Energy::MP2R12Energy(StateIn& si) : SavableState(si)
{
  r12eval_ << SavableState::restore_state(si);

  er12_aa_ << SavableState::restore_state(si);
  er12_ab_ << SavableState::restore_state(si);
  emp2r12_ab_ << SavableState::restore_state(si);
  emp2r12_ab_ << SavableState::restore_state(si);
  
  int stdapprox;
  si.get(stdapprox);
  stdapprox_ = (LinearR12::StandardApproximation) stdapprox;
  si.get(debug_);
  int evaluated;
  si.get(evaluated);
  evaluated_ = (bool) evaluated;
}

MP2R12Energy::~MP2R12Energy()
{
  r12eval_ = 0;
}

void MP2R12Energy::save_data_state(StateOut& so)
{
  SavableState::save_state(r12eval_.pointer(),so);

  er12_aa_.save(so);
  er12_ab_.save(so);
  emp2r12_aa_.save(so);
  emp2r12_ab_.save(so);
  
  so.put((int)stdapprox_);
  so.put(debug_);
  so.put((int)evaluated_);
}

void MP2R12Energy::obsolete()
{
  evaluated_ = false;
}

Ref<R12IntEval> MP2R12Energy::r12eval() const { return r12eval_; };
LinearR12::StandardApproximation MP2R12Energy::stdapp() const { return stdapprox_; };
void MP2R12Energy::set_debug(int debug) { debug_ = debug; };
int MP2R12Energy::get_debug() const { return debug_; };

double MP2R12Energy::energy()
{
  double value = emp2tot_aa_() + emp2tot_ab_() + er12tot_aa_() + er12tot_ab_();
  return value;
}

double MP2R12Energy::emp2tot_aa_() const
{
  RefSCVector emp2_aa = r12eval_->emp2_aa();
  int nij = emp2_aa.dim().n();
  double value = 0;
  for(int ij=0; ij<nij; ij++)
    value += emp2_aa.get_element(ij);

  return value;
}

double MP2R12Energy::emp2tot_ab_() const
{
  RefSCVector emp2_ab = r12eval_->emp2_ab();
  int nij = emp2_ab.dim().n();
  double value = 0;
  for(int ij=0; ij<nij; ij++)
    value += emp2_ab.get_element(ij);

  return value;
}

double MP2R12Energy::er12tot_aa_()
{
  compute();
  int nij = er12_aa_.dim().n();
  double value = 0;
  for(int ij=0; ij<nij; ij++)
    value += er12_aa_.get_element(ij);

  return value;
}

double MP2R12Energy::er12tot_ab_()
{
  compute();
  int nij = er12_ab_.dim().n();
  double value = 0;
  for(int ij=0; ij<nij; ij++)
    value += er12_ab_.get_element(ij);

  return value;
}

void MP2R12Energy::compute()
{
  if (evaluated_)
    return;
  
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MessageGrp> msg = r12info->msg();
  int me = msg->me();
  int ntasks = msg->n();

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //

  // Need eigenvalues
  int nocc = r12info->nocc();
  int nfzc = r12info->nfzc();
  int nocc_act = r12info->nocc_act();
  RefDiagSCMatrix evalmat = r12eval_->evals();
  double* evals = new double[nocc_act];
  for(int i=nfzc; i<nocc; i++)
    evals[i-nfzc] = evalmat(i);
  evalmat = 0;

  // Get the intermediates
  RefSCMatrix Vaa = r12eval_->V_aa();
  RefSCMatrix Xaa = r12eval_->X_aa();
  RefSCMatrix Baa = r12eval_->B_aa();
  RefSCMatrix Vab = r12eval_->V_ab();
  RefSCMatrix Xab = r12eval_->X_ab();
  RefSCMatrix Bab = r12eval_->B_ab();
  RefSCVector emp2_aa = r12eval_->emp2_aa();
  RefSCVector emp2_ab = r12eval_->emp2_ab();

  // Prepare total and R12 pairs
  Ref<SCMatrixKit> localkit = Vaa.kit();
  RefSCDimension dim_aa = r12eval_->dim_aa();
  RefSCDimension dim_ab = r12eval_->dim_ab();
  int naa = dim_aa.n();
  int nab = dim_ab.n();
  emp2r12_aa_ = localkit->vector(dim_aa);
  emp2r12_ab_ = localkit->vector(dim_ab);
  er12_aa_ = localkit->vector(dim_aa);
  er12_ab_ = localkit->vector(dim_ab);
  double* er12_aa_vec = new double[naa];
  double* er12_ab_vec = new double[nab];
  bzerofast(er12_aa_vec,naa);
  bzerofast(er12_ab_vec,nab);

  //
  // Alpha-alpha pairs
  //
  
  // Allocate the B matrix:
  // 1) in MP2-R12/A the B matrix is the same for all pairs
  // 2) int MP2-R12/A' the B matrix is pair-specific
  RefSCMatrix Baa_inv = Baa.clone();
  if (stdapprox_ == LinearR12::StdApprox_A) {
    Baa_inv->assign(Baa);
    Baa_inv->gen_invert_this();
    if (debug_ > 1)
      Baa_inv.print("Inverse alpha-alpha MP2-R12/A B matrix");
  }
  
  int ij=0;
  for(int i=0; i<nocc_act; i++)
    for(int j=0; j<i; j++, ij++) {

      if (ij%ntasks != me)
        continue;

      RefSCVector Vaa_ij = Vaa.get_column(ij);

      // In MP2-R12/A' matrices B are pair-specific:
      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      if (stdapprox_ == LinearR12::StdApprox_Ap) {
        Baa_inv.assign(Baa);
        int kl=0;
        for(int k=0; k<nocc_act; k++)
          for(int l=0; l<k; l++, kl++) {
            int ow=0;
            for(int o=0; o<nocc_act; o++)
              for(int w=0; w<o; w++, ow++) {
                double fx = 0.5 * (evals[k] + evals[l] + evals[o] + evals[w] - 2.0*evals[i] - 2.0*evals[j]) *
                Xaa.get_element(kl,ow);
                Baa_inv.accumulate_element(kl,ow,fx);
              }
          }
        if (debug_ > 1)
          Baa_inv.print("Alpha-alpha MP2-R12/A' B matrix");

        Baa_inv->gen_invert_this();

        if (debug_ > 1)
          Baa_inv.print("Inverse alpha-alpha MP2-R12/A' B matrix");
      }

      double eaa_ij = -2.0*Vaa_ij.dot(Baa_inv * Vaa_ij);
      er12_aa_vec[ij] = eaa_ij;
    }
  Baa_inv = 0;
  msg->sum(er12_aa_vec,naa,0,-1);
  er12_aa_->assign(er12_aa_vec);
  emp2r12_aa_->assign(emp2_aa);
  emp2r12_aa_->accumulate(er12_aa_);
  delete[] er12_aa_vec;

  //
  // Alpha-beta pairs
  //

  RefSCMatrix Bab_ij = Bab.clone();
  // In MP2-R12/A the B matrix is the same for all pairs
  if (stdapprox_ == LinearR12::StdApprox_A) {
    Bab_ij.assign(Bab);
    if (debug_ > 1)
      Bab_ij.print("Inverse alpha-beta MP2-R12/A B matrix");
    Bab_ij->gen_invert_this();
  }

  ij=0;
  for(int i=0; i<nocc_act; i++)
    for(int j=0; j<nocc_act; j++, ij++) {

      if (ij%ntasks != me)
        continue;

      RefSCVector Vab_ij = Vab.get_column(ij);

      // In MP2-R12/A' matrices B are pair-specific:
      // Form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      if (stdapprox_ == LinearR12::StdApprox_Ap) {
        Bab_ij.assign(Bab);
        int kl=0;
        for(int k=0; k<nocc_act; k++)
          for(int l=0; l<nocc_act; l++, kl++) {
            int ow=0;
            for(int o=0; o<nocc_act; o++)
              for(int w=0; w<nocc_act; w++, ow++) {
                double fx = 0.5 * (evals[k] + evals[l] + evals[o] + evals[w] - 2.0*evals[i] - 2.0*evals[j]) *
                Xab.get_element(kl,ow);
                Bab_ij.accumulate_element(kl,ow,fx);
              }
          }
        if (debug_ > 1)
	    Bab_ij.print("Alpha-beta MP2-R12/A' B matrix");

        Bab_ij->gen_invert_this();

        if (debug_ > 1)
          Bab_ij.print("Inverse alpha-beta MP2-R12/A' B matrix");
      }

      double eab_ij = -1.0*Vab_ij.dot(Bab_ij * Vab_ij);
      er12_ab_vec[ij] = eab_ij;
    }
  Bab_ij=0;
  msg->sum(er12_ab_vec,nab,0,-1);
  er12_ab_->assign(er12_ab_vec);
  emp2r12_ab_->assign(emp2_ab);
  emp2r12_ab_->accumulate(er12_ab_);
  delete[] er12_ab_vec;

  delete[] evals;

  evaluated_ = true;
  
  return;
}

void MP2R12Energy::print(std::ostream& so) const
{
} 

void MP2R12Energy::print_pair_energies(bool spinadapted, std::ostream& so)
{
  compute();
  
  char* SA_str;
  switch (stdapprox_) {
    case LinearR12::StdApprox_A:
      SA_str = strdup("A");
      break;

    case LinearR12::StdApprox_Ap:
      SA_str = strdup("A'");
      break;

    case LinearR12::StdApprox_B:
      SA_str = strdup("B");
      break;

    default:
      throw std::runtime_error("MP2R12Energy::print_pair_energies -- stdapprox_ is not valid");
  }

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  int nocc_act = r12info->nocc_act();
  double escf = r12info->ref()->energy();

  double emp2tot_aa = 0.0;
  double emp2tot_ab = 0.0;
  double er12tot_aa = 0.0;
  double er12tot_ab = 0.0;
  double emp2tot_0 = 0.0;
  double emp2tot_1 = 0.0;
  double er12tot_0 = 0.0;
  double er12tot_1 = 0.0;

  RefSCVector emp2_aa = r12eval_->emp2_aa();
  RefSCVector emp2_ab = r12eval_->emp2_ab();

  /*---------------------------------------
    Spin-adapt pair energies, if necessary
   ---------------------------------------*/
  if (!spinadapted) {

    so << endl << indent << "Alpha-alpha MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                 emp2_aa->get_element(ij),
                                 er12_aa_->get_element(ij),
                                 emp2r12_aa_->get_element(ij)) << endl;
      }
      
    so << endl << indent << "Alpha-beta MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<nocc_act;j++,ij++) {
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                 emp2_ab->get_element(ij),
                                 er12_ab_->get_element(ij),
                                 emp2r12_ab_->get_element(ij)) << endl;
      }

  }
  else {

    Ref<SCMatrixKit> localkit = er12_aa_.kit();
    RefSCVector emp2r12_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector emp2r12_1 = localkit->vector(r12eval_->dim_t());
    RefSCVector emp2_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector emp2_1 = localkit->vector(r12eval_->dim_t());
    RefSCVector er12_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector er12_1 = localkit->vector(r12eval_->dim_t());

    // Triplet pairs are easy
    emp2r12_1->assign(emp2r12_aa_);
    emp2r12_1->scale(1.5);
    emp2_1->assign(emp2_aa);
    emp2_1->scale(1.5);
    er12_1->assign(er12_aa_);
    er12_1->scale(1.5);

    // Singlet pairs are a bit trickier
    int ij_s = 0;
    for(int i=0; i<nocc_act; i++)
      for(int j=0; j<=i; j++, ij_s++) {
        int ij_ab = i*nocc_act + j;
        int ij_aa = i*(i-1)/2 + j;
        double eab, eaa, e_s;

        eab = emp2r12_ab_->get_element(ij_ab);
        if (i != j)
          eaa = emp2r12_aa_->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        emp2r12_0->set_element(ij_s,e_s);

        eab = emp2_ab->get_element(ij_ab);
        if (i != j)
          eaa = emp2_aa->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        emp2_0->set_element(ij_s,e_s);

        eab = er12_ab_->get_element(ij_ab);
        if (i != j)
          eaa = er12_aa_->get_element(ij_aa);
        else
          eaa = 0.0;
        e_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
        er12_0->set_element(ij_s,e_s);
      }

    // compute total singlet and triplet energies
    RefSCVector unit_0 = localkit->vector(r12eval_->dim_s());
    RefSCVector unit_1 = localkit->vector(r12eval_->dim_t());
    unit_0->assign(1.0);
    unit_1->assign(1.0);
    emp2tot_0 = emp2_0.dot(unit_0);
    emp2tot_1 = emp2_1.dot(unit_1);
    er12tot_0 = er12_0.dot(unit_0);
    er12tot_1 = er12_1.dot(unit_1);

    so << endl << indent << "Singlet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<=i;j++,ij++) {
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                 emp2_0->get_element(ij),
                                 er12_0->get_element(ij),
                                 emp2r12_0->get_element(ij)) << endl;
      }
      
    so << endl << indent << "Triplet MBPT2-R12/" << SA_str << " pair energies:" << endl;
    so << indent << scprintf("    i       j        mp2(ij)        r12(ij)      mp2-r12(ij)") << endl;
    so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
    for(int i=0,ij=0;i<nocc_act;i++)
      for(int j=0;j<i;j++,ij++) {
        so << indent << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf",i+1,j+1,
                                 emp2_1->get_element(ij),
                                 er12_1->get_element(ij),
                                 emp2r12_1->get_element(ij)) << endl;
      }

  }


  double mp2_corr_energy_ = emp2tot_aa_() + emp2tot_ab_();
  double r12_corr_energy_ = er12tot_aa_() + er12tot_ab_();
  
  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////

  if (spinadapted) {
    so <<endl<<indent
      <<scprintf("Singlet MP2 correlation energy [au]:          %17.12lf\n", emp2tot_0);
    so <<indent
      <<scprintf("Triplet MP2 correlation energy [au]:          %17.12lf\n", emp2tot_1);
    so <<indent
      <<scprintf("Singlet (MP2)-R12/%2s correlation energy [au]: %17.12lf\n", SA_str, er12tot_0);
    so <<indent
      <<scprintf("Triplet (MP2)-R12/%2s correlation energy [au]: %17.12lf\n", SA_str, er12tot_1);
    so <<indent
      <<scprintf("Singlet MP2-R12/%2s correlation energy [au]:   %17.12lf\n", SA_str,
                 emp2tot_0 + er12tot_0);
    so <<indent
      <<scprintf("Triplet MP2-R12/%2s correlation energy [au]:   %17.12lf\n", SA_str,
                 emp2tot_1 + er12tot_1);
  }
  
  double etotal = escf + mp2_corr_energy_ + r12_corr_energy_;
  so <<endl<<indent
    <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
  so <<indent
    <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", mp2_corr_energy_);
  so <<indent
    <<scprintf("(MBPT2)-R12/%2s correlation energy [au]:    %17.12lf\n", SA_str, r12_corr_energy_);
  so <<indent
    <<scprintf("MBPT2-R12/%2s correlation energy [au]:      %17.12lf\n", SA_str,
               mp2_corr_energy_ + r12_corr_energy_);
  so <<indent
    <<scprintf("MBPT2-R12/%2s energy [au]:                  %17.12lf\n", SA_str, etotal) << endl;

  so.flush();

  return;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
