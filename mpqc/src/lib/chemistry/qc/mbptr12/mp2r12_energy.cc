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

#include <ostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/svd.h>

using namespace std;
using namespace sc;
using namespace sc::exp;

inline int max(int a,int b) { return (a > b) ? a : b;}

#define USE_INVERT 0

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

  init_();
}

void
MP2R12Energy::init_()
{

  RefSCDimension dim_oo_aa = r12eval_->dim_oo_aa();
  RefSCDimension dim_oo_ab = r12eval_->dim_oo_ab();
  Ref<SCMatrixKit> kit = r12eval_->r12info()->matrixkit();
  er12_aa_ = kit->vector(dim_oo_aa);
  er12_ab_ = kit->vector(dim_oo_ab);
  emp2r12_aa_ = kit->vector(dim_oo_aa);
  emp2r12_ab_ = kit->vector(dim_oo_ab);

  RefSCDimension dim_vv_aa = r12eval_->dim_vv_aa();
  RefSCDimension dim_vv_ab = r12eval_->dim_vv_ab();
  Caa_ = kit->matrix(dim_oo_aa, dim_oo_aa);
  Cab_ = kit->matrix(dim_oo_ab, dim_oo_ab);

} 

MP2R12Energy::MP2R12Energy(StateIn& si) : SavableState(si)
{
  r12eval_ << SavableState::restore_state(si);

  init_();

  er12_aa_.restore(si);
  er12_ab_.restore(si);
  emp2r12_aa_.restore(si);
  emp2r12_ab_.restore(si);

  Caa_.restore(si);
  Cab_.restore(si);
  
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
  
  Caa_.save(so);
  Cab_.save(so);
  
  so.put((int)stdapprox_);
  so.put(debug_);
  so.put((int)evaluated_);
}

void MP2R12Energy::obsolete()
{
  evaluated_ = false;
}

Ref<R12IntEval> MP2R12Energy::r12eval() const { return r12eval_; };
const bool MP2R12Energy::ebc() const { return ebc_; };
const bool MP2R12Energy::gbc() const { return r12eval_->gbc(); };
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
  
  bool ebc = r12eval_->ebc();

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //

  // Need eigenvalues
  int nocc = r12info->ndocc();
  int nfzc = r12info->nfzc();
  int nocc_act = r12info->ndocc_act();
  int nvir_act = r12info->nvir_act();
  RefDiagSCMatrix evalmat = r12eval_->evals();
  vector<double> evals_act_occ(nocc_act);
  vector<double> evals_act_vir(nvir_act);
  for(int i=nfzc; i<nocc; i++)
    evals_act_occ[i-nfzc] = evalmat(i);
  for(int i=0; i<nvir_act; i++)
    evals_act_vir[i] = evalmat(i+nocc);
  evalmat = 0;

  // Get the intermediates
  RefSCMatrix Vaa = r12eval_->V_aa();
  RefSCMatrix Xaa = r12eval_->X_aa();
  RefSymmSCMatrix Baa = r12eval_->B_aa();
  RefSCMatrix Aaa = r12eval_->A_aa();
  RefSCMatrix Vab = r12eval_->V_ab();
  RefSCMatrix Xab = r12eval_->X_ab();
  RefSymmSCMatrix Bab = r12eval_->B_ab();
  RefSCMatrix Aab = r12eval_->A_ab();
  RefSCVector emp2_aa = r12eval_->emp2_aa();
  RefSCVector emp2_ab = r12eval_->emp2_ab();

  // Prepare total and R12 pairs
  Ref<SCMatrixKit> localkit = Vaa.kit();
  RefSCDimension dim_oo_aa = r12eval_->dim_oo_aa();
  RefSCDimension dim_oo_ab = r12eval_->dim_oo_ab();
  int naa = dim_oo_aa.n();
  int nab = dim_oo_ab.n();
  emp2r12_aa_ = localkit->vector(dim_oo_aa);
  emp2r12_ab_ = localkit->vector(dim_oo_ab);
  er12_aa_ = localkit->vector(dim_oo_aa);
  er12_ab_ = localkit->vector(dim_oo_ab);
  double* er12_aa_vec = new double[naa];
  double* er12_ab_vec = new double[nab];
  bzerofast(er12_aa_vec,naa);
  bzerofast(er12_ab_vec,nab);

  //
  // Alpha-alpha pairs
  //
  if (naa > 0) {
    if (debug_ > 1) {
      Vaa.print("Alpha-alpha V matrix");
      Baa.print("Alpha-alpha MP2-R12/A B matrix");
      if (ebc == false)
        Aaa.print("Alpha-alpha A matrix");
    }

    // Allocate the B matrix:
    // 1) in MP2-R12/A the B matrix is the same for all pairs
    // 2) int MP2-R12/A' the B matrix is pair-specific
    RefSymmSCMatrix Baa_ij = Baa.clone();
    if (stdapprox_ == LinearR12::StdApprox_A) {
#if USE_INVERT
      Baa_ij->assign(Baa);
      Baa_ij->gen_invert_this();
      if (debug_ > 1)
        Baa_ij.print("Inverse alpha-alpha MP2-R12/A B matrix");
#else
      // solve B * C = V
      RefSCMatrix Caa_kl_by_ij = Caa_.clone();
      sc::exp::lapack_linsolv_symmnondef(Baa, Caa_kl_by_ij, Vaa);
      Caa_kl_by_ij = Caa_kl_by_ij.t();
      Caa_.assign(Caa_kl_by_ij);  Caa_kl_by_ij = 0;
      Caa_.scale(-1.0);
#endif
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
          Baa_ij.assign(Baa);
          int kl=0;
          for(int k=0; k<nocc_act; k++)
            for(int l=0; l<k; l++, kl++) {
              int ow=0;
              for(int o=0; o<nocc_act; o++)
                for(int w=0; w<o; w++, ow++) {
                  
                  if (ow > kl)
                    continue;
                  
                  double fx = 0.5 * (evals_act_occ[k] + evals_act_occ[l] + evals_act_occ[o] + evals_act_occ[w]
                                     - 2.0*evals_act_occ[i] - 2.0*evals_act_occ[j]) *
                    Xaa.get_element(kl,ow);
                  
                  Baa_ij.accumulate_element(kl,ow,fx);
                  
                  // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ei-ej)
                  if (ebc == false) {
                    double fy = 0.0;
                    int cd=0;
                    for(int c=0; c<nvir_act; c++)
                      for(int d=0; d<c; d++, cd++) {
                        
                        fy -= Aaa.get_element(kl,cd)*Aaa.get_element(ow,cd)/(evals_act_vir[c] + evals_act_vir[d]
                                                                             - evals_act_occ[i] - evals_act_occ[j]);
                      }
                    
                    Baa_ij.accumulate_element(kl,ow,fy);
                  }
                  
                }
            }
          if (debug_ > 1)
            Baa_ij.print("Alpha-alpha MP2-R12/A' B matrix");
          
#if USE_INVERT
          Baa_ij->gen_invert_this();
          
          if (debug_ > 1)
            Baa_ij.print("Inverse alpha-alpha MP2-R12/A' B matrix");
#endif
          
        }
        
#if USE_INVERT
        // The r12 amplitudes B^-1 * V
        RefSCVector Cij = -1.0*(Baa_ij * Vaa_ij);
        const int nkl = Cij.dim().n();
        for(int kl=0; kl<nkl; kl++)
          Caa_.set_element(ij,kl,Cij.get_element(kl));
#else
        RefSCVector Cij = Vaa_ij.clone();
        if (stdapprox_ == LinearR12::StdApprox_A) {
          double* v = new double[Cij.n()];
          Caa_.get_row(ij).convert(v);
          Cij.assign(v);
          delete[] v;
        }
        else {
          // solve B * C = V
          Cij = Vaa_ij.clone();
          sc::exp::lapack_linsolv_symmnondef(Baa_ij, Cij, Vaa_ij);
          Cij.scale(-1.0);
          const int nkl = Cij.dim().n();
          for(int kl=0; kl<nkl; kl++)
            Caa_.set_element(ij,kl,Cij.get_element(kl));
        }
#endif
        double eaa_ij = 2.0*Vaa_ij.dot(Cij);
        er12_aa_vec[ij] = eaa_ij;
      }
    Baa_ij = 0;
    msg->sum(er12_aa_vec,naa,0,-1);
    er12_aa_->assign(er12_aa_vec);
    emp2r12_aa_->assign(emp2_aa);
    emp2r12_aa_->accumulate(er12_aa_);
    delete[] er12_aa_vec;
  }

  //
  // Alpha-beta pairs
  //
  if (nab > 0) {
    if (debug_ > 1) {
      Vab.print("Alpha-beta V matrix");
      Bab.print("Alpha-beta MP2-R12/A B matrix");
      if (ebc == false)
        Aab.print("Alpha-beta A matrix");
    }
    
    RefSymmSCMatrix Bab_ij = Bab.clone();
    // In MP2-R12/A the B matrix is the same for all pairs
    if (stdapprox_ == LinearR12::StdApprox_A) {
#if USE_INVERT
      Bab_ij.assign(Bab);
      if (debug_ > 1)
        Bab_ij.print("Inverse alpha-beta MP2-R12/A B matrix");
      Bab_ij->gen_invert_this();
#else
      // solve B * C = V
      RefSCMatrix Cab_kl_by_ij = Cab_.clone();
      sc::exp::lapack_linsolv_symmnondef(Bab, Cab_kl_by_ij, Vab);
      Cab_kl_by_ij = Cab_kl_by_ij.t();
      Cab_.assign(Cab_kl_by_ij);  Cab_kl_by_ij = 0;
      Cab_.scale(-1.0);
#endif
    }
    
    int ij=0;
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
                  
                  if (ow > kl)
                    continue;
                  
                  double fx = 0.5 * (evals_act_occ[k] + evals_act_occ[l] + evals_act_occ[o] + evals_act_occ[w]
                                     - 2.0*evals_act_occ[i] - 2.0*evals_act_occ[j]) *
                    Xab.get_element(kl,ow);
                  Bab_ij.accumulate_element(kl,ow,fx);
                  
                  // If EBC is not assumed add Akl,cd*Acd,ow/(ec+ed-ei-ej)
                  if (ebc == false) {
                    double fy = 0.0;
                    int cd=0;
                    for(int c=0; c<nvir_act; c++)
                      for(int d=0; d<nvir_act; d++, cd++) {
                        
                        fy -= Aab.get_element(kl,cd)*Aab.get_element(ow,cd)/(evals_act_vir[c] + evals_act_vir[d]
                                                                             - evals_act_occ[i] - evals_act_occ[j]);
                      }
                    
                    Bab_ij.accumulate_element(kl,ow,fy);
                  }
                  
                }
            }
          if (debug_ > 1)
	    Bab_ij.print("Alpha-beta MP2-R12/A' B matrix");
          
#if USE_INVERT
          Bab_ij->gen_invert_this();
          
          if (debug_ > 1)
            Bab_ij.print("Inverse alpha-beta MP2-R12/A' B matrix");
#endif
          
        }

#if USE_INVERT
        // the r12 amplitudes B^-1 * V
        RefSCVector Cij = -1.0*(Bab_ij * Vab_ij);
        const int nkl = Cij.dim().n();
        for(int kl=0; kl<nkl; kl++)
          Cab_.set_element(ij,kl,Cij.get_element(kl));
#else
        RefSCVector Cij = Vab_ij.clone();
        if (stdapprox_ == LinearR12::StdApprox_A) {
          double* v = new double[Cij.n()];
          Cab_.get_row(ij).convert(v);
          Cij.assign(v);
          delete[] v;
        }
        else {
          // solve B * C = V
          Cij = Vab_ij.clone();
          sc::exp::lapack_linsolv_symmnondef(Bab_ij, Cij, Vab_ij);
          Cij.scale(-1.0);
          const int nkl = Cij.dim().n();
          for(int kl=0; kl<nkl; kl++)
            Cab_.set_element(ij,kl,Cij.get_element(kl));
        }
#endif
        double eab_ij = 1.0*Vab_ij.dot(Cij);
        er12_ab_vec[ij] = eab_ij;
      }
    Bab_ij=0;
    msg->sum(er12_ab_vec,nab,0,-1);
    er12_ab_->assign(er12_ab_vec);
    emp2r12_ab_->assign(emp2_ab);
    emp2r12_ab_->accumulate(er12_ab_);
    delete[] er12_ab_vec;
  }

  evaluated_ = true;
  
  return;
}

static void print_psi_values(std::ostream& fout, const SCVector3& r1, const SCVector3& r2, double phi_0, double phi_1_mp2, double phi_1_r12)
{
  fout << scprintf("%9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf %12.8lf %12.8lf %12.8lf",
                   r1.x(),r1.y(),r1.z(),r2.x(),r2.y(),r2.z(),phi_0,phi_1_mp2,phi_1_r12) << endl;
}

double
MP2R12Energy::compute_pair_function_aa(int ij, const SCVector3& r1, const SCVector3& r2)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2aa = Amps->T2_aa();
  RefSCMatrix Rvv_aa = Amps->Rvv_aa();
  RefSCMatrix Roo_aa = Amps->Roo_aa();
  RefSCMatrix Rvo_aa = Amps->Rvo_aa();
  RefSCMatrix Rxo_aa = Amps->Rxo_aa();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Caa = localkit->matrix(Caa_.rowdim(),Caa_.coldim());
  double* caa = new double[Caa_.rowdim().n()*Caa_.coldim().n()];
  Caa_.convert(caa);
  Caa.assign(caa);
  delete[] caa;
  RefSCMatrix Cvv = Caa * Rvv_aa;
  RefSCMatrix Coo = Caa * Roo_aa;
  RefSCMatrix Cov = Caa * Rvo_aa;
  RefSCMatrix Cox = Caa * Rxo_aa;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->act_vir_space();
  Ref<MOIndexSpace> act_occ_space = r12info->act_occ_space();
  Ref<MOIndexSpace> occ_space = r12info->occ_space();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  RefSCVector phi_aa = compute_2body_values_(true,act_occ_space,act_occ_space,r1,r2);
  RefSCVector phi_vv = compute_2body_values_(true,act_vir_space,act_vir_space,r1,r2);
  RefSCVector phi_oo = compute_2body_values_(true,occ_space,occ_space,r1,r2);
  RefSCVector phi_ov = compute_2body_values_(true,occ_space,act_vir_space,r1,r2);
  RefSCVector phi_ox = compute_2body_values_(true,occ_space,ribs_space,r1,r2);

  double phi_t2 = T2aa.get_row(ij).dot(phi_vv);

  SCVector3 r12 = r1 - r2;
  const double dist12 = r12.norm();
  double phi_r12;
  phi_r12 = 0.5 * Caa.get_row(ij).dot(phi_aa) * dist12;
  phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
  phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
  phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
  phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

  print_psi_values(ExEnv::out0(),r1,r2,phi_aa.get_element(ij),phi_t2,phi_r12);

}

void
MP2R12Energy::compute_pair_function_aa(int ij, const Ref<TwoBodyGrid>& tbgrid)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2aa = Amps->T2_aa();
  RefSCMatrix Rvv_aa = Amps->Rvv_aa();
  RefSCMatrix Roo_aa = Amps->Roo_aa();
  RefSCMatrix Rvo_aa = Amps->Rvo_aa();
  RefSCMatrix Rxo_aa = Amps->Rxo_aa();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Caa = localkit->matrix(Caa_.rowdim(),Caa_.coldim());
  double* caa = new double[Caa_.rowdim().n()*Caa_.coldim().n()];
  Caa_.convert(caa);
  Caa.assign(caa);
  delete[] caa;
  RefSCMatrix Cvv = Caa * Rvv_aa;
  RefSCMatrix Coo = Caa * Roo_aa;
  RefSCMatrix Cov = Caa * Rvo_aa;
  RefSCMatrix Cox = Caa * Rxo_aa;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->act_vir_space();
  Ref<MOIndexSpace> act_occ_space = r12info->act_occ_space();
  Ref<MOIndexSpace> occ_space = r12info->occ_space();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  const int nelem = tbgrid->nelem();
  std::stringstream output_file_name;
  output_file_name << tbgrid->name() << ".ab.pair"
                   << ij << ".txt";
  ofstream ofile(output_file_name.str().c_str());

  for(int i=0; i<nelem; i++) {
    RefSCVector phi_aa = compute_2body_values_(true,act_occ_space,act_occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_vv = compute_2body_values_(true,act_vir_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_oo = compute_2body_values_(true,occ_space,occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ov = compute_2body_values_(true,occ_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ox = compute_2body_values_(true,occ_space,ribs_space,tbgrid->xyz1(i),tbgrid->xyz2(i));

    double phi_t2 = T2aa.get_row(ij).dot(phi_vv);

    SCVector3 r12 = tbgrid->xyz1(i) - tbgrid->xyz2(i);
    const double dist12 = r12.norm();
    double phi_r12;
    phi_r12 = 0.5 * Caa.get_row(ij).dot(phi_aa) * dist12;
    phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
    phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
    phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
    phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

    print_psi_values(ofile,tbgrid->xyz1(i),tbgrid->xyz2(i),phi_aa.get_element(ij),phi_t2,phi_r12);
  }
}

double
MP2R12Energy::compute_pair_function_ab(int ij, const SCVector3& r1, const SCVector3& r2)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2ab = Amps->T2_ab();
  RefSCMatrix Rvv_ab = Amps->Rvv_ab();
  RefSCMatrix Roo_ab = Amps->Roo_ab();
  RefSCMatrix Rvo_ab = Amps->Rvo_ab();
  RefSCMatrix Rxo_ab = Amps->Rxo_ab();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Cab = localkit->matrix(Cab_.rowdim(),Cab_.coldim());
  double* cab = new double[Cab_.rowdim().n()*Cab_.coldim().n()];
  Cab_.convert(cab);
  Cab.assign(cab);
  delete[] cab;
  RefSCMatrix Cvv = Cab * Rvv_ab;
  RefSCMatrix Coo = Cab * Roo_ab;
  RefSCMatrix Cov = Cab * Rvo_ab;
  RefSCMatrix Cox = Cab * Rxo_ab;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->act_vir_space();
  Ref<MOIndexSpace> act_occ_space = r12info->act_occ_space();
  Ref<MOIndexSpace> occ_space = r12info->occ_space();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  RefSCVector phi_aa = compute_2body_values_(false,act_occ_space,act_occ_space,r1,r2);
  RefSCVector phi_vv = compute_2body_values_(false,act_vir_space,act_vir_space,r1,r2);
  RefSCVector phi_oo = compute_2body_values_(false,occ_space,occ_space,r1,r2);
  RefSCVector phi_ov = compute_2body_values_(false,occ_space,act_vir_space,r1,r2);
  RefSCVector phi_ox = compute_2body_values_(false,occ_space,ribs_space,r1,r2);

  double phi_t2 = T2ab.get_row(ij).dot(phi_vv);

  SCVector3 r12 = r1 - r2;
  const double dist12 = r12.norm();
  double phi_r12;
  phi_r12 = 0.5*Cab.get_row(ij).dot(phi_aa) * dist12;
  phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
  phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
  phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
  phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

  print_psi_values(ExEnv::out0(),r1,r2,phi_aa.get_element(ij),phi_t2,phi_r12);
}

void
MP2R12Energy::compute_pair_function_ab(int ij, const Ref<TwoBodyGrid>& tbgrid)
{
  Ref<R12Amplitudes> Amps = r12eval_->amps();
  RefSCMatrix T2ab = Amps->T2_ab();
  RefSCMatrix Rvv_ab = Amps->Rvv_ab();
  RefSCMatrix Roo_ab = Amps->Roo_ab();
  RefSCMatrix Rvo_ab = Amps->Rvo_ab();
  RefSCMatrix Rxo_ab = Amps->Rxo_ab();

  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  RefSCMatrix Cab = localkit->matrix(Cab_.rowdim(),Cab_.coldim());
  double* cab = new double[Cab_.rowdim().n()*Cab_.coldim().n()];
  Cab_.convert(cab);
  Cab.assign(cab);
  delete[] cab;
  RefSCMatrix Cvv = Cab * Rvv_ab;
  RefSCMatrix Coo = Cab * Roo_ab;
  RefSCMatrix Cov = Cab * Rvo_ab;
  RefSCMatrix Cox = Cab * Rxo_ab;

  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MOIndexSpace> act_vir_space = r12info->act_vir_space();
  Ref<MOIndexSpace> act_occ_space = r12info->act_occ_space();
  Ref<MOIndexSpace> occ_space = r12info->occ_space();
  Ref<MOIndexSpace> ribs_space = r12info->ribs_space();

  const int nelem = tbgrid->nelem();
  std::stringstream output_file_name;
  output_file_name << tbgrid->name() << ".ab.pair"
                   << ij << ".txt";
  ofstream ofile(output_file_name.str().c_str());

  for(int i=0; i<nelem; i++) {
    RefSCVector phi_aa = compute_2body_values_(false,act_occ_space,act_occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_vv = compute_2body_values_(false,act_vir_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_oo = compute_2body_values_(false,occ_space,occ_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ov = compute_2body_values_(false,occ_space,act_vir_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    RefSCVector phi_ox = compute_2body_values_(false,occ_space,ribs_space,tbgrid->xyz1(i),tbgrid->xyz2(i));
    
    double phi_t2 = T2ab.get_row(ij).dot(phi_vv);
    
    SCVector3 r12 = tbgrid->xyz1(i) - tbgrid->xyz2(i);
    const double dist12 = r12.norm();
    double phi_r12;
    phi_r12 = 0.5*Cab.get_row(ij).dot(phi_aa) * dist12;
    phi_r12 -= 0.5 * Cvv.get_row(ij).dot(phi_vv);
    phi_r12 -= 0.5 * Coo.get_row(ij).dot(phi_oo);
    phi_r12 -= 1.0 * Cov.get_row(ij).dot(phi_ov);
    phi_r12 -= 1.0 * Cox.get_row(ij).dot(phi_ox);

    print_psi_values(ofile,tbgrid->xyz1(i),tbgrid->xyz2(i),phi_aa.get_element(ij),phi_t2,phi_r12);
  }
}

RefSCVector
MP2R12Energy::compute_2body_values_(bool equiv, const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                    const SCVector3& r1, const SCVector3& r2) const
{
  const Ref<Integral> ints = r12eval_->r12info()->integral();
  const Ref<GaussianBasisSet> bs1 = space1->basis();
  const Ref<GaussianBasisSet> bs2 = space2->basis();
  ints->set_basis(bs1,bs2);
  GaussianBasisSet::ValueData* vdata1 = new GaussianBasisSet::ValueData(bs1,ints);
  GaussianBasisSet::ValueData* vdata2 = new GaussianBasisSet::ValueData(bs2,ints);

  const bool space1_eq_space2 = (space1 == space2);
  const int nbasis1 = bs1->nbasis();
  const int nbasis2 = bs2->nbasis();
  const int rank1 = space1->rank();
  const int rank2 = space2->rank();

  const int npair = (space1_eq_space2 && equiv) ? rank1*(rank1-1)/2 : rank1*rank2;
  RefSCDimension pairdim = new SCDimension(npair);

  double* values11 = new double[nbasis1];
  double* values12 = new double[nbasis1];
  double* values21 = new double[nbasis2];
  double* values22 = new double[nbasis2];

  bs1->values(r1,vdata1,values11);
  bs1->values(r2,vdata1,values12);
  bs2->values(r1,vdata2,values21);
  bs2->values(r2,vdata2,values22);

  RefSCMatrix ao2mo_1 = space1->coefs().t();
  RefSCMatrix ao2mo_2 = space2->coefs().t();

  Ref<SCMatrixKit> kit = ao2mo_1.kit();
  RefSCVector vals11 = kit->vector(ao2mo_1.coldim());
  RefSCVector vals12 = kit->vector(ao2mo_1.coldim());
  RefSCVector vals21 = kit->vector(ao2mo_2.coldim());
  RefSCVector vals22 = kit->vector(ao2mo_2.coldim());
  vals11.assign(values11);
  vals12.assign(values12);
  vals21.assign(values21);
  vals22.assign(values22);
  delete[] values11;
  delete[] values12;
  delete[] values21;
  delete[] values22;

  RefSCVector movals11 = ao2mo_1 * vals11;
  RefSCVector movals12 = ao2mo_1 * vals12;
  RefSCVector movals21 = ao2mo_2 * vals21;
  RefSCVector movals22 = ao2mo_2 * vals22;

  kit = new LocalSCMatrixKit;
  RefSCVector vals = kit->vector(pairdim);
  
  MOPairIterFactory PIFactory;
  Ref<SpatialMOPairIter> ij_iter = PIFactory.mopairiter(space1,space2);
  for(ij_iter->start();int(*ij_iter.pointer());ij_iter->next()) {
    const int i = ij_iter->i();
    const int j = ij_iter->j();
    const int ij_aa = ij_iter->ij_aa();
    const int ij_ab = ij_iter->ij_ab();
    const int ij_ba = ij_iter->ij_ba();

    if (equiv) {
      if (ij_aa != -1) {
        const double value = movals11.get_element(i) * movals22.get_element(j) -
          movals12.get_element(i) * movals21.get_element(j);
        vals.set_element(ij_aa,value);
      }
    }
    else {
      const double value = movals11.get_element(i) * movals22.get_element(j);
      vals.set_element(ij_ab,value);
      if (space1_eq_space2 && ij_ab != ij_ba) {
        const double value = movals11.get_element(j) * movals22.get_element(i);
        vals.set_element(ij_ba,value);
      }
    }
        
  }

  vdata1->~ValueData();
  vdata2->~ValueData();

  return vals;
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
  int nocc_act = r12info->ndocc_act();
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
    RefSCVector emp2r12_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector emp2r12_1 = localkit->vector(r12eval_->dim_oo_t());
    RefSCVector emp2_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector emp2_1 = localkit->vector(r12eval_->dim_oo_t());
    RefSCVector er12_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector er12_1 = localkit->vector(r12eval_->dim_oo_t());

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
    RefSCVector unit_0 = localkit->vector(r12eval_->dim_oo_s());
    RefSCVector unit_1 = localkit->vector(r12eval_->dim_oo_t());
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

  free(SA_str);
  
  return;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
