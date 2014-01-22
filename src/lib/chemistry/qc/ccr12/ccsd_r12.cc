//
// ccsd_r12.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: TS
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
#include <util/misc/regtime.h>
#include <util/state/stateio.h>
#include <math/optimize/diis.h>
#include <chemistry/qc/ccr12/ccsd_r12.h>
#include <chemistry/qc/ccr12/ccsd_r12_e.h>
#include <chemistry/qc/ccr12/ccsd_r12_t1.h>
#include <chemistry/qc/ccr12/ccsd_r12_t2.h>

using namespace std;
using namespace sc;


static ClassDesc CCSD_R12_cd(
  typeid(CCSD_R12), "CCSD_R12", 1, "public CCR12",
  0,create<CCSD_R12>, create<CCSD_R12>);


CCSD_R12::CCSD_R12(StateIn& s): CCR12(s){
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructer not yet implemented", __FILE__, __LINE__);
}


CCSD_R12::CCSD_R12(const Ref<KeyVal>& keyval): CCR12(keyval){
  string theory("CCSD-R12");
  theory_ = theory;

  perturbative_ = keyval->stringvalue("perturbative", KeyValValuestring(""));
  ostringstream oss; oss.setf(ios_base::uppercase);
  oss << perturbative_;
  perturbative_ = oss.str();

  print_theory();
}


CCSD_R12::~CCSD_R12(){
}


void CCSD_R12::compute(){

  CCR12::compute();

  Ref<Tensor> e0 = new Tensor("e",mem_);
  ccr12_info_->offset_e(e0);

  Ref<Tensor> r1 = new Tensor("r1",mem_);
  ccr12_info_->offset_t1(r1, false);
  Ref<Tensor> r2 = new Tensor("r2",mem_);
  ccr12_info_->offset_t2(r2, false);

  CCSD_R12_E*  ccsd_r12_e  = new CCSD_R12_E( info());
  CCSD_R12_T1* ccsd_r12_t1 = new CCSD_R12_T1(info());
  CCSD_R12_T2* ccsd_r12_t2 = new CCSD_R12_T2(info());

  string theory_ = "CCSD-R12";
  print_iteration_header(theory_);

  // only fixed amplitude ansatz is supported at the moment
  assert(ccr12_info_->r12world()->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt);
  assert(ccr12_info_->r12world()->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_scaledfixed);

  timer_->enter("CCR12 iterations");
  double iter_start = 0.0;
  double iter_end   = timer_->get_wall_time();
  double energy;


  Ref<DIIS> t1diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);
  Ref<DIIS> t2diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);

  // qy = F12 * gt2. Tilde t in our paper
  ccr12_info_->update_qy();

  for (int iter = 0; iter < maxiter_; ++iter){
    iter_start = iter_end;

    e0->zero();
    r1->zero();
    r2->zero();

    ccsd_r12_e->compute_amp(e0);
    ccsd_r12_t1->compute_amp(r1);
    ccsd_r12_t2->compute_amp(r2);

    energy = ccr12_info_->get_e(e0);

    // compute new amplitudes from the residuals
    Ref<Tensor> t1_old = info()->t1()->copy();
    Ref<Tensor> t2_old = info()->t2()->copy();
    ccr12_info_->jacobi_t1(r1);
    ccr12_info_->jacobi_t2(r2);
    // compute errors
    Ref<Tensor> t1_err = t1_old;
    Ref<Tensor> t2_err = t2_old;
    t1_err->daxpy(info()->t1(), -1.0);
    t2_err->daxpy(info()->t2(), -1.0);

    const double r1norm = RMS(*t1_err);
    const double r2norm = RMS(*t2_err);
    const double rnorm  = std::sqrt(r1norm*r1norm+
                                    r2norm*r2norm);

    iter_end = timer_->get_wall_time();
    print_iteration(iter, energy, rnorm, iter_start, iter_end);

    // done? break free
    if (rnorm<ccthresh_) break;

    // extrapolate
    t1diis->extrapolate(info()->edata(info()->t1()), info()->eerr(t1_err));
    t2diis->extrapolate(info()->edata(info()->t2()), info()->eerr(t2_err));



  }
  timer_->exit("CCR12 iterations");

  print_iteration_footer();

  // if not fully optimized, we need to add the geminal-lambda contribution.
  if (ccr12_info_->r12world()->r12tech()->ansatz()->amplitudes() != R12Technology::GeminalAmplitudeAnsatz_fullopt) {
    assert(false);
    timer_->enter("Lambda contribution");
    iter_start = timer_->get_wall_time();

    e0->zero();
    ccr12_info_->guess_glambda2(); // glambda2 = gt2^dagger
// to be implemented >>
// qx = F^AA_ii * gt2
//
//  ccr12_info_->update_qx();
// <<<<<<<<<<<<<<<<<<<<

    // TODO implement ccsd_r12_c and use it
    // TODO not yet implemented
    const double lambda_corr = ccr12_info_->get_e(e0);
    print_correction(lambda_corr, energy, "Lambda contribution");
    energy += lambda_corr;

    print_timing(timer_->get_wall_time()-iter_start,"Lambda contribution");
    timer_->exit("Lambda contribution");
  }

  delete ccsd_r12_t2;
  delete ccsd_r12_t1;
  delete ccsd_r12_e;


// this was copied from CCSD::compute(); needs some modification;
// commented out for a while to use in the future development
#if 0
  if (perturbative_=="(T)") {
    timer_->enter("(T) correction");
    iter_start=timer_->get_wall_time();

    // ccsd(t) calls specific driver, because it is a bit
    // different from others (i.e. (2)t etc)
    Ref<CCSD_PT_LEFT>     eval_left=new CCSD_PT_LEFT(info());
    Ref<CCSD_PT_RIGHT>   eval_right=new CCSD_PT_RIGHT(info());
    Ref<CCSD_PT>            ccsd_pt=new CCSD_PT(info());
    const double ccsd_pt_correction=ccsd_pt->compute(eval_left,eval_right);
    print_correction(ccsd_pt_correction,energy,"CCSD(T)");

    print_timing(timer_->get_wall_time()-iter_start,"(T) correction");
    timer_->exit("(T) correction");
  }


/////////////////////////////////////////////////////////////////////////////////////////////////


  bool do_lambda=false; // will be judeged from input keywords
  // more will come; e.g. dipole, etc
  if (perturbative_=="(2)T" or perturbative_=="(2)TQ") do_lambda=true;

  if (do_lambda) {
    Ref<DIIS> l1diis=new DIIS(diis_start_,ndiis_,0.005,3,1,0.0);
    Ref<DIIS> l2diis=new DIIS(diis_start_,ndiis_,0.005,3,1,0.0);

    LAMBDA_CCSD_T1* lambda_ccsd_t1=new LAMBDA_CCSD_T1(info());
    LAMBDA_CCSD_T2* lambda_ccsd_t2=new LAMBDA_CCSD_T2(info());

    Ref<Tensor> lr1=new Tensor("lr1",mem_);
    Ref<Tensor> lr2=new Tensor("lr2",mem_);
    ccr12_info_->offset_l1(lr1);
    ccr12_info_->offset_l2(lr2);

    ccr12_info_->guess_lambda1();
    ccr12_info_->guess_lambda2();

    iter_end=timer_->get_wall_time();

    std::string headername="Lambda CCSD";
    print_iteration_header_short(headername);

    for (int iter=0;iter<maxiter_;++iter){
     iter_start=iter_end;

     lr1->zero();
     lr2->zero();

     lambda_ccsd_t1->compute_amp(lr1);
     lambda_ccsd_t2->compute_amp(lr2);

     // compute new amplitudes from the residuals
     Ref<Tensor> lambda1_old = info()->lambda1()->copy();
     Ref<Tensor> lambda2_old = info()->lambda2()->copy();
     ccr12_info_->jacobi_lambda1(lr1);
     ccr12_info_->jacobi_lambda2(lr2);
     // compute errors
     Ref<Tensor> lambda1_err = lambda1_old;
     Ref<Tensor> lambda2_err = lambda2_old;
     lambda1_err->daxpy(info()->lambda1(), -1.0);
     lambda2_err->daxpy(info()->lambda2(), -1.0);

     const double lr1norm = RMS(*lambda1_err);
     const double lr2norm = RMS(*lambda2_err);
     double rnorm =std::sqrt(lr1norm*lr1norm+lr2norm*lr2norm);

     iter_end=timer_->get_wall_time();
     print_iteration_short(iter,rnorm,iter_start,iter_end);

     if (rnorm<ccthresh_) break;

     // extrapolate
     l1diis->extrapolate(info()->edata(info()->lambda1()),info()->eerr(lambda1_err));
     l2diis->extrapolate(info()->edata(info()->lambda2()),info()->eerr(lambda2_err));
    }

    print_iteration_footer_short();

    if (perturbative_=="(2)T" or perturbative_=="(2)TQ") {
      timer_->enter("(2)_T correction");
      iter_start=timer_->get_wall_time();

      Ref<CCSD_2T_LEFT>    eval_left=new CCSD_2T_LEFT(info());
      Ref<CCSD_2T_RIGHT>  eval_right=new CCSD_2T_RIGHT(info());
      Ref<Parenthesis2t>     ccsd_2t=new Parenthesis2t(info());
      const double ccsd_2t_correction=ccsd_2t->compute(eval_left,eval_right);
      print_correction(ccsd_2t_correction,energy,"CCSD(2)_T");

      print_timing(timer_->get_wall_time()-iter_start,"(2)_T correction");
      timer_->exit("(2)_T correction");

      if (perturbative_=="(2)TQ") {
       timer_->enter("(2)_Q correction");
       Ref<CCSD_2Q_LEFT>    eval_left_q=new CCSD_2Q_LEFT(info());
       Ref<CCSD_2Q_RIGHT>  eval_right_q=new CCSD_2Q_RIGHT(info());
       Ref<Parenthesis2q>     ccsd_2q=new Parenthesis2q(info());
       const double ccsd_2q_correction=ccsd_2q->compute(eval_left_q,eval_right_q);
       print_correction(ccsd_2q_correction+ccsd_2t_correction,energy,"CCSD(2)_TQ");

       print_timing(timer_->get_wall_time()-iter_start,"(2)_Q correction");
       timer_->exit("(2)_Q correction");
      }
    }

    delete lambda_ccsd_t1;
    delete lambda_ccsd_t2;

  } // end of do_lambda
#endif

  set_energy(energy + this->ref()->energy());
  mem_->sync();
}

