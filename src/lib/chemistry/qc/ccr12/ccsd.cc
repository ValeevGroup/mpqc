//
// ccsd.cc
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

#include <util/misc/regtime.h>
#include <math/optimize/diis.h>
#include <chemistry/qc/ccr12/ccsd.h>
#include <chemistry/qc/ccr12/ccsd_e.h>
#include <chemistry/qc/ccr12/ccsd_t1.h>
#include <chemistry/qc/ccr12/ccsd_t2.h>
#include <chemistry/qc/ccr12/ccsd_pt.h>
#include <chemistry/qc/ccr12/ccsd_pt_left.h>
#include <chemistry/qc/ccr12/ccsd_pt_right.h>
#include <chemistry/qc/ccr12/lambda_ccsd_t1.h>
#include <chemistry/qc/ccr12/lambda_ccsd_t2.h>
#include <chemistry/qc/ccr12/parenthesis2t.h>
#include <chemistry/qc/ccr12/parenthesis2q.h>
#include <chemistry/qc/ccr12/ccsd_2t_left.h>
#include <chemistry/qc/ccr12/ccsd_2t_right.h>
#include <chemistry/qc/ccr12/ccsd_2q_left.h>
#include <chemistry/qc/ccr12/ccsd_2q_right.h>
#include <chemistry/qc/ccr12/ccsd_sub_r12_left.h>
#include <chemistry/qc/ccr12/ccsd_sub_r12_right.h>
#include <chemistry/qc/ccr12/ccr12_triples.h>
#include <chemistry/qc/ccr12/ccsd_sub_bar_r12.h> 
#include <chemistry/qc/ccr12/ccsd_sub_full_r12.h>

using namespace std;
using namespace sc;


static ClassDesc CCSD_cd(
  typeid(CCSD), "CCSD", 1, "public CCR12",
  0, create<CCSD>, create<CCSD>);


CCSD::CCSD(StateIn& s): CCR12(s){
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructer not yet implemented", __FILE__, __LINE__);
}


CCSD::CCSD(const Ref<KeyVal>& keyval): CCR12(keyval){
  string theory("CCSD");
  theory_ = theory;

  perturbative_ = keyval->stringvalue("perturbative", KeyValValuestring(""));
  ostringstream oss; oss.setf(ios_base::uppercase);
  oss << perturbative_;
  perturbative_ = oss.str();

  print_theory();

  // validate user input
  if (perturbative_ == "(2)R12FULL" || perturbative_ == "(2)TQR12" || perturbative_ == "(T)R12[DT]") {
    // diagonal ansatz not implemented
    if (this->r12world()->r12tech()->ansatz()->diag())
      throw FeatureNotImplemented("diagonal ansatz for perturbative R12 corrections in SMITH-based code (use Psi-based code)",
                                  __FILE__, __LINE__, this->class_desc());
  }
}


CCSD::~CCSD(){
}


void CCSD::compute(){

  CCR12::compute();

  Ref<Tensor> e0 = new Tensor("e", mem_);
  ccr12_info_->offset_e(e0);

  Ref<Tensor> r1 = new Tensor("r1", mem_);
  ccr12_info_->offset_t1(r1, false);
  Ref<Tensor> r2 = new Tensor("r2", mem_);
  ccr12_info_->offset_t2(r2, false);

  CCSD_E*  ccsd_e  = new CCSD_E( info());
  CCSD_T1* ccsd_t1 = new CCSD_T1(info());
  CCSD_T2* ccsd_t2 = new CCSD_T2(info());

  string theory_ = "CCSD";
  print_iteration_header(theory_);


  timer_->enter("CC iterations");
  double iter_start = 0.0;
  double iter_end   = timer_->get_wall_time();
  double energy;


  Ref<DIIS> t1diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);
  Ref<DIIS> t2diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);

  for (int iter = 0; iter < maxiter_; ++iter){
    iter_start = iter_end;

    e0->zero();
    r1->zero();
    r2->zero();

    ccsd_t1->compute_amp(r1);
    ccsd_t2->compute_amp(r2);
    ccsd_e->compute_amp(e0);

    energy = ccr12_info_->get_e(e0) + ccr12_info_->t1()->ddot(r1) + ccr12_info_->energy_lagrangian_r2(r2);

    // compute new amplitudes from the residuals
    Ref<Tensor> t1_old = info()->t1()->copy();
    Ref<Tensor> t2_old = info()->t2()->copy();
    ccr12_info_->jacobi_t1(r1);
    ccr12_info_->jacobi_t2(r2);
    // compute errors
    Ref<Tensor> t1_err = t1_old->copy();
    Ref<Tensor> t2_err = t2_old->copy();
    t1_err->daxpy(info()->t1(), -1.0);
    t2_err->daxpy(info()->t2(), -1.0);

    const double r1norm = RMS(*t1_err);
    const double r2norm = RMS(*t2_err);
    const double rnorm = std::sqrt(r1norm * r1norm + r2norm * r2norm);

    iter_end = timer_->get_wall_time();
    print_iteration(iter, energy, rnorm, iter_start, iter_end);

#if 0
    ExEnv::out0() << scprintf("t1norm = %20.15lf t2norm = %20.15lf", t1_old->norm(), t2_old->norm()) << std::endl;
    ExEnv::out0() << scprintf("r1norm = %20.15lf r2norm = %20.15lf", r1norm, r2norm) << std::endl;
#endif

    // done? break free
    if (rnorm < ccthresh_) break;

    // extrapolate
    t1diis->extrapolate(info()->edata(info()->t1()), info()->eerr(t1_err));
    t2diis->extrapolate(info()->edata(info()->t2()), info()->eerr(t2_err));



  }
  timer_->exit("CC iterations");

  print_iteration_footer();

  delete ccsd_t2;
  delete ccsd_t1;
  delete ccsd_e;


  // perturbative treatment of geminal doubles 
  if (perturbative_ == "(T)R12" || perturbative_ == "(T)R12[DT]" || perturbative_ == "(2)R12") {
    timer_->enter("R12 doubles correction");
    iter_start = timer_->get_wall_time();

    Ref<CCSD_Sub_Bar_R12> ccsd_sub_bar_r12 = new CCSD_Sub_Bar_R12(info()); 
    const double doubles_correction = ccsd_sub_bar_r12->compute();
    print_correction(doubles_correction, energy, "R12 doubles");

    print_timing(timer_->get_wall_time() - iter_start, "R12 doubles correction");
    timer_->exit("R12 doubles correction");
    energy += doubles_correction;

  }

  // conventional perturbative triples
  if (perturbative_ == "(T)" || perturbative_ == "(T)R12" || perturbative_ == "(T)R12[DT]") {
    timer_->enter("(T) correction");
    iter_start = timer_->get_wall_time();

    // ccsd(t) calls specific driver, because it is a bit
    // different from others (i.e. (2)t etc)
    Ref<CCSD_PT_LEFT> eval_left = new CCSD_PT_LEFT(info());
    Ref<CCSD_PT_RIGHT> eval_right = new CCSD_PT_RIGHT(info());
    Ref<CCSD_PT> ccsd_pt = new CCSD_PT(info());
    const double ccsd_pt_correction = ccsd_pt->compute_energy(eval_left,eval_right);
    print_correction(ccsd_pt_correction, energy, "CCSD(T)");

    print_timing(timer_->get_wall_time() - iter_start, "(T) correction");
    timer_->exit("(T) correction");
    energy += ccsd_pt_correction;
  }


  // the CCSD(T)_\mathrm{R12[DT]} correction comes here...
  if (perturbative_ == "(T)R12[DT]") {
    timer_->enter("R12 triples correction");
    iter_start = timer_->get_wall_time();

    if (!info()->restricted())
      throw ProgrammingError("Triples corrections are implemented only for RHF reference so far.", __FILE__, __LINE__);

    Ref<CCR12_Triples> triples = new CCR12_Triples(info());
    const double triples_correction = triples->compute();
    print_correction(triples_correction, energy, "R12 triples");

    print_timing(timer_->get_wall_time() - iter_start, "R12 triples correction");
    timer_->exit("R12 triples correction");
    energy += triples_correction;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////


  bool do_lambda = false; // will be judeged by input keywords
  // more will come; e.g. dipole, etc
  if (perturbative_ == "(2)T" || perturbative_ == "(2)TQ"
    || perturbative_ == "(2)TQR12" || perturbative_ == "(2)R12FULL") do_lambda = true;

  if (do_lambda) {
    Ref<DIIS> l1diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);
    Ref<DIIS> l2diis = new DIIS(diis_start_, ndiis_, 0.005, 3, 1, 0.0);

    LAMBDA_CCSD_T1* lambda_ccsd_t1 = new LAMBDA_CCSD_T1(info());
    LAMBDA_CCSD_T2* lambda_ccsd_t2 = new LAMBDA_CCSD_T2(info());

    Ref<Tensor> lr1 = new Tensor("lr1", mem_);
    Ref<Tensor> lr2 = new Tensor("lr2", mem_);
    ccr12_info_->offset_l1(lr1);
    ccr12_info_->offset_l2(lr2);

    ccr12_info_->guess_lambda1();
    ccr12_info_->guess_lambda2();

    iter_end=timer_->get_wall_time();

    string headername = "Lambda CCSD";
    print_iteration_header_short(headername);

    for (int iter = 0; iter < maxiter_; ++iter){
     iter_start = iter_end;

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
     double rnorm = ::sqrt(lr1norm * lr1norm + lr2norm * lr2norm);

     iter_end = timer_->get_wall_time();
     print_iteration_short(iter, rnorm, iter_start, iter_end);

     if (rnorm < ccthresh_) break;

     // extrapolate
     l1diis->extrapolate(info()->edata(info()->lambda1()), info()->eerr(lambda1_err));
     l2diis->extrapolate(info()->edata(info()->lambda2()), info()->eerr(lambda2_err));
    }

    print_iteration_footer_short();

    if (perturbative_ == "(2)T" || perturbative_ == "(2)TQ" || perturbative_ == "(2)TQR12") {
      timer_->enter("(2)T correction");
      iter_start = timer_->get_wall_time();

      Ref<CCSD_2T_LEFT>   eval_left = new CCSD_2T_LEFT(info());
      Ref<CCSD_2T_RIGHT> eval_right = new CCSD_2T_RIGHT(info());
      Ref<Parenthesis2t>    ccsd_2t = new Parenthesis2t(info());
      const double ccsd_2t_correction = ccsd_2t->compute_energy(eval_left, eval_right);
      print_correction(ccsd_2t_correction, energy, "CCSD(2)T");

      print_timing(timer_->get_wall_time() - iter_start, "(2)T correction");
      timer_->exit("(2)T correction");

      energy += ccsd_2t_correction;

      if (perturbative_ == "(2)TQ" || perturbative_ == "(2)TQR12") {
        timer_->enter("(2)Q correction");
        Ref<CCSD_2Q_LEFT>   eval_left_q = new CCSD_2Q_LEFT(info());
        Ref<CCSD_2Q_RIGHT> eval_right_q = new CCSD_2Q_RIGHT(info());
        Ref<Parenthesis2q>      ccsd_2q = new Parenthesis2q(info());
        const double ccsd_2q_correction = ccsd_2q->compute_energy(eval_left_q,eval_right_q);
        print_correction(ccsd_2t_correction + ccsd_2q_correction, energy - ccsd_2t_correction, "CCSD(2)TQ");

        print_timing(timer_->get_wall_time() - iter_start, "(2)Q correction");
        timer_->exit("(2)Q correction");

        energy += ccsd_2q_correction;
      }

    }
    if (perturbative_ == "(2)R12FULL" || perturbative_ == "(2)TQR12") {

      timer_->enter("(2)R12 correction (full)");
      iter_start = timer_->get_wall_time();

      Ref<CCSD_SUB_R12_RIGHT> eval_right = new CCSD_SUB_R12_RIGHT(info());
      Ref<CCSD_SUB_R12_LEFT> eval_left = new CCSD_SUB_R12_LEFT(info());
      Ref<Tensor> num_right = new Tensor("num_right", mem_);
      Ref<Tensor> num_left = new Tensor("num_left", mem_);
      ccr12_info_->offset_gt2(num_right, false);
      ccr12_info_->offset_gt2(num_left, false);

      eval_right->compute_amp(num_right);
      eval_left->compute_amp(num_left);

      // then evaluate the energy contribution
      Ref<CCSD_Sub_Full_R12> ccsd_sub_full_r12 = new CCSD_Sub_Full_R12(info(), num_right, num_left);
      const double ccsd_sub_r12_correction = ccsd_sub_full_r12->compute();
      print_correction(ccsd_sub_r12_correction, energy, "CCSD(2)R12");

      print_timing(timer_->get_wall_time() - iter_start, "(2)R12 correction (full)");
      timer_->exit("(2)R12 correction (full)");
      energy += ccsd_sub_r12_correction;
    }

    delete lambda_ccsd_t1;
    delete lambda_ccsd_t2;

  } // end of do_lambda

  set_energy(energy + this->ref()->energy());
  mem_->sync();
}

