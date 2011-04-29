//
// ccsdtq.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
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
#include <chemistry/qc/ccr12/ccsdtq.h>
#include <chemistry/qc/ccr12/ccsdt_t1.h>
#include <chemistry/qc/ccr12/ccsdtq_t2.h>
#include <chemistry/qc/ccr12/ccsdtq_t3.h>
#include <chemistry/qc/ccr12/ccsdtq_t4.h>
#include <chemistry/qc/ccr12/ccsd_e.h>
#include <chemistry/qc/ccr12/tensorextrap.h>

using std::vector;
using namespace std;
using namespace sc;



/*--------------------------------
  CCR12
 --------------------------------*/
static ClassDesc CCSDTQ_cd(
  typeid(CCSDTQ), "CCSDTQ",1, "public CCR12",
  0, create<CCSDTQ>, create<CCSDTQ>);


CCSDTQ::CCSDTQ(StateIn& s): CCR12(s){
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructer not yet implemented", __FILE__, __LINE__);
}


CCSDTQ::CCSDTQ(const Ref<KeyVal>& keyval): CCR12(keyval) {
  string theory("CCSDTQ");
  theory_ = theory;
  print_theory();
}


CCSDTQ::~CCSDTQ(){
}


void CCSDTQ::compute(){

  CCR12::compute();

  Ref<Tensor> e0 = new Tensor("e", mem_);
  ccr12_info_->offset_e(e0);

  Ref<Tensor> r1 = new Tensor("r1", mem_);
  ccr12_info_->offset_t1(r1, false);
  Ref<Tensor> r2 = new Tensor("r2", mem_);
  ccr12_info_->offset_t2(r2, false);
  Ref<Tensor> r3 = new Tensor("r3", mem_);
  ccr12_info_->offset_t3(r3, false);
  Ref<Tensor> r4 = new Tensor("r4", mem_);
  ccr12_info_->offset_t4(r4, false);

  CCSD_E*    ccsdtq_e  = new CCSD_E(   info());
  CCSDT_T1*  ccsdtq_t1 = new CCSDT_T1( info());
  CCSDTQ_T2* ccsdtq_t2 = new CCSDTQ_T2(info());
  CCSDTQ_T3* ccsdtq_t3 = new CCSDTQ_T3(info());
  CCSDTQ_T4* ccsdtq_t4 = new CCSDTQ_T4(info());

  string theory_ = "CCSDTQ";
  print_iteration_header(theory_);


  timer_->enter("iterations");
  double iter_start = 0.0;
  double iter_end   = timer_->get_wall_time();
  double energy;

  Ref<DIIS> t1diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);
  Ref<DIIS> t2diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);
  Ref<DIIS> t3diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);
  Ref<DIIS> t4diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);

  for (int iter = 0; iter < maxiter_; ++iter){
    iter_start = iter_end;

    e0->zero();
    r1->zero();
    r2->zero();
    r3->zero();
    r4->zero();

    ccsdtq_e->compute_amp(e0);
    ccsdtq_t1->compute_amp(r1);
    ccsdtq_t2->compute_amp(r2);
    ccsdtq_t3->compute_amp(r3);
    ccsdtq_t4->compute_amp(r4);

    // compute new amplitudes from the residuals
    Ref<Tensor> t1_old = info()->t1()->copy();
    Ref<Tensor> t2_old = info()->t2()->copy();
    Ref<Tensor> t3_old = info()->t3()->copy();
    Ref<Tensor> t4_old = info()->t4()->copy();
    ccr12_info_->jacobi_t1(r1);
    ccr12_info_->jacobi_t2(r2);
    ccr12_info_->jacobi_t3(r3);
    ccr12_info_->jacobi_t4(r4);
    // compute errors
    Ref<Tensor> t1_err = t1_old;
    Ref<Tensor> t2_err = t2_old;
    Ref<Tensor> t3_err = t3_old;
    Ref<Tensor> t4_err = t4_old;
    t1_err->daxpy(info()->t1(), -1.0);
    t2_err->daxpy(info()->t2(), -1.0);
    t3_err->daxpy(info()->t3(), -1.0);
    t4_err->daxpy(info()->t4(), -1.0);

    const double r1norm = RMS(*t1_err);
    const double r2norm = RMS(*t2_err);
    const double r3norm = RMS(*t3_err);
    const double r4norm = RMS(*t4_err);
    const double rnorm = std::sqrt(r1norm * r1norm
                                 + r2norm * r2norm
                                 + r3norm * r3norm
                                 + r4norm * r4norm);

    energy = ccr12_info_->get_e(e0);

    iter_end = timer_->get_wall_time();
    print_iteration(iter, energy, rnorm, iter_start, iter_end);

    if (rnorm < ccthresh_) break;

    // extrapolate
    t1diis->extrapolate(info()->edata(info()->t1()), info()->eerr(t1_err));
    t2diis->extrapolate(info()->edata(info()->t2()), info()->eerr(t2_err));
    t3diis->extrapolate(info()->edata(info()->t3()), info()->eerr(t3_err));
    t4diis->extrapolate(info()->edata(info()->t4()), info()->eerr(t4_err));

  }
  timer_->exit("iterations");

  print_iteration_footer();


  delete ccsdtq_t4;
  delete ccsdtq_t3;
  delete ccsdtq_t2;
  delete ccsdtq_t1;
  delete ccsdtq_e;

  set_energy(energy + this->ref()->energy());
  mem_->sync();
}

