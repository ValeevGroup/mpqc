//
// ccsdt.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <util/misc/regtime.h>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/optimize/diis.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/ccsdt.h>
#include <chemistry/qc/ccr12/ccsdt_t1.h>
#include <chemistry/qc/ccr12/ccsdt_t2.h>
#include <chemistry/qc/ccr12/ccsdt_t3.h>
#include <chemistry/qc/ccr12/ccsd_e.h>
#include <chemistry/qc/ccr12/tensorextrap.h>

using std::vector;
using namespace std;
using namespace sc;
using namespace sc::LinearR12;


/*--------------------------------
  CCR12
 --------------------------------*/
static ClassDesc CCSDT_cd(
  typeid(CCSDT), "CCSDT", 1, "public CCR12",
  0, create<CCSDT>, create<CCSDT>);


CCSDT::CCSDT(StateIn& s): CCR12(s){
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructer not yet implemented", __FILE__, __LINE__);
}


CCSDT::CCSDT(const Ref<KeyVal>& keyval): CCR12(keyval){

  keyval_ = keyval;
  string theory_ = "CCSDT";
  common_init(theory_);
}


CCSDT::~CCSDT(){
}


void CCSDT::compute(){

  CCR12::compute();

  Ref<Tensor> e0 = new Tensor("e", mem_);
  ccr12_info_->offset_e(e0);

  Ref<Tensor> r1 = new Tensor("r1", mem_);
  ccr12_info_->offset_t1(r1, false);
  Ref<Tensor> r2 = new Tensor("r2", mem_);
  ccr12_info_->offset_t2(r2, false);
  Ref<Tensor> r3 = new Tensor("r3", mem_);
  ccr12_info_->offset_t3(r3, false);

  CCSD_E*   ccsdt_e  = new CCSD_E(  info());
  CCSDT_T1* ccsdt_t1 = new CCSDT_T1(info());
  CCSDT_T2* ccsdt_t2 = new CCSDT_T2(info());
  CCSDT_T3* ccsdt_t3 = new CCSDT_T3(info());

  string theory_ = "CCSDT";
  print_iteration_header(theory_);


  timer_->enter("iterations");
  double iter_start = 0.0;
  double iter_end   = timer_->get_wall_time();
  double energy;

  Ref<DIIS> t1diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);
  Ref<DIIS> t2diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);
  Ref<DIIS> t3diis = new DIIS(diis_start_, ndiis_, 0.05, 2, 1, 0.0);

  for (int iter = 0; iter < maxiter_; ++iter){
    iter_start = iter_end;

    e0->zero();
    r1->zero();
    r2->zero();
    r3->zero();

    ccsdt_e->compute_amp(e0);
    ccsdt_t1->compute_amp(r1);
    ccsdt_t2->compute_amp(r2);
    ccsdt_t3->compute_amp(r3);

    // compute new amplitudes from the residuals
    Ref<Tensor> t1_old = info()->t1()->copy();
    Ref<Tensor> t2_old = info()->t2()->copy();
    Ref<Tensor> t3_old = info()->t3()->copy();
    ccr12_info_->jacobi_t1(r1);
    ccr12_info_->jacobi_t2(r2);
    ccr12_info_->jacobi_t3(r3);
    // compute errors
    Ref<Tensor> t1_err = t1_old;
    Ref<Tensor> t2_err = t2_old;
    Ref<Tensor> t3_err = t3_old;
    t1_err->daxpy(info()->t1(), -1.0);
    t2_err->daxpy(info()->t2(), -1.0);
    t3_err->daxpy(info()->t3(), -1.0);

    const double r1norm = RMS(*t1_err);
    const double r2norm = RMS(*t2_err);
    const double r3norm = RMS(*t3_err);
    const double rnorm =std::sqrt(r1norm * r1norm
                                + r2norm * r2norm
                                + r3norm * r3norm);

    energy=ccr12_info_->get_e(e0);

    iter_end=timer_->get_wall_time();
    print_iteration(iter, energy, rnorm, iter_start, iter_end);

    if (rnorm < ccthresh_) break;

    // extrapolate
    t1diis->extrapolate(info()->edata(info()->t1()), info()->eerr(t1_err));
    t2diis->extrapolate(info()->edata(info()->t2()), info()->eerr(t2_err));
    t3diis->extrapolate(info()->edata(info()->t3()), info()->eerr(t3_err));

  }
  timer_->exit("iterations");

  print_iteration_footer();


  delete ccsdt_t3;
  delete ccsdt_t2;
  delete ccsdt_t1;
  delete ccsdt_e;

  set_energy(energy);
  mem_->sync();
}

