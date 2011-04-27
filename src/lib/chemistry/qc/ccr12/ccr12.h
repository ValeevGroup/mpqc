//
// ccr12.h -- base class for CC and CC-R12 classes
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
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

#ifndef _chemistry_qc_ccr12_ccr12_h
#define _chemistry_qc_ccr12_ccr12_h

#include <string>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

/** CCR12 is the base class for CC and CC-R12 methods.  */
class CCR12: public Wavefunction {
  protected:
    Ref<ThreadGrp> thrgrp_;
    Ref<MessageGrp> msggrp_;
    Ref<MemoryGrp> mem_;

    Ref<R12IntEval> r12eval_;              // the R12 intermediates evaluator
    Ref<R12WavefunctionWorld> r12world_;   // parameters for r12eval_
    Ref<SCF> reference_;

    CCR12_Info* ccr12_info_;
    int nfzc_, nfzv_;
    long worksize_;
    long memorysize_;
    bool rhf_;
    int maxiter_;
    Ref<RegionTimer> timer_;
    std::string theory_;
    std::string perturbative_;
    double ccthresh_;

    int ndiis_;
    int diis_start_;

  public:
    CCR12(StateIn&);
    CCR12(const Ref<KeyVal>&);
    ~CCR12();
    CCR12_Info* info() const { return ccr12_info_;};

    const Ref<R12WavefunctionWorld>& r12world() const { return r12world_; }
    const Ref<R12IntEval>& r12eval() const { return r12eval_; }
    Ref<SCF> ref() { return reference_; }

    int spin_polarized() { return ccr12_info_->restricted(); };

    RefSymmSCMatrix density() { return 0; };
    int nelectron() { return ccr12_info_->naoa()+ccr12_info_->naob(); };
    void obsolete();
    int value_implemented() const { return 1; };

  protected:
    void compute();
    static double ref_to_ccr12_acc() { return 100.0; }

    // print utilities
    void print_theory();
    void print_iteration_header(std::string);
    void print_iteration_footer();
    void print_iteration(int,double,double,double,double);
    void print_iteration_header_short(std::string);
    void print_iteration_footer_short();
    void print_iteration_short(int,double,double,double);
    void print_correction(double,double,std::string);
    void print(std::ostream&);
    void print_timing(double,std::string);

};

}

#endif
