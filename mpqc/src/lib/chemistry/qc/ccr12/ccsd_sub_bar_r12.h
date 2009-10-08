//
// ccsd_sub_bar_r12.h : Torheyden & Valeev, Phys Chem Chem Phys 10, 106 (2008)
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

#pragma once
#ifndef __chemistry_qc_ccr12_ccsd_sub_bar_r12_h
#define __chemistry_qc_ccr12_ccsd_sub_bar_r12_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_Sub_Bar_R12 : public RefCount {
  protected:
    CCR12_Info* z;
    Ref<Tensor> tildeV_;
    Ref<Tensor> intermediate_;
    Ref<Tensor> energy_;

    void compute_amp();
    void smith_0_1(Ref<Tensor>&); 
    void smith_0_2(Ref<Tensor>&); 
    void denom_contraction();

  public:
    CCSD_Sub_Bar_R12(CCR12_Info* inz) : z(inz) {
      tildeV_ = new Tensor("tildeV", z->mem()); 
      intermediate_ = new Tensor("intermediate", z->mem()); 
      energy_ = new Tensor("energy", z->mem()); 
      z->offset_gt2(tildeV_, false);
      z->offset_gt2(intermediate_, false);
      z->offset_e(energy_);
    };

    ~CCSD_Sub_Bar_R12() {};

    double compute() {
      compute_amp();
#if 1
      denom_contraction();
#else    // just for checking...
      intermediate_ = tildeV_;
#endif
      z->prod_iiii(tildeV_, intermediate_, energy_, true); 
      return z->get_e(energy_);
    };

};

}

#endif

