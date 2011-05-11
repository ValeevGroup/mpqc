//
// ccsd_sub_bar_r12.h : Torheyden & Valeev, Phys Chem Chem Phys 10, 3410 (2008)
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

#pragma once
#ifndef __chemistry_qc_ccr12_ccsd_sub_bar_r12_h
#define __chemistry_qc_ccr12_ccsd_sub_bar_r12_h

#include <chemistry/qc/ccr12/ccsd_sub_r12.h>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <iostream>

namespace sc {

class CCSD_Sub_Bar_R12 : public CCSD_Sub_R12 {
  protected:

    void compute_amp();
    void smith_0_1(Ref<Tensor>&); 
    void smith_0_2(Ref<Tensor>&); 

  public:
    CCSD_Sub_Bar_R12(CCR12_Info* inz) : CCSD_Sub_R12(inz, true) { };

    ~CCSD_Sub_Bar_R12() {};

    double compute() {
      compute_amp();
      if (!z->r12world()->r12tech()->ansatz()->diag()) {
        // ijkl ansatz
        z->denom_contraction(tildeV_, intermediate_);
        z->prod_iiii(tildeV_, intermediate_, energy_, true); 
        return z->get_e(energy_);
      } else {
        // fixed diagonal ansatz
        const Ref<Tensor> gt2 = z->gt2();
        intermediate_->zero();
        z->denom_contraction(gt2, intermediate_);
        z->prod_iiii(gt2, intermediate_, energy_, true);
        const double bterm = z->get_e(energy_);
        energy_->zero();
        z->prod_iiii(tildeV_, gt2, energy_, false); 
        const double direct_en = 2.0*(z->get_e(energy_));
        return direct_en - bterm;
      }
    };

};

}

#endif

