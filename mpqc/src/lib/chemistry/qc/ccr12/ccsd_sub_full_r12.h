//
// ccsd_sub_full_r12.h : Valeev, Phys Chem Chem Phys 10, 106 (2008)
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
#ifndef __chemistry_qc_ccr12_ccsd_sub_full_r12_h
#define __chemistry_qc_ccr12_ccsd_sub_full_r12_h

#include <chemistry/qc/ccr12/ccsd_sub_r12.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_Sub_Full_R12 : public CCSD_Sub_R12 {
  protected:
    Ref<Tensor> tildeV_dagger_; // left_hand_side 

  public:
    CCSD_Sub_Full_R12(CCR12_Info* inz, Ref<Tensor> right, Ref<Tensor> left)
     : CCSD_Sub_R12(inz, false) {
      tildeV_ = right; 
      tildeV_dagger_ = left;
    };

    ~CCSD_Sub_Full_R12() {};

    double compute() {
      denom_contraction();
      z->prod_iiii(tildeV_dagger_, intermediate_, energy_, false); 
      return z->get_e(energy_);
    };

};

}

#endif

