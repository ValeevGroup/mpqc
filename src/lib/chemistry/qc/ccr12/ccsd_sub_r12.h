//
// ccsd_sub_r12.h : Base class for (2)R12 family of methods 
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

#pragma once
#ifndef __chemistry_qc_ccr12_ccsd_sub_r12_h
#define __chemistry_qc_ccr12_ccsd_sub_r12_h

#include <chemistry/qc/ccr12/ccr12_info.h>
#include <cassert>

namespace sc {

/** CCSD_Sub_R12 is the base class for some (2)R12 methods */
class CCSD_Sub_R12 : public RefCount {

  protected:
    CCR12_Info* z;
    Ref<Tensor> tildeV_;
    Ref<Tensor> intermediate_;
    Ref<Tensor> energy_;

    void denom_contraction() { z->denom_contraction(tildeV_, intermediate_); };

  public:
    CCSD_Sub_R12(CCR12_Info* inz, const bool do_init) : z(inz) {
      intermediate_ = new Tensor("intermediate", z->mem()); 
      energy_ = new Tensor("energy", z->mem()); 
      if (do_init) {
        tildeV_ = new Tensor("tildeV", z->mem()); 
        z->offset_gt2(tildeV_, false);
      }
      z->offset_gt2(intermediate_, false);
      z->offset_e(energy_);
    };

    virtual ~CCSD_Sub_R12() {};
    virtual double compute() = 0;
};

}

#endif

