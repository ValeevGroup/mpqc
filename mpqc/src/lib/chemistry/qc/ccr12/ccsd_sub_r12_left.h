//
// ccsd_sub_r12_left.h --- computes the lhs numerator of the CCSD(2)R12
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

#ifndef _chemistry_qc_ccr12_ccsd_sub_r12_left_h
#define _chemistry_qc_ccr12_ccsd_sub_r12_left_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_SUB_R12_LEFT : public RefCount { 

  protected:

    CCR12_Info* z;

    std::vector<Tensor*> in; 

    void smith_0_1(Ref<Tensor>&); //z->vr2()=>Ref<Tensor>&
    void offset_smith_0_2();
    void smith_0_2_0(); //z->vr2()=>in.at(1)
    void offset_smith_1_2();
    void smith_1_2_0(); //z->f1()=>in.at(2)
    void smith_2_5(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_2(); //z->fx(),in.at(2)=>in.at(1)
    void smith_0_2(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_3();
    void smith_1_3(); //z->l1(),z->fx()=>in.at(1)
    void smith_0_3(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_6();
    void offset_smith_1_6();
    void smith_2_6(); //z->l1(),z->fx()=>in.at(2)
    void smith_1_6(); //z->t1(),in.at(2)=>in.at(1)
    void offset_smith_1_21();
    void smith_2_21(); //z->t2(),z->l2()=>in.at(2)
    void smith_1_21(); //z->fx(),in.at(2)=>in.at(1)
    void smith_0_6(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_7();
    void smith_1_7(); //z->t1(),z->l1()=>in.at(1)
    void smith_1_23(); //z->t2(),z->l2()=>in.at(1)
    void smith_0_7(Ref<Tensor>&); //z->vr2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_8();
    void smith_0_8_0(); //z->vr2()=>in.at(1)
    void offset_smith_1_8();
    void smith_1_8_0(); //z->f1()=>in.at(2)
    void smith_2_12(); //z->t1(),z->v2()=>in.at(2)
    void smith_2_19(); //z->t2(),z->v2()=>in.at(2)
    void smith_1_8(); //z->fx(),in.at(2)=>in.at(1)
    void offset_smith_1_9();
    void offset_smith_2_9();
    void smith_2_9_0(); //z->f1()=>in.at(3)
    void smith_3_16(); //z->t1(),z->v2()=>in.at(3)
    void smith_2_9(); //z->fx(),in.at(3)=>in.at(2)
    void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
    void smith_0_8(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_10();
    void offset_smith_1_10();
    void smith_1_10_0(); //z->v2()=>in.at(2)
    void smith_2_14(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_10(); //z->l2(),in.at(2)=>in.at(1)
    void offset_smith_1_13();
    void smith_2_13(); //z->t1(),z->l2()=>in.at(2)
    void smith_1_13(); //z->v2(),in.at(2)=>in.at(1)
    void offset_smith_1_17();
    void offset_smith_2_17();
    void smith_3_17(); //z->t1(),z->l2()=>in.at(3)
    void smith_2_17(); //z->t1(),in.at(3)=>in.at(2)
    void smith_2_20(); //z->t2(),z->l2()=>in.at(2)
    void smith_1_17(); //z->v2(),in.at(2)=>in.at(1)
    void smith_0_10(Ref<Tensor>&); //z->fx(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_15();
    void smith_1_15(); //z->t1(),z->l2()=>in.at(1)
    void smith_0_15(Ref<Tensor>&); //z->vr2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_18();
    void offset_smith_1_18();
    void smith_2_18(); //z->t1(),z->l2()=>in.at(2)
    void smith_1_18(); //z->t1(),in.at(2)=>in.at(1)
    void smith_1_22(); //z->t2(),z->l2()=>in.at(1)
    void smith_0_18(Ref<Tensor>&); //z->vr2(),in.at(1)=>Ref<Tensor>&

  public:
    CCSD_SUB_R12_LEFT(CCR12_Info* info);
    ~CCSD_SUB_R12_LEFT();

    void compute_amp(Ref<Tensor>&);

};



}

#endif


