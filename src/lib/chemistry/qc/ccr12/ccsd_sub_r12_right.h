//
// ccsd_sub_r12_right.h --- computes the rhs numerator of the CCSD(2)R12
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

#ifndef _chemistry_qc_ccr12_ccsd_sub_r12_right_h
#define _chemistry_qc_ccr12_ccsd_sub_r12_right_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_SUB_R12_RIGHT : public RefCount { 

  protected:

    CCR12_Info* z;

    std::vector<Tensor*> in; 

    void offset_smith_0_1();
    void smith_0_1_0(); //z->vd2()=>in.at(1)
    void offset_smith_1_1();
    void smith_1_1_0(); //z->f1()=>in.at(2)
    void smith_2_10(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_1(); //z->fy(),in.at(2)=>in.at(1)
    void offset_smith_1_12();
    void smith_2_12(); //z->t1(),z->fy()=>in.at(2)
    void smith_1_12(); //z->v2(),in.at(2)=>in.at(1)
    void smith_0_1(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
    void smith_0_2(Ref<Tensor>&); //z->vd2()=>Ref<Tensor>&
    void offset_smith_0_3();
    void smith_1_3(); //z->t1(),z->fy()=>in.at(1)
    void smith_0_3(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_4();
    void smith_0_4_0(); //z->vd2()=>in.at(1)
    void offset_smith_1_5();
    void smith_2_5(); //z->t1(),z->fy()=>in.at(2)
    void smith_1_5(); //z->v2(),in.at(2)=>in.at(1)
    void smith_1_6(); //z->t1(),z->vd2()=>in.at(1)
    void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
    void offset_smith_0_7();
    void offset_smith_1_7();
    void smith_1_7_0(); //z->v2()=>in.at(2)
    void smith_2_11(); //z->t1(),z->v2()=>in.at(2)
    void smith_1_7(); //z->t2(),in.at(2)=>in.at(1)
    void offset_smith_1_9();
    void offset_smith_2_9();
    void smith_3_9(); //z->t1(),z->v2()=>in.at(3)
    void smith_2_9(); //z->t1(),in.at(3)=>in.at(2)
    void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
    void smith_0_7(Ref<Tensor>&); //z->fy(),in.at(1)=>Ref<Tensor>&

  public:
    CCSD_SUB_R12_RIGHT(CCR12_Info* info);
    ~CCSD_SUB_R12_RIGHT();

    void compute_amp(Ref<Tensor>&);

};



}

#endif


