//
// lambda_ccsdpr12_t1.h --- computes the L1 residual vector for Lambda-CCSD(R12)
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

#ifndef _chemistry_qc_ccr12_lambda_ccsdp12_t1_h
#define _chemistry_qc_ccr12_lambda_ccsdp12_t1_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class LAMBDA_CCSDPR12_T1 {

  protected:

   CCR12_Info* z;

   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->l1()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->l2()=>kn.at(2)
   void offset_k3();
   void smith_k3(); //z->t1(),z->ly()=>kn.at(3)
   void offset_k4();
   void smith_k4(); //z->t2(),z->l2()=>kn.at(4)
   void offset_k5();
   void smith_k5(); //z->t2(),z->l2()=>kn.at(5)
   void offset_k6();
   void smith_k6(); //z->t1(),z->v2()=>kn.at(6)
   void offset_k7();
   void smith_k7(); //z->t1(),kn.at(2)=>kn.at(7)
   void smith_0_1(Ref<Tensor>&); //kn.at(0)=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_5(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_9(); //z->t1(),kn.at(0)=>in.at(1)
   void smith_1_12(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_15(); //z->c2(),z->vr2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_3();
   void smith_0_3_0(); //z->f1()=>in.at(1)
   void smith_1_6(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_10(); //z->t1(),kn.at(0)=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
   void smith_0_4(Ref<Tensor>&); //z->l1(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_7();
   void smith_0_7_0(); //kn.at(1)=>in.at(1)
   void smith_0_7_1(); //kn.at(5)=>in.at(1)
   void smith_0_7(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_8();
   void smith_1_8(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_38(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_11();
   void offset_smith_1_11();
   void smith_1_11_0(); //kn.at(1)=>in.at(2)
   void smith_1_11_1(); //kn.at(5)=>in.at(2)
   void smith_1_11(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_14(); //z->t2(),z->l1()=>in.at(1)
   void smith_1_56(); //z->t2(),kn.at(2)=>in.at(1)
   void smith_0_11(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_13();
   void smith_1_13(); //z->t2(),z->l1()=>in.at(1)
   void offset_smith_1_49();
   void smith_1_49_0(); //kn.at(7)=>in.at(2)
   void smith_1_49_1(); //kn.at(4)=>in.at(2)
   void smith_1_49(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_53(); //z->t2(),kn.at(2)=>in.at(1)
   void smith_1_55(); //z->t1(),kn.at(5)=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_16();
   void smith_1_16(); //z->l1(),z->qy()=>in.at(1)
   void smith_1_61(); //z->qy(),kn.at(2)=>in.at(1)
   void smith_0_16(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_17();
   void smith_1_17(); //z->l1(),z->qy()=>in.at(1)
   void smith_1_62(); //z->qy(),kn.at(2)=>in.at(1)
   void smith_0_17(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_18();
   void smith_0_18_0(); //z->v2()=>in.at(1)
   void offset_smith_1_23();
   void smith_1_23_0(); //z->v2()=>in.at(2)
   void smith_2_30(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_23(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_29();
   void offset_smith_2_29();
   void smith_2_29_0(); //z->v2()=>in.at(3)
   void smith_2_29_1(); //kn.at(6)=>in.at(3)
   void smith_2_29(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_50(); //z->t2(),z->v2()=>in.at(2)
   void smith_2_58(); //z->c2(),z->vr2()=>in.at(2)
   void smith_1_29(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_33();
   void smith_1_33_0(); //z->v2()=>in.at(2)
   void smith_1_33_1(); //kn.at(6)=>in.at(2)
   void smith_1_33(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_34(); //z->t2(),z->v2()=>in.at(1)
   void offset_smith_1_39();
   void smith_1_39_0(); //z->v2()=>in.at(2)
   void smith_2_59(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_39(); //z->qy(),in.at(2)=>in.at(1)
   void smith_1_40(); //z->c2(),z->vr2()=>in.at(1)
   void smith_1_52(); //z->t2(),kn.at(0)=>in.at(1)
   void offset_smith_1_60();
   void smith_2_60(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_60(); //z->qy(),in.at(2)=>in.at(1)
   void smith_0_18(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void smith_0_19(Ref<Tensor>&); //z->l2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_20();
   void smith_0_20_0(); //z->v2()=>in.at(1)
   void offset_smith_1_26();
   void smith_1_26_0(); //z->v2()=>in.at(2)
   void smith_2_43(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_26(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_45(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_20(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void smith_0_21(Ref<Tensor>&); //z->lc2(),z->vd2()=>Ref<Tensor>&
   void smith_0_22(Ref<Tensor>&); //z->v2(),kn.at(2)=>Ref<Tensor>&
   void smith_0_24(Ref<Tensor>&); //z->v2(),kn.at(2)=>Ref<Tensor>&
   void offset_smith_0_25();
   void smith_1_25(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_25(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void smith_0_27(Ref<Tensor>&); //z->v2(),kn.at(3)=>Ref<Tensor>&
   void offset_smith_0_28();
   void smith_1_28(); //z->t1(),z->lc2()=>in.at(1)
   void smith_0_28(Ref<Tensor>&); //z->vd2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_31();
   void smith_0_31_0(); //kn.at(7)=>in.at(1)
   void smith_0_31_1(); //kn.at(4)=>in.at(1)
   void smith_0_31(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_32();
   void smith_1_32(); //z->t1(),kn.at(2)=>in.at(1)
   void smith_1_36(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_32(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_41();
   void smith_1_41(); //z->l2(),z->qy()=>in.at(1)
   void smith_0_41(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_42();
   void smith_1_42(); //z->l2(),z->qy()=>in.at(1)
   void smith_0_42(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_44();
   void smith_1_44(); //z->t1(),kn.at(3)=>in.at(1)
   void smith_1_46(); //z->t2(),z->ly()=>in.at(1)
   void smith_0_44(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_47();
   void smith_1_47(); //z->t2(),z->ly()=>in.at(1)
   void smith_0_47(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&

  public:
   LAMBDA_CCSDPR12_T1(CCR12_Info* info);
   ~LAMBDA_CCSDPR12_T1();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


