//
// lambda_ccsdpr12_t2.h --- computes the L2 residual vector for Lambda-CCSD(R12)
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

#ifndef _chemistry_qc_ccr12_lambda_ccsdp12_t2_h
#define _chemistry_qc_ccr12_lambda_ccsdp12_t2_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class LAMBDA_CCSDPR12_T2 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;

   void smith_0_1(Ref<Tensor>&); //z->v2()=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->v2()=>in.at(1)
   void smith_1_5(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
   void smith_0_3(Ref<Tensor>&); //z->l1(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_4();
   void smith_1_4(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->l1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_1_6(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_35(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_7();
   void smith_0_7_0(); //z->f1()=>in.at(1)
   void smith_1_15(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_24();
   void smith_2_24(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_24(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_29(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_36(); //z->gt2(),z->vr2()=>in.at(1)
   void smith_0_7(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_8();
   void smith_0_8_0(); //z->f1()=>in.at(1)
   void smith_1_17(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_26();
   void smith_2_26(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_26(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_31(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_38(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_9();
   void smith_0_9_0(); //z->f1()=>in.at(1)
   void smith_1_21(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_10();
   void smith_0_10_0(); //z->v2()=>in.at(1)
   void offset_smith_1_16();
   void smith_1_16_0(); //z->v2()=>in.at(2)
   void smith_2_25(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_16(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_30(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_37(); //z->gt2(),z->vr2()=>in.at(1)
   void smith_0_10(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_11();
   void smith_0_11_0(); //z->v2()=>in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_27();
   void smith_2_27(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_27(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_39(); //z->qy(),z->v2()=>in.at(1)
   void smith_0_11(Ref<Tensor>&); //z->l2(),in.at(1)=>Ref<Tensor>&
   void smith_0_12(Ref<Tensor>&); //z->l2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>in.at(1)
   void smith_1_22(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->ly(),in.at(1)=>Ref<Tensor>&
   void smith_0_14(Ref<Tensor>&); //z->glambda2(),z->vd2()=>Ref<Tensor>&
   void offset_smith_0_18();
   void smith_1_18(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_18(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_20();
   void smith_1_20(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_20(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_23();
   void smith_1_23(); //z->t1(),z->ly()=>in.at(1)
   void smith_0_23(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_28();
   void offset_smith_1_28();
   void smith_2_28(); //z->t1(),z->l2()=>in.at(2)
   void smith_1_28(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_34(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_28(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_32();
   void smith_1_32(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_32(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_33();
   void smith_1_33(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_33(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_40();
   void smith_1_40(); //z->l2(),z->qy()=>in.at(1)
   void smith_0_40(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&

  public:
   LAMBDA_CCSDPR12_T2(CCR12_Info* info);
   ~LAMBDA_CCSDPR12_T2();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


