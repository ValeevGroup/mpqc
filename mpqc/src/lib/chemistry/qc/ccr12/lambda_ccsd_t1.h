//
// lambda_ccsd_t1.h --- computes the L1 residual vector for Lambda-CCSD
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

#ifndef _chemistry_qc_ccr12_lambda_ccsd_t1_h
#define _chemistry_qc_ccr12_lambda_ccsd_t1_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class LAMBDA_CCSD_T1 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;

   void smith_0_1(Ref<Tensor>&); //z->f1()=>out
   void smith_0_2(Ref<Tensor>&); //z->t1(),z->v2()=>out
   void offset_smith_0_3();
   void smith_0_3_0(); //z->f1()=>in.at(1)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->f1()=>in.at(2)
   void smith_2_12(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_5(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_8(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_15(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->l1(),in.at(1)=>out
   void offset_smith_0_4();
   void smith_0_4_0(); //z->f1()=>in.at(1)
   void smith_1_9(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_13();
   void smith_2_13(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_13(); //z->t1(),in.at(2)=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->l1(),in.at(1)=>out
   void offset_smith_0_6();
   void smith_1_6(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_19(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->f1(),in.at(1)=>out
   void smith_0_7(Ref<Tensor>&); //z->l1(),z->v2()=>out
   void offset_smith_0_10();
   void smith_1_10(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_34(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_10(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_11();
   void smith_1_11(); //z->t1(),z->l1()=>in.at(1)
   void smith_1_35(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_11(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_14();
   void offset_smith_1_14();
   void smith_2_14(); //z->t1(),z->l1()=>in.at(2)
   void smith_2_45(); //z->t2(),z->l2()=>in.at(2)
   void smith_1_14(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_17(); //z->t2(),z->l1()=>in.at(1)
   void offset_smith_1_44();
   void smith_2_44(); //z->t1(),z->l2()=>in.at(2)
   void smith_1_44(); //z->t2(),in.at(2)=>in.at(1)
   void smith_0_14(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_16();
   void smith_1_16(); //z->t2(),z->l1()=>in.at(1)
   void offset_smith_1_37();
   void offset_smith_2_37();
   void smith_3_37(); //z->t1(),z->l2()=>in.at(3)
   void smith_2_37(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_42(); //z->t2(),z->l2()=>in.at(2)
   void smith_1_37(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_41();
   void smith_2_41(); //z->t1(),z->l2()=>in.at(2)
   void smith_1_41(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_43();
   void smith_2_43(); //z->t2(),z->l2()=>in.at(2)
   void smith_1_43(); //z->t1(),in.at(2)=>in.at(1)
   void smith_0_16(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_18();
   void smith_0_18_0(); //z->v2()=>in.at(1)
   void offset_smith_1_18();
   void smith_1_18_0(); //z->f1()=>in.at(2)
   void smith_2_40(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_18(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_23();
   void smith_1_23_0(); //z->v2()=>in.at(2)
   void smith_2_27(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_23(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_26();
   void offset_smith_2_26();
   void smith_2_26_0(); //z->v2()=>in.at(3)
   void smith_3_36(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_26(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_38(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_26(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_30();
   void smith_1_30_0(); //z->v2()=>in.at(2)
   void smith_2_39(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_30(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_31(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_18(Ref<Tensor>&); //z->l2(),in.at(1)=>out
   void smith_0_21(Ref<Tensor>&); //z->l2(),z->v2()=>out
   void offset_smith_0_22();
   void smith_1_22(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_22(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_24();
   void smith_1_24(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_24(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_25();
   void smith_1_25(); //z->t1(),z->l2()=>in.at(1)
   void smith_0_25(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_28();
   void offset_smith_1_28();
   void smith_2_28(); //z->t1(),z->l2()=>in.at(2)
   void smith_1_28(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_32(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_28(Ref<Tensor>&); //z->v2(),in.at(1)=>out
   void offset_smith_0_29();
   void offset_smith_1_29();
   void smith_2_29(); //z->t1(),z->l2()=>in.at(2)
   void smith_1_29(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_33(); //z->t2(),z->l2()=>in.at(1)
   void smith_0_29(Ref<Tensor>&); //z->v2(),in.at(1)=>out

  public:
   LAMBDA_CCSD_T1(CCR12_Info* info);
   ~LAMBDA_CCSD_T1();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


