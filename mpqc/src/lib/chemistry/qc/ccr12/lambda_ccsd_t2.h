//
// lambda_ccsd_t2.h --- computes the L2 residual vector for Lambda-CCSD
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

#ifndef _chemistry_qc_ccr12_lambda_ccsd_t2_h
#define _chemistry_qc_ccr12_lambda_ccsd_t2_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class LAMBDA_CCSD_T2 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>z->kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>z->kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->l2()=>z->kn.at(2)
   void smith_0_1(Ref<Tensor>&); //z->v2()=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>z->in.at(1)
   void smith_0_2_1(); //z->kn.at(0)=>z->in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->l1(),z->in.at(1)=>out
   void offset_smith_0_3();
   void smith_0_3_0(); //z->v2()=>z->in.at(1)
   void smith_0_3_1(); //z->kn.at(1)=>z->in.at(1)
   void smith_0_3(Ref<Tensor>&); //z->l1(),z->in.at(1)=>out
   void smith_0_4(Ref<Tensor>&); //z->l1(),z->v2()=>out
   void offset_smith_0_7();
   void smith_1_7(); //z->t1(),z->l1()=>z->in.at(1)
   void smith_1_32(); //z->t2(),z->l2()=>z->in.at(1)
   void smith_0_7(Ref<Tensor>&); //z->v2(),z->in.at(1)=>out
   void offset_smith_0_8();
   void smith_0_8_0(); //z->f1()=>z->in.at(1)
   void offset_smith_1_10();
   void smith_1_10_0(); //z->f1()=>z->in.at(2)
   void smith_1_10_1(); //z->kn.at(0)=>z->in.at(2)
   void smith_1_10(); //z->t1(),z->in.at(2)=>z->in.at(1)
   void smith_1_15(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_1_26(); //z->t2(),z->v2()=>z->in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->l2(),z->in.at(1)=>out
   void offset_smith_0_9();
   void smith_0_9_0(); //z->f1()=>z->in.at(1)
   void smith_1_17(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_1_23(); //z->t1(),z->kn.at(0)=>z->in.at(1)
   void smith_1_28(); //z->t2(),z->v2()=>z->in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->l2(),z->in.at(1)=>out
   void smith_0_11(Ref<Tensor>&); //z->f1(),z->kn.at(2)=>out
   void offset_smith_0_12();
   void smith_0_12_0(); //z->v2()=>z->in.at(1)
   void offset_smith_1_16();
   void smith_1_16_0(); //z->v2()=>z->in.at(2)
   void smith_1_16_1(); //z->kn.at(1)=>z->in.at(2)
   void smith_1_16(); //z->t1(),z->in.at(2)=>z->in.at(1)
   void smith_1_27(); //z->t2(),z->v2()=>z->in.at(1)
   void smith_0_12(Ref<Tensor>&); //z->l2(),z->in.at(1)=>out
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>z->in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>z->in.at(1)
   void smith_1_24(); //z->t1(),z->kn.at(1)=>z->in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->l2(),z->in.at(1)=>out
   void smith_0_14(Ref<Tensor>&); //z->l2(),z->v2()=>out
   void smith_0_18(Ref<Tensor>&); //z->v2(),z->kn.at(2)=>out
   void smith_0_20(Ref<Tensor>&); //z->v2(),z->kn.at(2)=>out
   void offset_smith_0_25();
   void smith_1_25(); //z->t1(),z->kn.at(2)=>z->in.at(1)
   void smith_1_31(); //z->t2(),z->l2()=>z->in.at(1)
   void smith_0_25(Ref<Tensor>&); //z->v2(),z->in.at(1)=>out
   void offset_smith_0_29();
   void smith_1_29(); //z->t2(),z->l2()=>z->in.at(1)
   void smith_0_29(Ref<Tensor>&); //z->v2(),z->in.at(1)=>out
   void offset_smith_0_30();
   void smith_1_30(); //z->t2(),z->l2()=>z->in.at(1)
   void smith_0_30(Ref<Tensor>&); //z->v2(),z->in.at(1)=>out

  public:
   LAMBDA_CCSD_T2(CCR12_Info* info);
   ~LAMBDA_CCSD_T2();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


