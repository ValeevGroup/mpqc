//
// ccsdtq_t3.h --- computes the T3 residual vector for CCSDTQ
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

#ifndef _chemistry_qc_ccr12_ccsdtq_t3_h
#define _chemistry_qc_ccr12_ccsdtq_t3_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12.h>
#include <chemistry/qc/ccr12/ccsdtq.h>

namespace sc {

class CCSDTQ_T3 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->v2()=>kn.at(2)
   void offset_k3();
   void smith_k3(); //z->t1(),z->v2()=>kn.at(3)
   void offset_k4();
   void smith_k4(); //z->t3(),z->v2()=>kn.at(4)
   void offset_k5();
   void smith_k5(); //z->t2(),z->v2()=>kn.at(5)
   void offset_k6();
   void smith_k6(); //z->t2(),z->v2()=>kn.at(6)
   void offset_k7();
   void smith_k7(); //z->t1(),kn.at(2)=>kn.at(7)
   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>in.at(1)
   void offset_smith_1_3();
   void smith_1_3_0(); //z->f1()=>in.at(2)
   void smith_1_3_1(); //kn.at(3)=>in.at(2)
   void smith_1_3(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_20(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_47(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_22(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_49(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_4();
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>in.at(2)
   void smith_1_4_1(); //kn.at(3)=>in.at(2)
   void smith_1_4(); //z->t3(),in.at(2)=>in.at(1)
   void offset_smith_1_11();
   void smith_1_11_0(); //z->v2()=>in.at(2)
   void smith_1_11_1(); //kn.at(1)=>in.at(2)
   void smith_1_11_2(); //kn.at(5)=>in.at(2)
   void smith_1_11(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_18();
   void smith_1_18_0(); //kn.at(4)=>in.at(2)
   void offset_smith_2_18();
   void smith_2_18_0(); //z->v2()=>in.at(3)
   void smith_2_18_1(); //kn.at(2)=>in.at(3)
   void smith_2_18(); //z->t2(),in.at(3)=>in.at(2)
   void smith_1_18(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_23();
   void smith_1_23_0(); //z->v2()=>in.at(2)
   void smith_1_23_1(); //kn.at(2)=>in.at(2)
   void smith_1_23(); //z->t3(),in.at(2)=>in.at(1)
   void smith_1_25(); //z->t3(),z->v2()=>in.at(1)
   void smith_1_46(); //z->t4(),z->v2()=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>out
   void offset_smith_0_5();
   void smith_0_5_0(); //z->v2()=>in.at(1)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->f1()=>in.at(2)
   void smith_1_5_1(); //kn.at(3)=>in.at(2)
   void smith_1_5(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_9();
   void smith_1_9_0(); //z->v2()=>in.at(2)
   void smith_1_9_1(); //kn.at(0)=>in.at(2)
   void smith_1_9_2(); //kn.at(7)=>in.at(2)
   void smith_1_9_3(); //kn.at(6)=>in.at(2)
   void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_10();
   void smith_1_10_0(); //z->v2()=>in.at(2)
   void smith_1_10_1(); //kn.at(1)=>in.at(2)
   void smith_1_10(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_27();
   void smith_1_27_0(); //z->v2()=>in.at(2)
   void smith_1_27_1(); //kn.at(2)=>in.at(2)
   void smith_1_27(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_29(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_53(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_6();
   void smith_0_6_0(); //z->f1()=>in.at(1)
   void smith_0_6_1(); //kn.at(3)=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void offset_smith_0_8();
   void smith_0_8_0(); //z->v2()=>in.at(1)
   void smith_1_12(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_26();
   void smith_1_26_0(); //z->v2()=>in.at(2)
   void smith_1_26_1(); //kn.at(2)=>in.at(2)
   void smith_1_26(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_28(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_51(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>in.at(1)
   void smith_0_13_1(); //kn.at(0)=>in.at(1)
   void smith_0_13_2(); //kn.at(7)=>in.at(1)
   void smith_0_13_3(); //kn.at(6)=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_14();
   void smith_0_14_0(); //z->v2()=>in.at(1)
   void smith_0_14_1(); //kn.at(1)=>in.at(1)
   void smith_0_14_2(); //kn.at(5)=>in.at(1)
   void smith_0_14(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void smith_0_15(Ref<Tensor>&); //z->t3(),z->v2()=>out
   void offset_smith_0_30();
   void smith_0_30_0(); //z->v2()=>in.at(1)
   void smith_0_30_1(); //kn.at(2)=>in.at(1)
   void smith_0_30(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void smith_0_31(Ref<Tensor>&); //z->t4(),z->v2()=>out
   void smith_0_52(Ref<Tensor>&); //z->t2(),kn.at(4)=>out

  public:
   CCSDTQ_T3(CCR12_Info* info);
    
   ~CCSDTQ_T3();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


