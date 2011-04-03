//
// ccsdtq_t2.h --- computes the T2 residual vector for CCSDTQ
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

#ifndef _chemistry_qc_ccr12_ccsdtq_t2_h
#define _chemistry_qc_ccr12_ccsdtq_t2_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSDTQ_T2 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;

   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>in.at(1)
   void offset_smith_1_3();
   void smith_1_3_0(); //z->f1()=>in.at(2)
   void smith_2_26(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_3(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_17(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_37(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_19(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_35(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_4();
   void smith_0_4_0(); //z->v2()=>in.at(1)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>in.at(2)
   void smith_2_28(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_4(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_9();
   void smith_1_9_0(); //z->v2()=>in.at(2)
   void offset_smith_2_15();
   void smith_2_15_0(); //z->v2()=>in.at(3)
   void smith_3_25(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_15(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_30(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_10();
   void smith_1_10_0(); //z->v2()=>in.at(2)
   void smith_2_16(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_10(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_20();
   void smith_1_20_0(); //z->v2()=>in.at(2)
   void smith_2_29(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_20(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_22(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_33(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>out
   void offset_smith_0_5();
   void smith_0_5_0(); //z->f1()=>in.at(1)
   void smith_1_31(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void smith_0_6(Ref<Tensor>&); //z->v2()=>out
   void offset_smith_0_8();
   void smith_0_8_0(); //z->v2()=>in.at(1)
   void smith_1_11(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->t1(),in.at(1)=>out
   void offset_smith_0_12();
   void smith_0_12_0(); //z->v2()=>in.at(1)
   void offset_smith_1_18();
   void smith_1_18_0(); //z->v2()=>in.at(2)
   void smith_2_27(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_18(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_36(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_12(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_13();
   void smith_0_13_0(); //z->v2()=>in.at(1)
   void smith_1_21(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_34(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_13(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void smith_0_14(Ref<Tensor>&); //z->t2(),z->v2()=>out
   void offset_smith_0_23();
   void smith_0_23_0(); //z->v2()=>in.at(1)
   void smith_1_32(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_23(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void smith_0_24(Ref<Tensor>&); //z->t3(),z->v2()=>out
   void smith_0_38(Ref<Tensor>&); //z->t4(),z->v2()=>out

  public:
   CCSDTQ_T2(CCR12_Info* info);
    
   ~CCSDTQ_T2();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


