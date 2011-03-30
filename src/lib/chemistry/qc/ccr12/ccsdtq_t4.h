//
// ccsdtq_t4.h --- computes the T4 residual vector for CCSDTQ
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
  
#ifndef _chemistry_qc_ccr12_ccsdtq_t4_h
#define _chemistry_qc_ccr12_ccsdtq_t4_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSDTQ_T4 {

  protected:
   CCR12_Info* z;

   std::vector<Tensor*> in;

   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>in.at(1)
   void offset_smith_1_3();
   void smith_1_3_0(); //z->f1()=>in.at(2)
   void smith_2_43(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_3(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_23(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_64(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_25(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_66(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void offset_smith_0_4();
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>in.at(2)
   void smith_2_45(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_4(); //z->t4(),in.at(2)=>in.at(1)
   void offset_smith_1_11();
   void smith_1_11_0(); //z->v2()=>in.at(2)
   void smith_2_22(); //z->t1(),z->v2()=>in.at(2)
   void smith_2_54(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_11(); //z->t3(),in.at(2)=>in.at(1)
   void offset_smith_1_21();
   void offset_smith_2_21();
   void smith_2_21_0(); //z->v2()=>in.at(3)
   void smith_3_42(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_21(); //z->t3(),in.at(3)=>in.at(2)
   void smith_2_47(); //z->t4(),z->v2()=>in.at(2)
   void offset_smith_2_58();
   void smith_3_58(); //z->t2(),z->v2()=>in.at(3)
   void smith_2_58(); //z->t2(),in.at(3)=>in.at(2)
   void smith_1_21(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_26();
   void smith_1_26_0(); //z->v2()=>in.at(2)
   void smith_2_46(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_26(); //z->t4(),in.at(2)=>in.at(1)
   void smith_1_28(); //z->t4(),z->v2()=>in.at(1)
   void offset_smith_1_30();
   void smith_2_30(); //z->t2(),z->v2()=>in.at(2)
   void smith_2_56(); //z->t3(),z->v2()=>in.at(2)
   void smith_1_30(); //z->t2(),in.at(2)=>in.at(1)
   void smith_0_4(Ref<Tensor>&); //z->t1(),in.at(1)=>out
   void offset_smith_0_5();
   void smith_0_5_0(); //z->v2()=>in.at(1)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->f1()=>in.at(2)
   void smith_2_48(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_5(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_9();
   void smith_1_9_0(); //z->v2()=>in.at(2)
   void offset_smith_2_19();
   void smith_2_19_0(); //z->v2()=>in.at(3)
   void smith_3_41(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_19(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_50(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_9(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_10();
   void smith_1_10_0(); //z->v2()=>in.at(2)
   void smith_2_20(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_10(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_33();
   void smith_1_33_0(); //z->v2()=>in.at(2)
   void smith_2_49(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_33(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_34(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_74(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_6();
   void offset_smith_1_6();
   void smith_1_6_0(); //z->f1()=>in.at(2)
   void smith_2_52(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_6(); //z->t3(),in.at(2)=>in.at(1)
   void offset_smith_1_13();
   void smith_1_13_0(); //z->v2()=>in.at(2)
   void offset_smith_2_31();
   void smith_2_31_0(); //z->v2()=>in.at(3)
   void smith_3_59(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_31(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_63(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_13(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_14();
   void smith_1_14_0(); //z->v2()=>in.at(2)
   void smith_2_32(); //z->t1(),z->v2()=>in.at(2)
   void smith_2_62(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_14(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_29();
   void offset_smith_2_29();
   void smith_2_29_0(); //z->v2()=>in.at(3)
   void smith_3_60(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_29(); //z->t2(),in.at(3)=>in.at(2)
   void smith_2_57(); //z->t3(),z->v2()=>in.at(2)
   void smith_1_29(); //z->t1(),in.at(2)=>in.at(1)
   void offset_smith_1_38();
   void smith_1_38_0(); //z->v2()=>in.at(2)
   void smith_2_55(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_38(); //z->t3(),in.at(2)=>in.at(1)
   void smith_1_40(); //z->t3(),z->v2()=>in.at(1)
   void smith_1_70(); //z->t4(),z->v2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_8();
   void smith_0_8_0(); //z->v2()=>in.at(1)
   void smith_1_12(); //z->t1(),z->v2()=>in.at(1)
   void offset_smith_1_36();
   void smith_1_36_0(); //z->v2()=>in.at(2)
   void smith_2_53(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_36(); //z->t2(),in.at(2)=>in.at(1)
   void smith_1_37(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_72(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_15();
   void offset_smith_1_15();
   void smith_1_15_0(); //z->v2()=>in.at(2)
   void smith_2_61(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_15(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_35();
   void smith_1_35_0(); //z->v2()=>in.at(2)
   void smith_2_51(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_35(); //z->t3(),in.at(2)=>in.at(1)
   void smith_1_39(); //z->t3(),z->v2()=>in.at(1)
   void smith_1_68(); //z->t4(),z->v2()=>in.at(1)
   void smith_0_15(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_16();
   void smith_0_16_0(); //z->v2()=>in.at(1)
   void offset_smith_1_24();
   void smith_1_24_0(); //z->v2()=>in.at(2)
   void smith_2_44(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_24(); //z->t1(),in.at(2)=>in.at(1)
   void smith_1_65(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_16(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void offset_smith_0_17();
   void smith_0_17_0(); //z->v2()=>in.at(1)
   void smith_1_27(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_67(); //z->t2(),z->v2()=>in.at(1)
   void smith_0_17(Ref<Tensor>&); //z->t4(),in.at(1)=>out
   void smith_0_18(Ref<Tensor>&); //z->t4(),z->v2()=>out
   void offset_smith_0_69();
   void smith_1_69(); //z->t4(),z->v2()=>in.at(1)
   void smith_0_69(Ref<Tensor>&); //z->t2(),in.at(1)=>out
   void offset_smith_0_71();
   void smith_1_71(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_71(Ref<Tensor>&); //z->t3(),in.at(1)=>out
   void offset_smith_0_73();
   void smith_1_73(); //z->t3(),z->v2()=>in.at(1)
   void smith_0_73(Ref<Tensor>&); //z->t3(),in.at(1)=>out

  public:
   CCSDTQ_T4(CCR12_Info* info);
    
   ~CCSDTQ_T4();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


