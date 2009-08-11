//
// ccsdpr12_c.h --- computes the GT2 residual vector of CCSD(R12)
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

#ifndef _chemistry_qc_ccr12_ccsdpr12_c_h
#define _chemistry_qc_ccr12_ccsdpr12_c_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSDPR12_C {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_smith_0_1();
   void smith_0_1_0(); //z->vd2()=>in.at(1)
   void offset_smith_1_1();
   void smith_1_1_0(); //z->f1()=>in.at(2)
   void smith_2_12(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_1(); //z->fy(),in.at(2)=>in.at(1)
   void offset_smith_1_14();
   void smith_2_14(); //z->t1(),z->fy()=>in.at(2)
   void smith_1_14(); //z->v2(),in.at(2)=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_1_2(); //z->c2(),z->f1()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->xs2(),in.at(1)=>Ref<Tensor>&
   void smith_0_3(Ref<Tensor>&); //z->c2(),z->bs2()=>Ref<Tensor>&
   void smith_0_4(Ref<Tensor>&); //z->vd2()=>Ref<Tensor>&
   void offset_smith_0_5();
   void smith_1_5(); //z->t1(),z->fy()=>in.at(1)
   void smith_0_5(Ref<Tensor>&); //z->v2(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_0_6_0(); //z->vd2()=>in.at(1)
   void offset_smith_1_7();
   void smith_2_7(); //z->t1(),z->fy()=>in.at(2)
   void smith_1_7(); //z->v2(),in.at(2)=>in.at(1)
   void smith_1_8(); //z->t1(),z->vd2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_9();
   void offset_smith_1_9();
   void smith_1_9_0(); //z->v2()=>in.at(2)
   void smith_2_13(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_9(); //z->t2(),in.at(2)=>in.at(1)
   void offset_smith_1_11();
   void offset_smith_2_11();
   void smith_3_11(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_11(); //z->t1(),in.at(3)=>in.at(2)
   void smith_1_11(); //z->t1(),in.at(2)=>in.at(1)
   void smith_0_9(Ref<Tensor>&); //z->fy(),in.at(1)=>Ref<Tensor>&

  public:
   CCSDPR12_C(CCR12_Info* info);
   ~CCSDPR12_C();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


