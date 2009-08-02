//
// ccsd_t2.h --- computes the T2 residual vector of CCSD
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

#ifndef _chemistry_qc_ccr12_ccsd_t2_h
#define _chemistry_qc_ccr12_ccsd_t2_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_T2 {

  protected:
   CCR12_Info* z;
 
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>z->kn[0]
   void offset_k1();
   void smith_k1(); //z->t1(),z->v2()=>z->kn[1]
   void offset_k2();
   void smith_k2(); //z->t1(),z->v2()=>z->kn[2]
   void offset_k3();
   void smith_k3(); //z->t1(),z->v2()=>z->kn[3]
   void offset_k4();
   void smith_k4(); //z->t2(),z->v2()=>z->kn[4]
   void offset_k5();
   void smith_k5(); //z->t1(),z->kn[2]=>z->kn[5]
   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>z->in[1]
   void offset_smith_1_3();
   void smith_1_3_0(); //z->f1()=>z->in[2]
   void smith_1_3_1(); //z->kn[3]=>z->in[2]
   void smith_1_3(); //z->t1(),z->in[2]=>z->in[1]
   void smith_1_16(); //z->t1(),z->v2()=>z->in[1]
   void smith_1_31(); //z->t2(),z->v2()=>z->in[1]
   void smith_0_1(Ref<Tensor>& out); //z->t2(),z->in[1]=>out
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>z->in[1]
   void smith_1_18(); //z->t1(),z->v2()=>z->in[1]
   void smith_1_29(); //z->t2(),z->v2()=>z->in[1]
   void smith_0_2(Ref<Tensor>& out); //z->t2(),z->in[1]=>out
   void offset_smith_0_4();
   void smith_0_4_0(); //z->v2()=>z->in[1]
   void offset_smith_1_4();
   void smith_1_4_0(); //z->f1()=>z->in[2]
   void smith_1_4_1(); //z->kn[3]=>z->in[2]
   void smith_1_4(); //z->t2(),z->in[2]=>z->in[1]
   void offset_smith_1_8();
   void smith_1_8_0(); //z->v2()=>z->in[2]
   void smith_1_8_1(); //z->kn[0]=>z->in[2]
   void smith_1_8_2(); //z->kn[5]=>z->in[2]
   void smith_1_8_3(); //z->kn[4]=>z->in[2]
   void smith_1_8(); //z->t1(),z->in[2]=>z->in[1]
   void offset_smith_1_9();
   void smith_1_9_0(); //z->v2()=>z->in[2]
   void smith_1_9_1(); //z->kn[1]=>z->in[2]
   void smith_1_9(); //z->t1(),z->in[2]=>z->in[1]
   void offset_smith_1_19();
   void smith_1_19_0(); //z->v2()=>z->in[2]
   void smith_1_19_1(); //z->kn[2]=>z->in[2]
   void smith_1_19(); //z->t2(),z->in[2]=>z->in[1]
   void smith_1_21(); //z->t2(),z->v2()=>z->in[1]
   void smith_0_4(Ref<Tensor>& out); //z->t1(),z->in[1]=>out
   void smith_0_5(Ref<Tensor>& out); //z->v2()=>out
   void offset_smith_0_7();
   void smith_0_7_0(); //z->v2()=>z->in[1]
   void smith_1_10(); //z->t1(),z->v2()=>z->in[1]
   void smith_0_7(Ref<Tensor>& out); //z->t1(),z->in[1]=>out
   void offset_smith_0_11();
   void smith_0_11_0(); //z->v2()=>z->in[1]
   void smith_0_11_1(); //z->kn[0]=>z->in[1]
   void smith_0_11_2(); //z->kn[5]=>z->in[1]
   void smith_0_11_3(); //z->kn[4]=>z->in[1]
   void smith_0_11(Ref<Tensor>& out); //z->t2(),z->in[1]=>out
   void offset_smith_0_12();
   void smith_0_12_0(); //z->v2()=>z->in[1]
   void smith_0_12_1(); //z->kn[1]=>z->in[1]
   void smith_1_28(); //z->t2(),z->v2()=>z->in[1]
   void smith_0_12(Ref<Tensor>& out); //z->t2(),z->in[1]=>out
   void smith_0_13(Ref<Tensor>& out); //z->t2(),z->v2()=>out

  public:
   CCSD_T2(CCR12_Info* info);
    
   ~CCSD_T2();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


