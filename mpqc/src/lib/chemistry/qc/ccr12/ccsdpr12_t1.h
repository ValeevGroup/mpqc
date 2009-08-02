//
// ccsdpr12_t1.h --- computes the T1 residual vector of CCSD(R12)
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

#ifndef _chemistry_qc_ccr12_ccsdpr12_t1_h
#define _chemistry_qc_ccr12_ccsdpr12_t1_h

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


namespace sc {

class CCSDPR12_T1 {

  protected:
   CCR12_Info* z;
   std::vector<Tensor*> in;
   std::vector<Tensor*> kn;

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>in.at(1)
   void smith_1_4(); //z->t1(),z->v2()=>in.at(1)
   void smith_1_10(); //z->t1(),kn.at(0)=>in.at(1)
   void smith_1_13(); //z->t2(),z->v2()=>in.at(1)
   void smith_1_16(); //z->c2(),z->vr2()=>in.at(1)
   void smith_0_1(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1)
   void smith_1_5(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_2(Ref<Tensor>&); //z->t1(),in.at(1)=>Ref<Tensor>&
   void smith_0_3(Ref<Tensor>&); //z->t1(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_6();
   void smith_0_6_0(); //z->v2()=>in.at(1)
   void smith_1_12(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_6(Ref<Tensor>&); //z->t2(),in.at(1)=>Ref<Tensor>&
   void smith_0_7(Ref<Tensor>&); //z->t2(),z->v2()=>Ref<Tensor>&
   void offset_smith_0_8();
   void smith_0_8_0(); //z->v2()=>in.at(1)
   void smith_1_15(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_8(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&
   void smith_0_9(Ref<Tensor>&); //z->c2(),z->vr2()=>Ref<Tensor>&
   void smith_0_11(Ref<Tensor>&); //z->t2(),kn.at(0)=>Ref<Tensor>&
   void offset_smith_0_14();
   void smith_1_14(); //z->t1(),z->v2()=>in.at(1)
   void smith_0_14(Ref<Tensor>&); //z->qy(),in.at(1)=>Ref<Tensor>&

  public:
   CCSDPR12_T1(CCR12_Info* info);
   ~CCSDPR12_T1();
    
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


