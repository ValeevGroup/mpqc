//
// ccsd_2t_pr12_right.h --- a numerator of the (2)T correction to CCSD(R12)
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
  
#ifndef _chemistry_qc_ccr12_ccsd_2t_pr12_right_h
#define _chemistry_qc_ccr12_ccsd_2t_pr12_right_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

namespace sc {

class CCSD_2T_PR12_RIGHT : public Parenthesis2tNum {

  protected:

   void offset_k0();
   void smith_k0(); //z->t1(),z->v2()=>kn.at(0)
   void offset_k1();
   void smith_k1(); //z->qy(),z->v2()=>kn.at(1)
   void offset_k2();
   void smith_k2(); //z->t1(),z->v2()=>kn.at(2)
   void offset_k3();
   void smith_k3(); //z->t1(),z->v2()=>kn.at(3)
   void offset_k4();
   void smith_k4(); //z->qy(),kn.at(3)=>kn.at(4)
   void offset_smith_0_1();
   void smith_0_1_0(); //z->v2()=>in.at(1x0)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->v2()=>in.at(2)
   void offset_smith_2_10();
   void smith_2_10_0(); //z->v2()=>in.at(3)
   void smith_2_10_1(); //kn.at(2)=>in.at(3)
   void smith_2_10(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_32(); //z->t2(),z->v2()=>in.at(2)
   void smith_2_38(); //z->c2(),z->vr2()=>in.at(2)
   void smith_1_4(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->v2()=>in.at(2)
   void smith_1_5_1(); //kn.at(0)=>in.at(2)
   void smith_1_5(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_17();
   void smith_1_17_0(); //z->v2()=>in.at(2)
   void smith_1_17_1(); //kn.at(2)=>in.at(2)
   void smith_1_17(); //z->t2(),in.at(2)=>in.at(1x0)
   void smith_1_19(); //z->t2(),z->v2()=>in.at(1x0)
   void offset_smith_1_22();
   void smith_1_22_0(); //z->v2()=>in.at(2)
   void smith_1_22_1(); //kn.at(3)=>in.at(2)
   void smith_1_22(); //z->qy(),in.at(2)=>in.at(1x0)
   void smith_1_24(); //z->c2(),z->vr2()=>in.at(1x0)
   void offset_smith_1_29();
   void smith_2_29(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_29(); //z->t2(),in.at(2)=>in.at(1x0)
   void offset_smith_1_35();
   void smith_2_35(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_35(); //z->qy(),in.at(2)=>in.at(1x0)
   void smith_0_1(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_2();
   void smith_0_2_0(); //z->v2()=>in.at(1x1)
   void smith_1_7(); //z->t1(),z->v2()=>in.at(1x1)
   void offset_smith_1_16();
   void smith_1_16_0(); //z->v2()=>in.at(2)
   void smith_1_16_1(); //kn.at(2)=>in.at(2)
   void smith_1_16(); //z->t2(),in.at(2)=>in.at(1x1)
   void smith_1_18(); //z->t2(),z->v2()=>in.at(1x1)
   void smith_1_23(); //z->qy(),z->v2()=>in.at(1x1)
   void smith_0_2(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_3();
   void smith_0_3_0(); //z->v2()=>in.at(1x2)
   void smith_1_9(); //z->t1(),z->v2()=>in.at(1x2)
   void smith_1_21(); //z->t2(),z->v2()=>in.at(1x2)
   void smith_0_3(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_6();
   void offset_smith_1_6();
   void smith_1_6_0(); //z->v2()=>in.at(2)
   void smith_1_6_1(); //kn.at(0)=>in.at(2)
   void smith_2_30(); //z->t2(),z->v2()=>in.at(2)
   void smith_2_36(); //z->qy(),z->v2()=>in.at(2)
   void smith_1_6(); //z->t2(),in.at(2)=>in.at(1x3)
   void offset_smith_1_8();
   void smith_1_8_0(); //z->v2()=>in.at(2)
   void smith_2_15(); //z->t1(),z->v2()=>in.at(2)
   void smith_2_33(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_8(); //z->qy(),in.at(2)=>in.at(1x3)
   void offset_smith_1_12();
   void smith_1_12_0(); //kn.at(1)=>in.at(2)
   void smith_1_12_1(); //kn.at(4)=>in.at(2)
   void offset_smith_2_12();
   void smith_2_12_0(); //z->v2()=>in.at(3)
   void smith_2_12_1(); //kn.at(2)=>in.at(3)
   void smith_2_12(); //z->t2(),in.at(3)=>in.at(2)
   void smith_1_12(); //z->t1(),in.at(2)=>in.at(1x3)
   void smith_0_6(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_20();
   void smith_0_20_0(); //kn.at(1)=>in.at(1x4)
   void smith_0_20_1(); //kn.at(4)=>in.at(1x4)
   void smith_0_20(double*,const long,const long,const long,const long,const long,const long);

  public:
   CCSD_2T_PR12_RIGHT(CCR12_Info* info);

   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


