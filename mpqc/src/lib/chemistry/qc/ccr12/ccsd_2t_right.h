//
// ccsd_2t_right.h -- a numerator of the (2)T correction to CCSD
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

#ifndef _chemistry_qc_ccr12_ccsd_2t_right_h
#define _chemistry_qc_ccr12_ccsd_2t_right_h

#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

namespace sc {

class CCSD_2T_RIGHT : public Parenthesis2tNum {

  protected:

   void offset_smith_0_1();
   void smith_0_1_0(); //z->v2()=>in.at(1x0)
   void offset_smith_1_1();
   void smith_1_1_0(); //z->f1()=>in.at(2)
   void smith_2_19(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_1(); //z->t2(),in.at(2)=>in.at(1x0)
   void offset_smith_1_4();
   void smith_1_4_0(); //z->v2()=>in.at(2)
   void offset_smith_2_8();
   void smith_2_8_0(); //z->v2()=>in.at(3)
   void smith_3_16(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_8(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_22(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_4(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_5();
   void smith_1_5_0(); //z->v2()=>in.at(2)
   void smith_2_9(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_5(); //z->t1(),in.at(2)=>in.at(1x0)
   void offset_smith_1_13();
   void smith_1_13_0(); //z->v2()=>in.at(2)
   void smith_2_21(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_13(); //z->t2(),in.at(2)=>in.at(1x0)
   void smith_1_15(); //z->t2(),z->v2()=>in.at(1x0)
   void smith_0_1(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_3();
   void smith_0_3_0(); //z->v2()=>in.at(1x1)
   void offset_smith_1_6();
   void smith_1_6_0(); //z->v2()=>in.at(2)
   void offset_smith_2_10();
   void smith_2_10_0(); //z->v2()=>in.at(3)
   void smith_3_17(); //z->t1(),z->v2()=>in.at(3)
   void smith_2_10(); //z->t1(),in.at(3)=>in.at(2)
   void smith_2_11(); //z->t1(),z->v2()=>in.at(2)
   void smith_2_20(); //z->t2(),z->v2()=>in.at(2)
   void smith_1_6(); //z->t1(),in.at(2)=>in.at(1x1)
   void smith_1_7(); //z->t1(),z->v2()=>in.at(1x1)
   void offset_smith_1_12();
   void smith_1_12_0(); //z->v2()=>in.at(2)
   void smith_2_18(); //z->t1(),z->v2()=>in.at(2)
   void smith_1_12(); //z->t2(),in.at(2)=>in.at(1x1)
   void smith_1_14(); //z->t2(),z->v2()=>in.at(1x1)
   void smith_0_3(double*,const long,const long,const long,const long,const long,const long);


  public:
   CCSD_2T_RIGHT(CCR12_Info* info);

   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


