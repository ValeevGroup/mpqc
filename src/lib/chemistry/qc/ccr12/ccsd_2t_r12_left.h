//
// ccsd_2t_r12_left.h --- a numerator of the (2)T correction to CCSD(R12) or CCSD-R12
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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
  
#ifndef _chemistry_qc_ccr12_ccsd_2t_r12_left_h
#define _chemistry_qc_ccr12_ccsd_2t_r12_left_h

#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

namespace sc {

class CCSD_2T_R12_LEFT : public Parenthesis2tNum {

  protected:

   void smith_0_1(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_2();
   void smith_0_2_0(); //z->f1()=>in.at(1x0)
   void smith_1_6(); //z->t1(),z->v2()=>in.at(1x0)
   void smith_0_2(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_3();
   void smith_0_3_0(); //z->v2()=>in.at(1x1)
   void smith_1_7(); //z->t1(),z->v2()=>in.at(1x1)
   void smith_0_3(double*,const long,const long,const long,const long,const long,const long);
   void smith_0_4(double*,const long,const long,const long,const long,const long,const long);
   void smith_0_5(double*,const long,const long,const long,const long,const long,const long);
   void offset_smith_0_8();
   void smith_1_8(); //z->t1(),z->l2()=>in.at(1x2)
   void smith_0_8(double*,const long,const long,const long,const long,const long,const long);

  public:
   CCSD_2T_R12_LEFT(CCR12_Info* info);

   void compute_amp(double*,const long,const long,const long,const long,const long,const long,const long);

};



}

#endif


