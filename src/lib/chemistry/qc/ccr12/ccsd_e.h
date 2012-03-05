//
// ccsd_e.h --- energy evaluator for CCSD, CCSDT, and CCSDTQ
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

#ifndef _chemistry_qc_ccr12_ccsd_e_h
#define _chemistry_qc_ccr12_ccsd_e_h

#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCSD_E {

  protected:
   CCR12_Info* z;
 
   std::vector<Tensor*> in;

   void offset_smith_0_1();
   void smith_0_1_0(); //z->f1()=>z->in[1]
   void smith_1_2(); //z->t1(),z->v2()=>z->in[1]
   void smith_0_1(Ref<Tensor>& out); //z->t1(),z->in[1]=>out
   void smith_0_3(Ref<Tensor>& out); //z->t2(),z->v2()=>out

  public:
   CCSD_E(CCR12_Info* info);
    
   ~CCSD_E();
   void compute_amp(Ref<Tensor>& out);

};



}

#endif


