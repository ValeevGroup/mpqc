//
// ccsd_pt.h --- (T) correction to CCSD
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

#ifndef _chemistry_qc_ccr12_ccsd_pt_h
#define _chemistry_qc_ccr12_ccsd_pt_h

#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/ptnum.h>

namespace sc {

class CCSD_PT : public RefCount {

  protected:
   CCR12_Info* z;

  public:
   CCSD_PT(CCR12_Info* info) : z(info) {};
   ~CCSD_PT() {};
   
   double compute_energy(Ref<PTNum>, Ref<PTNum>);

};



}

#endif


