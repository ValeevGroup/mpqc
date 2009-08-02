//
// ccsd_r12.h --- the full CCSD-R12 method
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

#ifndef _chemistry_qc_ccr12_ccsd_r12_h
#define _chemistry_qc_ccr12_ccsd_r12_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/stateio.h>
#include <chemistry/qc/ccr12/ccr12.h>

namespace sc {

//class R12IntEval;
class R12IntEvalInfo;

class CCSD_R12: public CCR12 {

  public:
    CCSD_R12(StateIn&);
    CCSD_R12(const Ref<KeyVal>&);
    ~CCSD_R12();
    void compute();

};

}

#endif


