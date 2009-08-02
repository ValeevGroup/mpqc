//
// ccsdt.h -- the coupled-cluster singles, doubles, and triples 
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

#ifndef _chemistry_qc_ccr12_ccsdt_h
#define _chemistry_qc_ccr12_ccsdt_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/stateio.h>
#include <chemistry/qc/ccr12/ccr12.h>

namespace sc {

class CCSDT: public CCR12 {

  public:
    CCSDT(StateIn&);
    CCSDT(const Ref<KeyVal>&);
    ~CCSDT();
    void compute();

};

}

#endif


