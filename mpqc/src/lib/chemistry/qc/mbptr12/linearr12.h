//
// linearr12.h
//
// Copyright (C) 2003 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_linearr12_h
#define _chemistry_qc_mbptr12_linearr12_h

namespace sc {
  namespace LinearR12 {
    enum StandardApproximation {StdApprox_A = 0,
				StdApprox_Ap = 1,
				StdApprox_B = 2};
    enum ABSMethod {ABS_ABS = 0,
		    ABS_ABSPlus = 1,
		    ABS_CABS = 2,
		    ABS_CABSPlus = 3};
  }

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


