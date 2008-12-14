//
// print_scmat_norms.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <string>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <math/scmat/matrix.h>

#ifndef _chemistry_qc_mbptr12_printscmatnorms_h
#define _chemistry_qc_mbptr12_printscmatnorms_h

namespace sc {

  /// Compute and print out neatly various matrix norms of A
  template <class RefSCMat>
    void print_scmat_norms(const RefSCMat& A, const std::string& label, std::ostream& os = ExEnv::out0())
    {
      Ref<SCElementMaxAbs> maxabs_op(new SCElementMaxAbs);
      A.element_op(maxabs_op);
      const double maxabs = maxabs_op->result();

      Ref<SCElementKNorm> onenorm_op(new SCElementKNorm(1.0));
      A.element_op(onenorm_op);
      const double onenorm = onenorm_op->result();

      Ref<SCElementKNorm> twonorm_op(new SCElementKNorm(2.0));
      A.element_op(twonorm_op);
      const double twonorm = twonorm_op->result();

      os << indent << "Norms of " << label << endl;
      os << indent << "------------------------" << endl;
      os << indent << "||A||_{\\infty} = " << scprintf("%10.5lf",maxabs) << endl;
      os << indent << "||A||_1        = " << scprintf("%10.5lf",onenorm) << endl;
      os << indent << "||A||_2        = " << scprintf("%10.5lf",twonorm) << endl << endl;
    }


};

#endif

