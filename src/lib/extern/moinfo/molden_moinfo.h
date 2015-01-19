//
// moinfo.h
//
// Copyright (C) 2011 Edward Valeev
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

#ifndef _mpqc_src_lib_extern_moinfo_molden_moinfo_h
#define _mpqc_src_lib_extern_moinfo_molden_moinfo_h

#include <vector>
#include <string>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/wfn/rdm.h>
#include <math/scmat/abstract.h>

namespace sc {

  /// Reads MO information from a text MOLDEN file
  class MOLDEN_ExternReadMOInfo {
    public:
      MOLDEN_ExternReadMOInfo(const std::string& filename);
      ~MOLDEN_ExternReadMOInfo() {}
      Ref<GaussianBasisSet> basis() const;
      RefSCMatrix coefs() const;
      std::vector<unsigned int> orbsym() const;

    private:
      Ref<GaussianBasisSet> basis_;
      RefSCMatrix coefs_;
      std::vector<unsigned int> orbsym_;
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
