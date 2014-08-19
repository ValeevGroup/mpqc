//
// gto.h
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

#ifndef _mpqc_src_lib_chemistry_qc_libint2_gto_h
#define _mpqc_src_lib_chemistry_qc_libint2_gto_h

#include <libint2.h>
#include <libint2/config.h>
#include <vector>

namespace sc {

  /// Provides precomputed information about Gaussian basis functions
  class GTOInfo : public RefCount {
    public:
    static const Ref<GTOInfo>& instance();
    ~GTOInfo();

    const double* norm(unsigned int l);
    const double* fp1();

    private:
      static Ref<GTOInfo> instance_;
#ifdef LIBINT_MAX_AM
      static const unsigned int lmax_ = LIBINT_MAX_AM;
#else
#  ifdef LIBINT2_MAX_AM_ERI
      static const unsigned int lmax_ = LIBINT2_MAX_AM_ERI;
#  else
      static const unsigned int lmax_ = LIBINT2_CARTGAUSS_MAX_AM;
#  endif
#endif

      GTOInfo();

      std::vector< std::vector<double> > norm_;
      std::vector<double> df_; // double factorials
      std::vector<double> fp1_; // vector of 1.0s
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
