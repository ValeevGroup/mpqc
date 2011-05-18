//
// fbclhf.h --- definition of the closed shell Hartree-Fock SCF class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_scf_fbclhf_h
#define _chemistry_qc_scf_fbclhf_h

#include <string>

#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/fockbuild.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/// CLHF is a Hartree-Fock specialization of CLSCF.
class FockBuildCLHF: public CLHF {
  protected:
    Ref<FockDistribution> fockdist_;
    Ref<FockBuild> fb_;
    std::string fockbuildmatrixtype_;
    void ao_fock(double accuracy);
    bool prefetch_blocks_;
  public:
    FockBuildCLHF(StateIn&);
    FockBuildCLHF(const Ref<KeyVal>&);
    ~FockBuildCLHF();
    void save_data_state(StateOut&);
    void init_threads();
    void done_threads();
    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End: