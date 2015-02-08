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
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/lcao/df.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/// FockBuildCLHF is a specialization of CLHF that uses FockBuild class for computing fock matrices
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

/// DFCLHF is a specialization of CLHF that uses a density-fitting FockBuild class for computing fock matrices
class DFCLHF: public CLHF {
  protected:
    void ao_fock(double accuracy);
    void reset_density();
  public:
    DFCLHF(StateIn&);
    /** Accepts all keywords of CLHF class + the following keywords:
        <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          <tr><td><tt>world</tt><td>WavefunctionWorld<td>see notes<td>the WavefunctionWorld object that
          this wave function belongs to. If not given, this object will live in its own WavefunctionWorld.
          It is recommended then that correlated Wavefunction objects that use this as a reference provide their
          own worlds to this.

          </table>

     * N.B. <tt>density_reset_freq</tt> is ignored by this method -- full Fock matrix is always constructed.
     */
    DFCLHF(const Ref<KeyVal>&);
    ~DFCLHF();
    void save_data_state(StateOut&);
    void print(std::ostream&o=ExEnv::out0()) const;
    const Ref<WavefunctionWorld>& world() const { return world_; }
    Ref<DensityFittingInfo> dfinfo() const;
#ifdef MPQC_NEW_FEATURES
    virtual boost::property_tree::ptree& write_xml(
        boost::property_tree::ptree& parent, const XMLWriter& writer
    );
#endif // MPQC_NEW_FEATURES
  private:
    RefSymmSCMatrix gmat_;
    Ref<WavefunctionWorld> world_;
    static ClassDesc cd_;
#ifdef MPQC_NEW_FEATURES
    bool xml_debug_ = false;
#endif
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
