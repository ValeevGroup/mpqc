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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_bin_pt2r12_moinfo_h
#define _mpqc_src_bin_pt2r12_moinfo_h

#include <vector>
#include <string>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/wfn/rdm.h>
#include <math/scmat/abstract.h>

namespace sc {

  /// Reads MO information from a text file
  class ExternReadMOInfo {
    public:
      ExternReadMOInfo(const std::string& filename);
      ~ExternReadMOInfo() {}
      Ref<GaussianBasisSet> basis() const;
      RefSCMatrix coefs() const;
      unsigned int nfzc() const;
      unsigned int nfzv() const;
      unsigned int nocc() const;
      std::vector<unsigned int> orbsym() const;

    private:
      Ref<GaussianBasisSet> basis_;
      RefSCMatrix coefs_;
      std::vector<unsigned int> orbsym_;
      unsigned int nocc_;
      unsigned int nfzc_;
      unsigned int nfzv_;
  };

  /// Reads 1-RDM from a text file
  class ExternReadRDMOne : public RDM<One>{
    public:
      /// reads 1-rdm from filename
      /// assumes that 1-rdm is expressed in orbs
      ExternReadRDMOne(const std::string& filename,
                       const Ref<OrbitalSpace>& orbs);
      /// receives 1-rdm as a constructor argument
      /// assumes that 1-rdm is expressed in orbs
      ExternReadRDMOne(const RefSymmSCMatrix& rdm,
                       const Ref<OrbitalSpace>& orbs);
      virtual ~ExternReadRDMOne();

      /// cannot be obsoleted
      void obsolete() {}
      /// already computed
      void compute() {}
      /// the orbital space of spincase s in which the density is reported
      Ref<OrbitalSpace> orbs(SpinCase1 s) const { return orbs_; }
      /// density matrix
      RefSymmSCMatrix scmat(SpinCase1 spin) const { return rdm_; }

    private:
      static ClassDesc class_desc_;

    private:
      Ref<OrbitalSpace> orbs_;
      // spin-independent
      RefSymmSCMatrix rdm_;
  };

  /// Reads 2-RDM from a text file
  class ExternReadRDMTwo : public RDM<Two>{
    public:
      typedef RDMCumulant<Two> cumulant_type;

      /// reads 2-rdm from filename
      /// assumes that 2-rdm is expressed in orbs
      ExternReadRDMTwo(const std::string& filename,
                       const Ref<OrbitalSpace>& orbs);
      virtual ~ExternReadRDMTwo();

      /// cannot be obsoleted
      void obsolete() {}
      /// already computed
      void compute() {}
      /// the orbital space of spincase s in which the density is reported
      Ref<OrbitalSpace> orbs(SpinCase1 s) const { return orbs_; }
      /// density matrix
      RefSymmSCMatrix scmat(SpinCase2 spin) const { return rdm_; }
      Ref<cumulant_type> cumulant() const;
      Ref< RDM<One> > rdm_m_1() const;

    private:
      static ClassDesc class_desc_;

    private:
      Ref<OrbitalSpace> orbs_;
      // spin-independent
      RefSymmSCMatrix rdm_;
      std::string filename_; // filename from which this was constructed -- may be useful to find rdm1 file
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
