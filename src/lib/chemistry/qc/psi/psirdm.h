//
// psirdm.h
//
// Copyright (C) 2009 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_psi_psirdm_h
#define _mpqc_src_lib_chemistry_qc_psi_psirdm_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/mbptr12/rdm.h>

namespace sc {

  class PsiRDMCumulantTwo;

  /// PsiRDMTwo is a 2-RDM from a PsiWavefunction
  class PsiRDMTwo : public RDM<Two> {
      typedef RDMCumulant<Two> cumulant_type;
    public:
    /** A KeyVal constructor is used to generate a PsiRDMTwo
        object from the input. The full list of keywords
        that are accepted is below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>wfn</tt><td>PsiWavefunction<td>none<td>the PsiWavefunction object

        </table>
     */
      PsiRDMTwo(const Ref<KeyVal>& kv);
      PsiRDMTwo(StateIn& si);
      ~PsiRDMTwo();
      void save_data_state(StateOut& so);

      Ref<PsiWavefunction> wfn() const { return wfn_; }
      RefSymmSCMatrix scmat(SpinCase2 spincase) const;
      Ref<cumulant_type> cumulant() const;
      Ref< RDM<One> > rdm_m_1() const;
      Ref<OrbitalSpace> orbs(SpinCase1 spin) const;

    private:
      Ref<PsiWavefunction> wfn_;

      static ClassDesc class_desc_;
  };

#if 0
  /// PsiRDMCumulantTwo is the cumulant of PsiRDMTwo
  class PsiRDMCumulantTwo : public RDMCumulant<Two> {
    public:
      PsiRDMCumulantTwo(const Ref<PsiRDMTwo>& density);
      PsiRDMCumulantTwo(StateIn& si);
      ~PsiRDMCumulantTwo();
      void save_data_state(StateOut& so);

      RefSymmSCMatrix scmat(SpinCase2 spincase) const;

    private:
      Ref<PsiRDMTwo> density_;

      static ClassDesc class_desc_;
  };
#endif

  /// PsiRDMOne is a 1-RDM from a PsiWavefunction
  class PsiRDMOne : public RDM<One> {
    public:
    /** A KeyVal constructor is used to generate a PsiRDMOne
        object from the input. The full list of keywords
        that are accepted is below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>wfn</tt><td>PsiWavefunction<td>none<td>the PsiWavefunction object

        </table>
     */
      PsiRDMOne(const Ref<KeyVal>& kv);
      PsiRDMOne(StateIn& si);
      PsiRDMOne(const Ref<PsiWavefunction>& wfn);
      ~PsiRDMOne();
      void save_data_state(StateOut& so);

      Ref<OrbitalSpace> orbs(SpinCase1 spin) const;
      RefSymmSCMatrix scmat(SpinCase1 spin) const;

    private:
      Ref<PsiWavefunction> wfn_;

      static ClassDesc class_desc_;
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
