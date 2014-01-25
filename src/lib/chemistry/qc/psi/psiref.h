//
// psiref.h
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

#ifndef _mpqc_src_lib_chemistry_qc_psi_psiref_h
#define _mpqc_src_lib_chemistry_qc_psi_psiref_h

#include <chemistry/qc/nbody/ref.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psici.h>

namespace sc {

  /// RefWavefunction specialization initialized with a PsiSCF wave function
  class PsiSCF_RefWavefunction : public RefWavefunction {
    public:
      /// construct from a PsiSCF object
      /// @param world The WavefunctionWorld in which this objects lives.
      /// @param scf The PsiSCF object that specifies the orbitals
      /// @param spin_restricted If false and obwfn is a PsiHSOSHF object,
      ///        will use semicanonical orbitals. The value of this parameter will be ignored for closed-shell
      ///        and spin-unrestricted open-shell references.
      /// @param nfzc The number of lowest-energy occupied orbitals to be kept inactive
      /// @param nfzv The number of highest-energy unoccupied orbitals to be kept inactive
      /// @param vir_space The space describing the unoccupied orbitals. Default is 0, which
      ///        means use unoccupied orbitals from obwfn.
      ///
      /// N.B. This will feed the FockBuildRuntime in world with the density matrices from obwfn!
      PsiSCF_RefWavefunction(const Ref<WavefunctionWorld>& world,
                               const Ref<PsiSCF>& scf,
                               bool spin_restricted = true,
                               unsigned int nfzc = 0,
                               unsigned int nfzv = 0,
                               Ref<OrbitalSpace> vir_space = 0);
      PsiSCF_RefWavefunction(StateIn&);
      ~PsiSCF_RefWavefunction();
      void save_data_state(StateOut&);
      const Ref<PsiSCF>& scf() const { return scf_; }
      const Ref<OrbitalSpace>& vir_space() const { return vir_space_; }
      const Ref<GaussianBasisSet>& uocc_basis() const {
        if (vir_space_) return vir_space_->basis();
        else return this->basis();
      }

      void print(std::ostream& os = ExEnv::out0()) const;

      bool sdref() const { return true; }
      double energy() { return scf()->energy(); }
      double actual_value_accuracy () const { return scf()->actual_value_accuracy(); }
      double desired_value_accuracy() const { return scf()->desired_value_accuracy(); }
      bool desired_value_accuracy_set_to_default() const { return scf()->desired_value_accuracy_set_to_default(); }
      int nelectron() const { return scf()->nelectron(); }
      double magnetic_moment() const { return scf_->magnetic_moment(); }
      bool spin_restricted() const { return spin_restricted_; }
      int dk() const { return 0; }
      Ref<GaussianBasisSet> momentum_basis() const { return 0; }
      RefSymmSCMatrix core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                 const Ref<GaussianBasisSet> &p_basis);
      unsigned int nfzc() const { return nfzc_; }
      unsigned int nfzv() const { return nfzv_; }
      RefSymmSCMatrix ordm(SpinCase1 spin) const;
      Ref<DensityFittingInfo> dfinfo() const;
    private:
      Ref<PsiSCF> scf_;
      Ref<OrbitalSpace> vir_space_;
      bool spin_restricted_;
      unsigned int nfzc_;
      unsigned int nfzv_;
      void init_spaces();
      void init_spaces_restricted();
      void init_spaces_unrestricted();
      void _set_desired_value_accuracy(double eps) { scf_->set_desired_value_accuracy(eps); }
  };

  /// RefWavefunction specialization for a general restricted-active-space multiconfiguration wave function
  class PsiRASCI_RefWavefunction : public RefWavefunction {
    public:
      /// construct from a PsiRASCI object
      /// @param world The WavefunctionWorld in which this objects lives.
      /// @param wfn The PsiRASCI object that specifies the general CI wavefunction
      /// @param spin_restricted If false and wfn is an open-shell wavefunction (wfn->spin_polarized() == true),
      ///        will use semicanonical orbitals. The value of this parameter will be ignored for closed-shell
      ///        wfn.
      /// @param nfzc The number of lowest-energy occupied orbitals to be kept inactive
      /// @param nfzv The number of highest-energy unoccupied orbitals to be kept inactive
      /// @param omit_uocc If true, omit all unoccupied orbitals (i.e. make the unoccupied space empty). N.B. This is
      ///                  not the same as "freezing" the unoccupieds.
      ///
      /// N.B. This will feed the FockBuildRuntime in world with the density matrices from wfn!
      PsiRASCI_RefWavefunction(const Ref<WavefunctionWorld>& world,
                               const Ref<PsiRASCI>& wfn,
                               bool spin_restricted = true,
                               unsigned int nfzc = 0,
                               unsigned int nfzv = 0,
                               bool omit_uocc = false);
      PsiRASCI_RefWavefunction(StateIn&);
      ~PsiRASCI_RefWavefunction();
      void save_data_state(StateOut&);
      const Ref<PsiRASCI>& wfn() const { return wfn_; }

      void print(std::ostream& os = ExEnv::out0()) const;

      bool sdref() const { return false; }
      double energy() { return wfn()->energy(); }
      double actual_value_accuracy () const { return wfn()->actual_value_accuracy(); }
      double desired_value_accuracy() const { return wfn()->desired_value_accuracy(); }
      bool desired_value_accuracy_set_to_default() const { return wfn()->desired_value_accuracy_set_to_default(); }
      int nelectron() const { return wfn()->nelectron(); }
      double magnetic_moment() const { return wfn_->magnetic_moment(); }
      bool spin_restricted() const { return spin_restricted_; }
      int dk() const { return 0; }
      Ref<GaussianBasisSet> momentum_basis() const { return 0; }
      RefSymmSCMatrix core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                 const Ref<GaussianBasisSet> &p_basis);
      unsigned int nfzc() const { return nfzc_; }
      unsigned int nfzv() const { return nfzv_; }
      bool omit_uocc() const { return omit_uocc_; }
      RefSymmSCMatrix ordm(SpinCase1 spin) const;
      Ref<DensityFittingInfo> dfinfo() const;
    private:
      Ref<PsiRASCI> wfn_;
      bool spin_restricted_;
      unsigned int nfzc_;
      unsigned int nfzv_;
      bool omit_uocc_;
      void init_spaces();
      void _set_desired_value_accuracy(double eps) { wfn_->set_desired_value_accuracy(eps); }
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
