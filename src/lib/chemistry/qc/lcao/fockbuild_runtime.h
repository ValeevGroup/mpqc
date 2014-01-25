//
// fockbuild_runtime.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_fockbuildruntime_h
#define _mpqc_src_lib_chemistry_qc_lcao_fockbuildruntime_h

#include <util/ref/ref.h>
#include <util/group/thread.h>
#include <util/group/message.h>
#include <math/scmat/matrix.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <util/misc/registry.h>
#include <chemistry/qc/wfn/spin.h>
#include <chemistry/qc/lcao/df_runtime.h>

namespace sc {

  /// Build Fock matrices using some combination of FockBuilder objects
  class FockBuildRuntime : virtual public SavableState {
    public:
      FockBuildRuntime(const Ref<OrbitalSpaceRegistry>& oreg,
                       const Ref<AOSpaceRegistry>& aoreg,
                       const Ref<GaussianBasisSet>& refbasis,
                       const RefSymmSCMatrix& aodensity_alpha,
                       const RefSymmSCMatrix& aodensity_beta,
                       const Ref<Integral>& integral,
                       const RefSCVector& electric_field,
                       Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp(),
                       Ref<ThreadGrp> thr = ThreadGrp::get_default_threadgrp());
      FockBuildRuntime(StateIn& si);
      void save_data_state(StateOut& so);

      /// obsoletes this object
      void obsolete();

      /** Returns true if the given matrix is available
        */
      bool exists(const std::string& key) const;

      /** Returns the matrix corresponding to key.

          key must be in format recognized by ParsedOneBodyIntKey.
          If this key is not known, the matrix will be computed by an appropriate FockMatrixBuild object.
        */
      RefSCMatrix get(const std::string& key);   // non-const: can add transforms

      const Ref<Integral>& integral() const { return integral_; }
      const Ref<MessageGrp>& msg() const { return msg_; }
      const Ref<ThreadGrp>& thr() const { return thr_; }
      const Ref<GaussianBasisSet>& basis() const { return basis_; }
      const Ref<DensityFittingInfo>& dfinfo() const { return dfinfo_; }
      void dfinfo(const Ref<DensityFittingInfo>& d) { dfinfo_ = d; }
      const Ref<OrbitalSpaceRegistry>& orbital_registry() const { return oreg_; }
      /**
       * @return prec \f$ \log_2(\epsilon) \f$, where \f$ \epsilon \f$ is the abolute numerical precision of the integrals
       *         requested from the produced operator matrices
       */
      double log2_precision() const { return log2_precision_; }

      /// sets AO densities. Unless these are identical to the current densities, density-dependent contents will be cleared.
      void set_densities(const RefSymmSCMatrix& aodensity_alpha,
                         const RefSymmSCMatrix& aodensity_beta);
      /// return total density in AO basis
      RefSymmSCMatrix P() const { return P_; }
      /// return open-shell density in AO basis
      RefSymmSCMatrix Po() const { return Po_; }

      /// returns the uniform electric field (may be a null reference)
      const RefSCVector& electric_field() const { return efield_; }
      /// sets uniform electric field. In presence of electric field the core and total Fock matrices will
      /// not be cached (i.e. only field-free Fock matrices and their components are cached; assembly
      /// of in-field Fock matrices is done each time)
      void set_electric_field(const RefSCVector& efield);

      /**
       * Specifies the precision of the computed operator matrices.
       * The default precision assumed by a newly constructed FockBuildRuntime
       * is -50 ( \f$ 2^{-50} \approx 10^{-15} \f$ ).
       * Using this function to increase the precision \em may cause some
       * matrices of lower precision to be purged from cache.
       *
       * @param prec \f$ \log_2(\epsilon) \f$, where \f$ \epsilon \f$ is the absolute numerical precision of the integrals
       *        requested from the produced Fock matrices.
       */
      void set_log2_precision(double prec);

    private:

      /// obsoletes all density-dependent data
      void obsolete_density_dependents();

      // set to 1 to debug this class
      static int debug() { return 0; }

      Ref<DensityFittingInfo> dfinfo_;
      bool use_density_fitting() { return dfinfo_; }

      Ref<OrbitalSpaceRegistry> oreg_;
      Ref<AOSpaceRegistry> aoreg_;
      Ref<Integral> integral_;
      Ref<MessageGrp> msg_;
      Ref<ThreadGrp> thr_;
      Ref<GaussianBasisSet> basis_;
      RefSCVector efield_;
      bool spin_polarized_;
      double log2_precision_;

      // Densities
      RefSymmSCMatrix P_, Po_;

      // Registry of known Fock matrices
      typedef Registry<std::string, RefSCMatrix, detail::NonsingletonCreationPolicy> FockMatrixRegistry;
      Ref<FockMatrixRegistry> registry_;

      /// throws if key is not parsable by ParsedOneBodyIntKey
      void validate_key(const std::string& key) const;

      /** computes the electric_field contribution to the core hamiltonian
       *
       * @param bra_key bra key
       * @param ket_key ket key
       * @return - mu . E
       */
      RefSCMatrix electric_field_contribution(std::string bra_key,
                                              std::string ket_key);

    public:
      /// the way I compute exchange matrices is by computing square root of the density (P)
      /// this Registry keeps track of P->sqrt(P) mapping
      typedef Registry<RefSymmSCMatrix, Ref<OrbitalSpace>,
                       detail::NonsingletonCreationPolicy,
                       RefSymmSCMatrixEqual, RefObjectEqual<OrbitalSpace> > PSqrtRegistry;
    private:
      Ref<PSqrtRegistry> psqrtregistry_;

  };

  /// Parsed representation of a string key that represents a set of one-body integrals
  class ParsedOneBodyIntKey {
    public:
      ParsedOneBodyIntKey(const std::string& key);

      const std::string& key() const { return key_; }
      const std::string& bra() const { return bra_; }
      const std::string& ket() const { return ket_; }
      const std::string& oper() const { return oper_; }
      SpinCase1 spin() const { return spin_; }

      /// computes key from its components
      static std::string key(const std::string& bra,
                             const std::string& ket,
                             const std::string& oper,
                             SpinCase1 spin = AnySpinCase1);

    private:
      std::string key_;
      std::string bra_, ket_;
      std::string oper_;
      SpinCase1 spin_;
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
