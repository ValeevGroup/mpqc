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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_fockbuildruntime_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_fockbuildruntime_h

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
      const RefSCVector& electric_field() const { return efield_; }
      const Ref<DensityFittingInfo>& dfinfo() const { return dfinfo_; }
      void dfinfo(const Ref<DensityFittingInfo>& d) { dfinfo_ = d; }

      /// sets AO densities. Unless these are identical to the current densities, contents will be cleared.
      void set_densities(const RefSymmSCMatrix& aodensity_alpha,
                         const RefSymmSCMatrix& aodensity_beta);

    private:

      // set to 1 to debug this class
      static int debug() { return 0; }

      Ref<DensityFittingInfo> dfinfo_;
      bool use_density_fitting() { return dfinfo_.nonnull(); }

      Ref<OrbitalSpaceRegistry> oreg_;
      Ref<AOSpaceRegistry> aoreg_;
      Ref<Integral> integral_;
      Ref<MessageGrp> msg_;
      Ref<ThreadGrp> thr_;
      Ref<GaussianBasisSet> basis_;
      RefSCVector efield_;
      bool spin_polarized_;

      // Densities
      RefSymmSCMatrix P_, Po_;

      // Registry of known Fock matrices
      typedef Registry<std::string, RefSCMatrix, detail::NonsingletonCreationPolicy> FockMatrixRegistry;
      Ref<FockMatrixRegistry> registry_;

      /// throws if key is not parsable by ParsedOneBodyIntKey
      void validate_key(const std::string& key) const;

    public:
      /// this functor compares RefSymmSCMatrix objects. Such objects are same if every element in one
      /// differs by less than DBL_EPSILON from the corresponding element in the other.
      struct RefSymmSCMatrixEqual {
        bool operator()(const RefSymmSCMatrix& mat1, const RefSymmSCMatrix& mat2) {
          if (mat1.pointer() == mat2.pointer()) return true;
          const int n = mat1.n();
          if (n != mat2.n()) return false;
          for(int r=0; r<n; ++r)
            for(int c=0; c<=r; ++c)
              if ( fabs(mat1(r,c) - mat2(r,c)) > DBL_EPSILON)
                return false;
          return true;
        }
      };
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
