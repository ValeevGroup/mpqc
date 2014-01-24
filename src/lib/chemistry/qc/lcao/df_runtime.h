//
// df_runtime.h
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

#ifndef _mpqc_src_lib_chemistry_qc_df_df_runtime_h
#define _mpqc_src_lib_chemistry_qc_df_df_runtime_h

#include <chemistry/qc/lcao/df.h>
#include <chemistry/qc/lcao/tbint_runtime.h>

namespace sc {

  /** Parsed representation of a string key that represents fitting of a product of space1 and space2 into fspace.
      kernel must be parsable by ParsedTwoBodyOperKey. */
  class ParsedDensityFittingKey {
    public:
      ParsedDensityFittingKey(const std::string& key);

      const std::string& key() const { return key_; }
      const std::string& space1() const { return space1_; }
      const std::string& space2() const { return space2_; }
      const std::string& fspace() const { return fspace_; }
      const std::string& kernel() const { return kernel_pkey_.key(); }
      const std::string& kernel_oper() const { return kernel_pkey_.oper(); }
      const std::string& kernel_param() const { return kernel_pkey_.params(); }

      /// computes key from its components
      static std::string key(const std::string& space1,
                             const std::string& space2,
                             const std::string& fspace,
                             const std::string& kernel);

    private:
      std::string key_;
      std::string space1_, space2_, fspace_;
      ParsedTwoBodyOperKey kernel_pkey_;
  };

  class DensityFittingParams;

  /**
   *    Smart runtime support for managing DensityFitting objects
   */
  class DensityFittingRuntime : virtual public SavableState {
    public:
      typedef DensityFittingRuntime this_type;
      typedef DistArray4 Result;
      typedef Ref<Result> ResultRef;
      typedef ParsedDensityFittingKey ParsedResultKey;
      typedef DensityFitting::MOIntsRuntime MOIntsRuntime;

      // uses MOIntsRuntime to evaluate integrals
      DensityFittingRuntime(const Ref<MOIntsRuntime>& moints_runtime,
                            const DensityFittingParams* dfparams);
      DensityFittingRuntime(StateIn& si);
      void save_data_state(StateOut& so);

      /// obsoletes this object
      void obsolete();

      /// which density fitting solver will be used to compute?
      const DensityFittingParams* dfparams() const { return dfparams_; }

      /** Returns true if the given DensityFitting is available
        */
      bool exists(const std::string& key) const;

      /** Returns the DistArray4 object corresponding to this key.

          key must be in format recognized by ParsedDensityFittingKey.
          If this key is not known, the DistArray4 object will be created
          and (possibly) computed.
        */
      ResultRef get(const std::string& key);   // non-const: can compute something

      /// returns the runtime used to compute results
      const Ref<MOIntsRuntime>& moints_runtime() const { return moints_runtime_; }

      /// removes all entries that contain this space
      void remove_if(const std::string& space_key);

      /**
       * tries to translate a library basis set label to the corresponding default value for the DF basis
       * @param obs_name orbital basis set name; to be useful must be a canonical library name
       * @param incX optional parameter to raise the cardinal number of the density fitting basis. This may be useful if using density-fitting
       *             slightly outside their intended area of application.
       * @param force_aug optional parameter similar to incX to force inclusion of diffuse components of the density fitting basis
       * @return canonical library name of the default density-fitting basis set; if not able to suggest the default basis, returns an empty string
       */
      static std::string default_dfbs_name(const std::string& obs_name,
                                           int incX = 0,
                                           bool force_aug = false);

    private:
      Ref<MOIntsRuntime> moints_runtime_;
      const DensityFittingParams* dfparams_;

      typedef Registry<std::string, ResultRef, detail::NonsingletonCreationPolicy > ResultRegistry;
      Ref<ResultRegistry> results_;

      // creates the result for a given key
      const ResultRef& create_result(const std::string& key);

      static ClassDesc class_desc_;

  };

  /** DensityFittingParams defines parameters used by DensityFittingRuntime and other runtime components
      to compute density fitting objects.
    */
  class DensityFittingParams : virtual public SavableState {
    public:
    /**
     * @param basis  The GaussianBasisSet object used to fit product densities. There is no default.
     *               @note DensityFittingRuntime does not use this, but other runtime objects may use it
     *               to set the global density fitting basis.
     * @param kernel_key A string describing the kernel_key. It must be parsable by ParsedTwoBodyOperKey, or be empty (the default).
     *               @note DensityFittingRuntime does not use this, but other runtime objects may use it to set the global density fitting method.
     * @param solver A string describing the method of solving the density fitting equations. This is used by DensityFittingRuntime
     *               to produce density fitting objects.
     */
      DensityFittingParams(const Ref<GaussianBasisSet>& basis,
                           const std::string& kernel = std::string(),
                           const std::string& solver = std::string("cholesky_inv"));
      DensityFittingParams(StateIn&);
      ~DensityFittingParams();
      void save_data_state(StateOut&);

      const Ref<GaussianBasisSet>& basis() const { return basis_; }
      const std::string& kernel_key() const { return kernel_; }
      DensityFitting::SolveMethod solver() const { return solver_; }

      void print(std::ostream& o) const;
    private:
      static ClassDesc class_desc_;

      Ref<GaussianBasisSet> basis_;
      std::string kernel_;
      DensityFitting::SolveMethod solver_;
  };

   inline bool operator==(const DensityFittingParams& A, const DensityFittingParams& B) {
    return A.basis()->equiv(B.basis()) && A.kernel_key() == B.kernel_key() && A.solver() == B.solver();
  }

  /// this class encapsulates objects needed to perform density fitting of a 4-center integral
  struct DensityFittingInfo : virtual public SavableState {
    public:
      typedef DensityFittingInfo this_type;

      DensityFittingInfo(const Ref<DensityFittingParams>& p,
                         const Ref<DensityFittingRuntime>& r) :
                           params_(p), runtime_(r) {}
      DensityFittingInfo(StateIn& si) {
        params_ << SavableState::restore_state(si);
        runtime_ << SavableState::restore_state(si);
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(params_.pointer(),so);
        SavableState::save_state(runtime_.pointer(),so);
      }

      const Ref<DensityFittingParams>& params() const { return params_; }
      const Ref<DensityFittingRuntime>& runtime() const { return runtime_; }

      /// obsoletes all dependent runtimes
      void obsolete();

    private:
      Ref<DensityFittingParams> params_;
      Ref<DensityFittingRuntime> runtime_;

      static ClassDesc class_desc_;
  };

  inline bool operator==(const DensityFittingInfo& A, const DensityFittingInfo& B) {
    // dfinfo objects are equivalent if they use equivalent params and same runtime object
    return (*A.params() == *B.params()) && A.runtime() == B.runtime();
  }


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
