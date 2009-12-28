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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_df_df_runtime_h
#define _mpqc_src_lib_chemistry_qc_df_df_runtime_h

#include <chemistry/qc/df/df.h>
#include <chemistry/qc/mbptr12/tbint_runtime.h>

namespace sc {

  /** Parsed representation of a string key that represents fitting of a product of space1 and space2 into fspace
      Coulomb fitting kernel is the default. */
  class ParsedDensityFittingKey {
    public:
      ParsedDensityFittingKey(const std::string& key);

      const std::string& key() const { return key_; }
      const std::string& space1() const { return space1_; }
      const std::string& space2() const { return space2_; }
      const std::string& fspace() const { return fspace_; }

      /// computes key from its components
      static std::string key(const std::string& space1,
                             const std::string& space2,
                             const std::string& fspace);

    private:
      std::string key_;
      std::string space1_, space2_, fspace_;
  };


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
      DensityFittingRuntime(const Ref<MOIntsRuntime>& moints_runtime);
      DensityFittingRuntime(StateIn& si);
      void save_data_state(StateOut& so);

      /// obsoletes this object
      void obsolete();

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

    private:
      Ref<MOIntsRuntime> moints_runtime_;

      typedef Registry<std::string, ResultRef, detail::NonsingletonCreationPolicy > ResultRegistry;
      Ref<ResultRegistry> results_;

      // creates the result for a given key
      const ResultRef& create_result(const std::string& key);

      static ClassDesc class_desc_;

  };

  /// this class encapsulates objects needed to perform density fitting of a 4-center integral
  struct DensityFittingInfo : virtual public SavableState {
    public:
      typedef DensityFittingInfo this_type;

      DensityFittingInfo(const Ref<GaussianBasisSet>& b1,
                         const Ref<GaussianBasisSet>& b2,
                         const Ref<DensityFittingRuntime>& r) :
                           basis1_(b1), basis2_(b2), runtime_(r) {}
      DensityFittingInfo(StateIn& si) {
        basis1_ << SavableState::restore_state(si);
        basis2_ << SavableState::restore_state(si);
        runtime_ << SavableState::restore_state(si);
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(basis1_.pointer(),so);
        SavableState::save_state(basis2_.pointer(),so);
        SavableState::save_state(runtime_.pointer(),so);
      }

      const Ref<GaussianBasisSet>& basis1() const { return basis1_; }
      const Ref<GaussianBasisSet>& basis2() const { return basis2_; }
      const Ref<DensityFittingRuntime>& runtime() const { return runtime_; }

    private:
      Ref<GaussianBasisSet> basis1_;
      Ref<GaussianBasisSet> basis2_;
      Ref<DensityFittingRuntime> runtime_;

      static ClassDesc class_desc_;
  };



} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
