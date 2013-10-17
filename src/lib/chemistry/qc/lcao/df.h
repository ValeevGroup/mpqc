//
// df.h
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


#ifndef _mpqc_src_lib_chemistry_qc_df_df_h
#define _mpqc_src_lib_chemistry_qc_df_df_h

#include <util/state/state.h>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <math/scmat/result.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/lcao/tbint_runtime.h>

namespace sc {

  class OrbitalSpace;
  class DistArray4;

  /** Decomposition by density fitting with respect to some kernel. The definition is as follows:
   \f[
   |rs) = \sum_A C^{rs}_A |A)  \quad ,
   \f]
   where \f$ |A) \f$ belong to a fitting basis and \f$ C^{rs}_A \f$ are determined by solving the linear system
   \f[
   (A|\hat{W}|rs) = \sum_B C^{rs}_B (A|\hat{W}|B)  \quad .
   \f]
   The kernel \f$\hat{W} \f$ is any (positive-definite) operator.

   \param runtime is the MOIntsRuntime object used to compute MO/AO integrals
   \param kernel_key denotes the kernel operator to be used for fitting. Currently only "1/r_{12}" is supported.
   \param space1 is space r in |rs) product space
   \param space2 is space s in |rs) product space
   \param fitting_basis is the basis set in which to perform the fitting

   */
  class DensityFitting: virtual public SavableState {
    public:
      enum SolveMethod {
        SolveMethod_InverseCholesky = 0,
        SolveMethod_Cholesky = 1,
        SolveMethod_RefinedCholesky = 2,
        SolveMethod_InverseBunchKaufman = 3,
        SolveMethod_BunchKaufman = 4,
        SolveMethod_RefinedBunchKaufman = 5
      };

      typedef TwoBodyMOIntsRuntimeUnion23 MOIntsRuntime;

      ~DensityFitting();
      /**
       * Determines density fitting of product (space1 space2) in terms of fitting basis.
       * @param rtime
       * @param kernel_key
       * @param solver_key
       * @param space1
       * @param space2
       * @param fitting_basis
       * @return
       */
      DensityFitting(const Ref<MOIntsRuntime>& rtime,
                     const std::string& kernel_key,
                     SolveMethod solver,
                     const Ref<OrbitalSpace>& space1,
                     const Ref<OrbitalSpace>& space2,
                     const Ref<GaussianBasisSet>& fitting_basis);
      DensityFitting(StateIn&);
      void save_data_state(StateOut&);

      const Ref<MOIntsRuntime>& runtime() const { return runtime_; }
      const Ref<Integral>& integral() const {
        return runtime()->factory()->integral();
      }
      const Ref<OrbitalSpace>& space1() const { return space1_; }
      const Ref<OrbitalSpace>& space2() const { return space2_; }
      const Ref<GaussianBasisSet>& fbasis() const { return fbasis_; }
      const std::string& kernel_key() const { return kernel_key_; }
      SolveMethod solver() const { return solver_; }
      /// returns the fitting coeffcients
      const Ref<DistArray4>& C() const {
        return C_;
      }
      /// produces RefSCDimension that corresponds to the product space
      RefSCDimension product_dimension() const;

      virtual void compute();

      /**
       * Test whether a density-fitting kernel has definite sign.
       * @note this is a preliminary implementation that makes common assumptions about the <tt>kernel_params</tt> argument
       * @param kernel_type kernel type; only some operators are accepted as density fitting kernel
       * @param kernel_params kernel parameters; must be in ParamsRegistry
       * @return +1 if the kernel is positive-definite, -1 if the kernel is negative-definite, 0 otherwise
       */
      static int definite_kernel(TwoBodyOper::type kernel_type, const Ref<IntParams>& kernel_params);
      static int definite_kernel(TwoBodyOper::type kernel_type, const std::string& kernel_params_key);

    private:

      Ref<MOIntsRuntime> runtime_;
      Ref<GaussianBasisSet> fbasis_;
      Ref<OrbitalSpace> space1_;
      Ref<OrbitalSpace> space2_;
      std::string kernel_key_;
      SolveMethod solver_;

      bool evaluated_;

      static ClassDesc class_desc_;

    protected:

      Ref<DistArray4> C_;
      RefSymmSCMatrix kernel_;
      Ref<DistArray4> cC_;

      /// returns the kernel matrix in the fitting basis, i.e. \f$ (A|\hat{W}|B) \f$ .
      const RefSymmSCMatrix& kernel() const {
        return kernel_;
      }
      /** returns the conjugate expansion matrix defined as \f$ \tilde{C}^{rs}_A \equiv (A|\hat{W}|rs) = \sum_B C^{rs}_B (A|\hat{W}|B) \f$ .
          \f$ \tilde{C} \f$ is conjugate to \f$ C \f$ in the sense that their
          product returns reconstructed integrals:
          \f$ (pq|rs) \approx \sum_A  \tilde{C}^A_{pq} C^{rs}_A \f$
        */
      const Ref<DistArray4>& conjugateC() const {
        return cC_;
      }

  };

  /// Computes density fitting for |ij) density from fitting of |iq) DensityFitting where q is the AO space supporting j
  class TransformedDensityFitting: public DensityFitting {
    public:
      typedef DensityFitting::MOIntsRuntime MOIntsRuntime;

      ~TransformedDensityFitting();
      /// compute density fitting for |mo1 mo2) from |mo1 ao2)
      TransformedDensityFitting(const Ref<MOIntsRuntime>& rtime,
                                const std::string& kernel_key,
                                SolveMethod solver,
                                const Ref<OrbitalSpace>& space1,
                                const Ref<OrbitalSpace>& space2,
                                const Ref<GaussianBasisSet>& fitting_basis,
                                const Ref<DistArray4>& mo1_ao2_df);
      TransformedDensityFitting(StateIn&);
      void save_data_state(StateOut&);

    private:

      Ref<DistArray4> mo1_ao2_df_;

      // overloads DensityFitting::compute()
      void compute();

      static ClassDesc class_desc_;

  };

  /// Computes density fitting for |ij) density from fitting of |ji) DensityFitting
  class PermutedDensityFitting: public DensityFitting {
    public:
      ~PermutedDensityFitting();
      /// compute density fitting for |mo1 mo2) from |mo1 ao2)
      PermutedDensityFitting(const Ref<DensityFitting>& df21);
      /// compute density fitting for |mo1 mo2) from |mo1 ao2)
      PermutedDensityFitting(const Ref<MOIntsRuntime>& rtime,
                             const std::string& kernel_key,
                             SolveMethod solver,
                             const Ref<OrbitalSpace>& space1,
                             const Ref<OrbitalSpace>& space2,
                             const Ref<GaussianBasisSet>& fitting_basis,
                             const Ref<DistArray4>& df21);
      PermutedDensityFitting(StateIn&);
      void save_data_state(StateOut&);

    private:
      Ref<DistArray4> df21_;

      // overloads DensityFitting::compute()
      void compute();

      static ClassDesc class_desc_;

  };

} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
