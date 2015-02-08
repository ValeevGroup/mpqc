//
// fockbuilder.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_fockbuilder_h
#define _mpqc_src_lib_chemistry_qc_lcao_fockbuilder_h

#include <cassert>
#include <util/ref/ref.h>
#include <util/group/thread.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/gpetite.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/lcao/fockbuild.h>
#include <chemistry/qc/lcao/clhfcontrib.h>
#include <chemistry/qc/lcao/hsoshfcontrib.h>
#include <chemistry/qc/lcao/moints_runtime.h>
#include <chemistry/qc/lcao/fockbuild_runtime.h>

namespace sc {

  namespace detail {
    template<bool bra_eq_ket> struct FockMatrixType;
    template<> struct FockMatrixType<true> {
      typedef RefSymmSCMatrix value;
      struct Factory {
        static value create(const Ref<SCMatrixKit>& kit, const RefSCDimension& bradim, const RefSCDimension& ketdim) {
          return kit->symmmatrix(bradim);
        }
        static void transform(const value& result, const value& original,
                               const RefSCMatrix& Ubra, const RefSCMatrix& Uket,
                               SCMatrix::Transform transform_type = SCMatrix::NormalTransform) {
          result.assign(0.0);
          result.accumulate_transform(Ubra,original,transform_type);
        }
        static void convert(const value& result, const value& original) {
          return result->convert(original);
        }
        template <typename T> static void convert(const value& result, const T& original) {
          abort();
        }
      };
    };
    template<> struct FockMatrixType<false> {
      typedef RefSCMatrix value;
      struct Factory {
        static value create(const Ref<SCMatrixKit>& kit, const RefSCDimension& bradim, const RefSCDimension& ketdim) {
          return kit->matrix(bradim,ketdim);
        }
        static void transform(const value& result, const value& original,
                               const RefSCMatrix& Ubra, const RefSCMatrix& Uket,
                               SCMatrix::Transform transform_type = SCMatrix::NormalTransform) {
          if (transform_type == SCMatrix::NormalTransform)
            result.accumulate_product(Ubra, original * Uket.t());
          else
            result.accumulate_product(Ubra.t(), original * Uket);
        }
        static void convert(const value& result, const value& original) {
          return result->convert(original);
        }
        template <typename T> static void convert(const value& result, const T& original) {
          abort();
        }
      };
    };

    typedef Ref<OneBodyInt> (Integral::*OneBodyIntCreator)();
    typedef Ref<OneBodyInt> (Integral::*MultipoleOneBodyIntCreator)(const Ref<IntParamsOrigin>&);

    template <OneBodyIntCreator OpEval>
    RefSymmSCMatrix onebodyint(const Ref<GaussianBasisSet>& bas,
                               const Ref<Integral>& integral) {
      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      // form skeleton operator matrix in AO basis
      RefSymmSCMatrix opmat_ao(bas->basisdim(), bas->matrixkit());
      opmat_ao.assign(0.0);
      Ref<SCElementOp> elemop =
          new OneBodyIntOp(new SymmOneBodyIntIter( (localints->*OpEval)(), pl));
      opmat_ao.element_op(elemop);
      elemop= 0;

      // now symmetrize
      RefSymmSCMatrix opmat(pl->SO_basisdim(), bas->so_matrixkit());
      pl->symmetrize(opmat_ao, opmat);

      return opmat;
    }

    template <OneBodyIntCreator OpEval>
    RefSymmSCMatrix onebodyint_ao(const Ref<GaussianBasisSet>& bas,
                                  const Ref<Integral>& integral) {
      RefSymmSCMatrix opmat_so = onebodyint<OpEval>(bas, integral);
      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      RefSymmSCMatrix opmat_ao = pl->to_AO_basis(opmat_so);
      return opmat_ao;
    }

    template <OneBodyIntCreator OpEval>
    RefSCMatrix onebodyint_ao(const Ref<GaussianBasisSet>& brabas,
                           const Ref<GaussianBasisSet>& ketbas,
                           const Ref<Integral>& integral) {

      Ref<Integral> localints = integral->clone();
      Ref<GPetiteList2> pl12 = GPetiteListFactory::plist2(brabas,ketbas);
      localints->set_basis(brabas,ketbas);
      Ref<OneBodyInt> ov_ints = (localints->*OpEval)();

      // form overlap in AO basis
      RefSCMatrix sao(brabas->basisdim(), ketbas->basisdim(), brabas->matrixkit());
      sao.assign(0.0);

      // TODO make more efficient and parallelizable by switching to operators
      const Ref<GaussianBasisSet>& bs1 = brabas;
      const Ref<GaussianBasisSet>& bs2 = ketbas;
      const int nshell1 = bs1->nshell();
      const int nshell2 = bs2->nshell();
      for(int sh1=0; sh1<nshell1; sh1++) {
        const int bf1_offset = bs1->shell_to_function(sh1);
        const int nbf1 = bs1->shell(sh1).nfunction();

        for(int sh2=0; sh2<nshell2; sh2++) {
          const int bf2_offset = bs2->shell_to_function(sh2);
          const int nbf2 = bs2->shell(sh2).nfunction();

          ov_ints->compute_shell(sh1,sh2);
          const double *ovintsptr = ov_ints->buffer();

          int bf1_index = bf1_offset;
          for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, ovintsptr+=nbf2) {
            int bf2_index = bf2_offset;
            const double *ptr = ovintsptr;
            for(int bf2=0; bf2<nbf2; bf2++, bf2_index++) {

              sao.set_element(bf1_index, bf2_index, *(ptr++));

            }
          }
        }
      }

      return sao;

    }

    template <MultipoleOneBodyIntCreator OpEval>
    void onebodyint_ao(const Ref<GaussianBasisSet>& brabas,
                       const Ref<GaussianBasisSet>& ketbas,
                       const Ref<Integral>& integral, const Ref<IntParamsOrigin>& multipole_origin,
                       std::vector<RefSCMatrix>& result) {

      Ref<Integral> localints = integral->clone();
      Ref<GPetiteList2> pl12 = GPetiteListFactory::plist2(brabas,ketbas);
      localints->set_basis(brabas,ketbas);
      Ref<OneBodyInt> ob_ints = (localints->*OpEval)(multipole_origin);

      // form obints in AO basis
      const size_t nopers = result.size();
      for(int i=0; i<nopers; ++i) {
        result[i] = brabas->matrixkit()->matrix(brabas->basisdim(), ketbas->basisdim());
        result[i].assign(0.0);
      }

      const Ref<GaussianBasisSet>& bs1 = brabas;
      const Ref<GaussianBasisSet>& bs2 = ketbas;
      const int nshell1 = bs1->nshell();
      const int nshell2 = bs2->nshell();
      for(int sh1=0; sh1<nshell1; sh1++) {
        const int bf1_offset = bs1->shell_to_function(sh1);
        const int nbf1 = bs1->shell(sh1).nfunction();

        for(int sh2=0; sh2<nshell2; sh2++) {
          const int bf2_offset = bs2->shell_to_function(sh2);
          const int nbf2 = bs2->shell(sh2).nfunction();

          ob_ints->compute_shell(sh1,sh2);
          const double *ovintsptr = ob_ints->buffer();

          int bf1_index = bf1_offset;
          for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, ovintsptr+=nopers*nbf2) {
            int bf2_index = bf2_offset;
            const double *ptr = ovintsptr;
            for(int bf2=0; bf2<nbf2; bf2++, bf2_index++) {

              for(size_t oper=0; oper<nopers; ++oper)
                result[oper].set_element(bf1_index, bf2_index, *(ptr++));

            }
          }
        }
      }

    }

    template <OneBodyIntCreator OpEval>
    RefSCMatrix onebodyint(const Ref<GaussianBasisSet>& brabas,
                           const Ref<GaussianBasisSet>& ketbas,
                           const Ref<Integral>& integral) {

      RefSCMatrix sao = onebodyint_ao<OpEval>(brabas, ketbas, integral);
      Ref<Integral> braints = integral->clone();  braints->set_basis(brabas);
      Ref<PetiteList> brapl = braints->petite_list();
      Ref<Integral> ketints = integral->clone();  ketints->set_basis(ketbas);
      Ref<PetiteList> ketpl = ketints->petite_list();

      // SO basis is always blocked, so first make sure a is blocked
      RefSCMatrix sao_blk = dynamic_cast<BlockedSCMatrix*>(sao.pointer());
      if (sao_blk.null()) {
        sao_blk = brabas->so_matrixkit()->matrix(brapl->AO_basisdim(),ketpl->AO_basisdim());
        sao_blk->convert(sao);
      }
      // if C1, then do nothing
      if (brabas->molecule()->point_group()->order() == 1)
        return sao_blk;

      // convert to SO basis
      RefSCMatrix ssoao = brapl->aotoso().t() * sao_blk;
      RefSCMatrix result = ssoao * ketpl->aotoso();

      return result;
    }

    /// computes overlap matrix in SO basis
    RefSymmSCMatrix overlap(const Ref<GaussianBasisSet>& bas,
                            const Ref<Integral>& integral);
    /// computes overlap matrix in SO basis
    RefSCMatrix overlap(const Ref<GaussianBasisSet>& brabs,
                        const Ref<GaussianBasisSet>& ketbs,
                        const Ref<Integral>& integral);
    /// computes electric dipole matrix in SO basis (minus sign not included)
    void edipole(const Ref<GaussianBasisSet>& basis,
                 const Ref<Integral>& integral,
                 RefSymmSCMatrix& mu_x,
                 RefSymmSCMatrix& mu_y,
                 RefSymmSCMatrix& mu_z
                 );
    /// computes 2-center 2-body Coulomb operator matrix in SO matrix
    RefSymmSCMatrix twobody_twocenter_coulomb(const Ref<GaussianBasisSet>& bas,
                                              const Ref<Integral>& integral);
    /// computes nonrelativistic Hamiltonian matrix in SO basis
    RefSymmSCMatrix nonrelativistic(const Ref<GaussianBasisSet>& bas,
                                    const Ref<Integral>& integral);
    /// computes Pauli Hamiltonian matrix (= nonrelativistic + MV + 1e Darwin) in SO basis
    RefSymmSCMatrix pauli(const Ref<GaussianBasisSet>& bas,
                          const Ref<Integral>& integral);
    /// computes DK Hamiltonian matrix in SO basis
    RefSymmSCMatrix dk(int dklev, const Ref<GaussianBasisSet>& bs,
                       const Ref<GaussianBasisSet>& pbs,
                       const Ref<Integral>& integral);

    RefSCMatrix coulomb(const Ref<TwoBodyFourCenterMOIntsRuntime>& ints_rtime,
                        const Ref<OrbitalSpace>& occspace,
                        const Ref<OrbitalSpace>& braspace,
                        const Ref<OrbitalSpace>& ketspace);
    RefSCMatrix exchange(const Ref<TwoBodyFourCenterMOIntsRuntime>& ints_rtime,
                         const Ref<OrbitalSpace>& occspace,
                         const Ref<OrbitalSpace>& braspace,
                         const Ref<OrbitalSpace>& ketspace);

    RefSCMatrix coulomb_df(const Ref<DensityFittingInfo>& df_info,
                           const RefSymmSCMatrix& P,
                           const Ref<GaussianBasisSet>& brabs,
                           const Ref<GaussianBasisSet>& ketbs,
                           const Ref<GaussianBasisSet>& obs);
#ifdef MPQC_NEW_RUNTIME
    RefSCMatrix coulomb_df_local(const Ref<DensityFittingInfo>& df_info,
                           const RefSymmSCMatrix& P,
                           const Ref<GaussianBasisSet>& brabs,
                           const Ref<GaussianBasisSet>& ketbs,
                           const Ref<GaussianBasisSet>& obs);
#endif // MPQC_NEW_RUNTIME
    RefSCMatrix exchange_df(const Ref<DensityFittingInfo>& df_info,
                            const RefSymmSCMatrix& P,
                            SpinCase1 spin,
                            const Ref<GaussianBasisSet>& brabs,
                            const Ref<GaussianBasisSet>& ketbs,
                            const Ref<GaussianBasisSet>& obs,
                            const Ref<FockBuildRuntime::PSqrtRegistry>& psqrtregistry);
#ifdef MPQC_NEW_RUNTIME
    RefSCMatrix exchange_df_local(const Ref<DensityFittingInfo>& df_info,
                            const RefSymmSCMatrix& P,
                            SpinCase1 spin,
                            const Ref<GaussianBasisSet>& brabs,
                            const Ref<GaussianBasisSet>& ketbs,
                            const Ref<GaussianBasisSet>& obs,
                            const Ref<FockBuildRuntime::PSqrtRegistry>& psqrtregistry);
#endif // MPQC_NEW_RUNTIME
  } // end of namespace detail


  /// Builds the two-body part of the Fock matrix in AO basis using integral-direct algorithm
  template<bool bra_eq_ket> class TwoBodyFockMatrixBuilder: public RefCount {

    public:

      typedef typename detail::FockMatrixType<bra_eq_ket>::value ResultType;
      typedef typename detail::FockMatrixType<bra_eq_ket>::Factory ResultFactory;

      TwoBodyFockMatrixBuilder(bool compute_F,
                        bool compute_J,
                        bool compute_K,
                        const Ref<GaussianBasisSet>& brabasis,
                        const Ref<GaussianBasisSet>& ketbasis,
                        const Ref<GaussianBasisSet>& densitybasis,
                        const RefSymmSCMatrix& density,
                        const RefSymmSCMatrix& openshelldensity,
                        const Ref<Integral>& integral,
                        const Ref<MessageGrp>& msg,
                        const Ref<ThreadGrp>& thr,
                        double accuracy = 1e-20) :
                          compute_F_(compute_F),
                          compute_J_(compute_J),
                          compute_K_(compute_K)
      {

        if (brabasis->equiv(ketbasis) != bra_eq_ket)
          throw ProgrammingError("FockMatrixBuilder::FockMatrixBuilder -- inconsistent constructor and template arguments",
                                 __FILE__, __LINE__);

        Ref<Integral> localints = integral->clone();

        Ref<FockContribution> fc;
        const bool openshell = openshelldensity;
        if (openshell) {
          fc = new HSOSHFContribution(brabasis, ketbasis, densitybasis, std::string("replicated"));
          ntypes_ = 2;
          fc->set_pmat(0, density);
          fc->set_pmat(1, openshelldensity);
        } else {
          fc = new CLHFContribution(brabasis, ketbasis, densitybasis, std::string("replicated"));
          ntypes_ = 1;
          fc->set_pmat(0, density);
        }

        // FockBuild can compute either J and K separately, or the total F.
        // Determine whether we really need to compute J and K, or F
        const bool really_compute_F = compute_F_ && !compute_J_ && !compute_K_;
        const bool really_compute_J = !really_compute_F && compute_J_;
        const bool really_compute_K = !really_compute_F && compute_K_;

        // FockBuild only accepts matrices created with appropriate basisdim(), which are not blocked, hence must
        // use non-blocked kit
        ResultType G_ao_skel[2][3];
        for(int t=0; t<ntypes_; ++t) {
          if (really_compute_J) {
            G_ao_skel[t][0] = ResultFactory::create(brabasis->matrixkit(),brabasis->basisdim(),ketbasis->basisdim());
            G_ao_skel[t][0].assign(0.0);
            fc->set_jmat(t, G_ao_skel[t][0]);
          }
          if (really_compute_K) {
            G_ao_skel[t][1] = ResultFactory::create(brabasis->matrixkit(),brabasis->basisdim(),ketbasis->basisdim());
            G_ao_skel[t][1].assign(0.0);
            fc->set_kmat(t, G_ao_skel[t][1]);
          }
          if (really_compute_F) {
            G_ao_skel[t][2] = ResultFactory::create(brabasis->matrixkit(),brabasis->basisdim(),ketbasis->basisdim());
            G_ao_skel[t][2].assign(0.0);
            fc->set_fmat(t, G_ao_skel[t][2]);
          }
        }

        const bool prefetch_blocks = false;
        Ref<FockDistribution> fd = new FockDistribution;
        fb_ = new FockBuild(fd, fc, prefetch_blocks, brabasis, ketbasis, densitybasis, msg, thr, localints);
        fb_->set_accuracy(accuracy);
        fb_->set_compute_J(really_compute_J);
        fb_->set_compute_K(really_compute_K);
        fb_->build();

        Ref<GPetiteList2> pl = GPetiteListFactory::plist2(brabasis,ketbasis);
        localints->set_basis(brabasis);
        Ref<PetiteList> brapl = localints->petite_list();
        localints->set_basis(ketbasis);
        Ref<PetiteList> ketpl = localints->petite_list();
        const int ng = pl->point_group()->char_table().order();

        /// convert skeleton matrices computed by FockBuild to the full matrices
        for(int t=0; t<ntypes_; ++t) {
          for(int c=0; c<3; ++c) {

            if (really_compute_J == false && c == 0) continue;
            if (really_compute_K == false && c == 1) continue;
            if (really_compute_F == false && c == 2) continue;
            // J matrices are spin-independent
            if (c == 0 && t == 1) continue;

            // if C1 -- nothing else needs to be done, return the result
            // same holds for brabasis != ketbasis since FockBuild does not use symmetry in that case
            if (ng == 1 || !bra_eq_ket) {
              RefSCDimension braaodim = brapl->AO_basisdim();
              RefSCDimension ketaodim = ketpl->AO_basisdim();
              result_[t][c] = ResultFactory::create(brabasis->so_matrixkit(),brapl->AO_basisdim(),ketpl->AO_basisdim());
              result_[t][c]->convert(G_ao_skel[t][c]);
            }
            else { // if not C1 and , symmetrize the skeleton G matrix to produce the SO basis G matrix
              G_ao_skel[t][c].scale(1.0/(double)ng);
              ResultType G_so = ResultFactory::create(brabasis->so_matrixkit(),brapl->SO_basisdim(),ketpl->SO_basisdim());
              symmetrize(pl,localints,G_ao_skel[t][c],G_so);
              G_ao_skel[t][c] = 0;

              // and convert back to AO basis, but this time make a blocked matrix
              RefSCDimension braaodim = brapl->AO_basisdim();
              RefSCDimension ketaodim = ketpl->AO_basisdim();
              result_[t][c] = ResultFactory::create(brabasis->so_matrixkit(),brapl->AO_basisdim(),ketpl->AO_basisdim());
              ResultFactory::transform(result_[t][c], G_so, brapl->sotoao(), ketpl->sotoao(), SCMatrix::TransposeTransform);

            }

            // FockBuild computes -K, hence change the sign
            if (c == 1) result_[t][1].scale(-1.0);

          }
        }

        if (compute_F_ && !really_compute_F) {
          // F = J - K = J - 1/2 (Ka + Kb)
          result_[0][2] = result_[0][0] - result_[0][1];
          if (ntypes_ == 2) { // Fo = -Ko = -1/2 (Ka - Kb)
            result_[1][2] = result_[1][1].copy();
            result_[1][2].scale(-1.0);
          }
        }

      }

      const Ref<FockBuild>& builder() const { return fb_; }
      double nints() const { return builder()->contrib()->nint(); }
      const ResultType& F(unsigned int t = 0) const {
        MPQC_ASSERT(compute_F_ && t < ntypes_);
        return result_[t][2];
      }
      const ResultType& J(unsigned int t = 0) const {
        MPQC_ASSERT(compute_J_ && t < ntypes_ && t == 0);
        return result_[t][0];
      }
      const ResultType& K(unsigned int t = 0) const {
        MPQC_ASSERT(compute_K_ && t < ntypes_);
        return result_[t][1];
      }
      ResultType F(SpinCase1 spin) const {
        if (ntypes_ == 1)
          return F(0);
        else {
          return (spin == Alpha) ? F(0) + F(1) : F(0) - F(1);
        }
      }
      ResultType K(SpinCase1 spin) const {
        if (ntypes_ == 1)
          return K(0);
        else {
          return (spin == Alpha) ? K(0) + K(1) : K(0) - K(1);
        }
      }

    private:

      Ref<FockBuild> fb_;
      bool compute_J_;
      bool compute_K_;
      bool compute_F_;
      int ntypes_;
      ResultType result_[2][3];

  }; // class TwoBodyFockMatrixBuilder

  /// Builds the two-body part of the Fock matrix in MO basis using AO->MO transforms
  class TwoBodyFockMatrixTransformBuilder: public RefCount {

    public:

      TwoBodyFockMatrixTransformBuilder(bool compute_J,
                        bool compute_K,
                        SpinCase1 spin,
                        const Ref<OrbitalSpace>& braspace,
                        const Ref<OrbitalSpace>& ketspace,
                        const Ref<OrbitalSpace>& occspace_A,
                        const Ref<OrbitalSpace>& occspace_B,
                        const Ref<TwoBodyFourCenterMOIntsRuntime>& ints_rtime) :
                          compute_J_(compute_J),
                          compute_K_(compute_K),
                          compute_F_(compute_J && compute_K)
      {

        for(int t=0; t<3; ++t) {

          if (compute_J == false && t == 0) continue;
          if (compute_K == false && t == 1) continue;

          if (t == 0) { // coulomb
            result_[t] = detail::coulomb(ints_rtime,occspace_A,braspace,ketspace);
            if (*occspace_A == *occspace_B) { // alpha and beta spins equivalent? Scale by 2, else compute
              result_[t].scale(2.0);
            }
            else {
              result_[t].accumulate( detail::coulomb(ints_rtime,occspace_B,braspace,ketspace) );
            }
          }

          if (t == 1) { // exchange
            result_[t] = detail::exchange(ints_rtime,(spin == Alpha ? occspace_A : occspace_B),braspace,ketspace);
          }

        }

        if (compute_F_)
          result_[2] = result_[0] - result_[1];

      }

      const RefSCMatrix& F() const {
        MPQC_ASSERT(compute_F_);
        return result_[2];
      }
      const RefSCMatrix& J() const {
        MPQC_ASSERT(compute_J_);
        return result_[0];
      }
      const RefSCMatrix& K() const {
        MPQC_ASSERT(compute_K_);
        return result_[1];
      }

    private:

      bool compute_J_;
      bool compute_K_;
      bool compute_F_;
      RefSCMatrix result_[3];

  }; // class TwoBodyFockMatrixTransformBuilder

  /// Builds the two-body part of the Fock matrix in AO basis using DF-based algorithm
  class TwoBodyFockMatrixDFBuilder: public RefCount {

    typedef RefSCMatrix ResultType;

    public:

      TwoBodyFockMatrixDFBuilder(bool compute_F,
                                 bool compute_J,
                                 bool compute_K,
                           const Ref<GaussianBasisSet>& brabasis,
                           const Ref<GaussianBasisSet>& ketbasis,
                           const Ref<GaussianBasisSet>& densitybasis,
                           const RefSymmSCMatrix& density,
                           const RefSymmSCMatrix& openshelldensity,
                           const Ref<DensityFittingInfo>& df_info,
                           const Ref<FockBuildRuntime::PSqrtRegistry>& psqrtregistry) :
                          compute_J_(compute_J),
                          compute_K_(compute_K),
                          compute_F_(compute_F),
                          ntypes_(openshelldensity ? 2 : 1)
      {

        // DF-based builds are separate for J and K separately
        // Determine whether we really need to compute J and K, or F
        const bool really_compute_F = compute_F_;
        const bool really_compute_J = compute_J_ || compute_F_;
        const bool really_compute_K = compute_K_ || compute_F_;

        for(int t=0; t<ntypes_; ++t) {
          const SpinCase1 spin = static_cast<SpinCase1>(t);
          for(int c=0; c<3; ++c) {

            if (really_compute_J == false && c == 0) continue;
            if (really_compute_K == false && c == 1) continue;
            if (really_compute_F == false && c == 2) continue;
            // J matrices are spin-independent
            if (c == 0 && t == 1) continue;

            if (c == 0) { // coulomb
              if(df_info->params()->local_coulomb()) {
#ifdef MPQC_NEW_RUNTIME
                result_[0][c] = detail::coulomb_df_local(df_info, density, brabasis, ketbasis, densitybasis);
#else // MPQC_NEW_RUNTIME
                throw sc::FeatureNotImplemented("Local DF Coulomb build requires MPQC3 runtime: compile in boost, tiledarray, etc.",
                                                __FILE__,__LINE__);
#endif // MPQC_NEW_RUNTIME
              }
              else
                result_[0][c] = detail::coulomb_df(df_info, density, brabasis, ketbasis, densitybasis);
            }

            if (c == 1) { // exchange
              RefSymmSCMatrix Pspin;
              SpinCase1 spincase;
              if (openshelldensity) {
                Pspin = (spin == Alpha) ? density + openshelldensity : density - openshelldensity;
                spincase = spin;
              }
              else {
                Pspin = density.copy();
                spincase = AnySpinCase1;
              }
              Pspin.scale(0.5);
              if(df_info->params()->local_exchange()) {
#ifdef MPQC_NEW_RUNTIME
                result_[t][c] = detail::exchange_df_local(df_info, Pspin, spincase, brabasis, ketbasis, densitybasis,
                                                    psqrtregistry);
#else // MPQC_NEW_RUNTIME
                throw sc::FeatureNotImplemented("Local DF exchange build requires MPQC3 runtime: compile in boost, tiledarray, etc.",
                                                __FILE__,__LINE__);
#endif // MPQC_NEW_RUNTIME
              }
              else
                result_[t][c] = detail::exchange_df(df_info, Pspin, spincase, brabasis, ketbasis, densitybasis,
                                                    psqrtregistry);
            }

            if (c == 2) { // fock
              result_[t][2] = result_[0][0] - result_[t][1];
            }
          }
        }

      }

      const ResultType& F(unsigned int t = 0) const {
        MPQC_ASSERT(compute_F_ && t < ntypes_);
        return result_[t][2];
      }
      const ResultType& J(unsigned int t = 0) const {
        MPQC_ASSERT(compute_J_ && t < ntypes_);
        return result_[t][0];
      }
      const ResultType& K(unsigned int t = 0) const {
        MPQC_ASSERT(compute_K_ && t < ntypes_);
        return result_[t][1];
      }
      ResultType F(SpinCase1 spin) const {
        if (ntypes_ == 1)
          return F(0);
        else {
          return (spin == Alpha) ? F(0) : F(1);
        }
      }
      ResultType K(SpinCase1 spin) const {
        if (ntypes_ == 1)
          return K(0);
        else {
          // DF-based computed alpha and beta exchange matrices
          return (spin == Alpha) ? K(0) : K(1);
        }
      }

    private:

      bool compute_J_;
      bool compute_K_;
      bool compute_F_;
      int ntypes_;
      ResultType result_[2][3];

  }; // class TwoBodyFockMatrixDFBuilder

  /// Builds the one-body part of the Fock matrix in AO basis
  template<bool bra_eq_ket> class OneBodyFockMatrixBuilder: public RefCount {

    public:

      typedef enum { NonRelativistic, Pauli, DK1, DK2 } OneBodyHamiltonianType;

      typedef typename detail::FockMatrixType<bra_eq_ket>::value ResultType;
      typedef typename detail::FockMatrixType<bra_eq_ket>::Factory ResultFactory;

      OneBodyFockMatrixBuilder(OneBodyHamiltonianType type,
                        const Ref<GaussianBasisSet>& brabasis,
                        const Ref<GaussianBasisSet>& ketbasis,
                        const Ref<GaussianBasisSet>& pbasis,
                        const Ref<Integral>& integral,
                        double accuracy = 1e-12)
      {
        if (brabasis->equiv(ketbasis) != bra_eq_ket)
          throw ProgrammingError("OneBodyFockMatrixBuilder::OneBodyFockMatrixBuilder -- inconsistent constructor and template arguments",
                                 __FILE__, __LINE__);

        const Ref<GaussianBasisSet>& bs1 = brabasis;
        const Ref<GaussianBasisSet>& bs2 = ketbasis;
        const bool bs1_eq_bs2 = bra_eq_ket;

        Ref<GaussianBasisSet> hcore_basis = (bs1_eq_bs2) ? bs1 : bs1+bs2;

        Ref<GaussianBasisSet> p_basis = pbasis;
        if (type == DK1 || type == DK2) {
          // momentum basis in DKH calculations must span both bra and ket basis sets.
          // easiest to achieve if both are included in p_basis
          p_basis = p_basis + hcore_basis;
        }

        RefSymmSCMatrix hsymm;
        switch (type) {
          case NonRelativistic: hsymm = detail::nonrelativistic(hcore_basis,integral); break;
          case Pauli: hsymm = detail::pauli(hcore_basis,integral); break;

          int dklev;
          case DK1:  dklev = 1;
          case DK2:  dklev = 2;
          hsymm = detail::dk(dklev, hcore_basis, p_basis, integral);
          break;

          default:
            throw ProgrammingError("Unrecognized Hamiltonian type",__FILE__,__LINE__);
        }

        // convert hsymm to the AO basis
        Ref<Integral> localints = integral->clone();
        localints->set_basis(hcore_basis,hcore_basis);
        Ref<PetiteList> hcore_pl = localints->petite_list();
        RefSymmSCMatrix hsymm_ao = hcore_pl->to_AO_basis(hsymm);
        hsymm = 0;

        localints->set_basis(brabasis);
        Ref<PetiteList> brapl = localints->petite_list();
        localints->set_basis(ketbasis);
        Ref<PetiteList> ketpl = localints->petite_list();
        RefSCDimension braaodim = brapl->AO_basisdim();
        RefSCDimension ketaodim = ketpl->AO_basisdim();
        if (bs1_eq_bs2) {
          result_ = ResultFactory::create(brabasis->so_matrixkit(),braaodim,ketaodim);
          ResultFactory::convert(result_,hsymm_ao);
        }
        else {
          RefSCMatrix result_rect = brabasis->so_matrixkit()->matrix(braaodim,ketaodim);
          RefSCMatrix hrect_ao = brabasis->so_matrixkit()->matrix(hsymm_ao.dim(),hsymm_ao.dim());
          hrect_ao.assign(0.0);
          hrect_ao.accumulate(hsymm_ao);

          // extract the bs1 by bs2 block:
          //   loop over all bs1 fblocks
          //     loop over all bs2 fblocks
          //       copy block to h
          //     end loop
          //   end loop
          GaussianBasisSetMap bs1_to_hbs(bs1, hcore_basis);
          GaussianBasisSetMap bs2_to_hbs(bs2, hcore_basis);
          const int nfblock1 = bs1_to_hbs.nfblock();
          const int nfblock2 = bs2_to_hbs.nfblock();
          for (int rb = 0; rb < nfblock1; ++rb) {
            const int rf_start1 = bs1_to_hbs.fblock_to_function(rb);
            const int rf_end1 = rf_start1 + bs1_to_hbs.fblock_size(rb) - 1;

            const int rf_start12 = bs1_to_hbs.map_function(rf_start1);
            const int rf_end12 = bs1_to_hbs.map_function(rf_end1);

            for (int cb = 0; cb < nfblock2; ++cb) {
              const int cf_start2 = bs2_to_hbs.fblock_to_function(cb);
              const int cf_end2 = cf_start2 + bs2_to_hbs.fblock_size(cb) - 1;

              const int cf_start12 =
                  bs2_to_hbs.map_function(cf_start2);
              const int cf_end12 = bs2_to_hbs.map_function(cf_end2);

              // assign row-/col-subblock to h
              result_rect.assign_subblock(hrect_ao, rf_start1, rf_end1, cf_start2, cf_end2,
                                                    rf_start12, cf_start12);
            }
          }

          result_ = ResultFactory::create(brabasis->so_matrixkit(),braaodim,ketaodim);
          ResultFactory::convert(result_,result_rect);

        }

      }

      const ResultType& result() const { return result_; }

    private:

      ResultType result_;

  }; // class OneBodyFockMatrixBuilder

} // end of namespace sc

#endif // end of header guard

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
