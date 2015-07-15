//
// singlereference_r12_intermediates.h
//
// Copyright (C) 2013 Edward Valeev
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
#pragma implementation
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediates_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediates_h

#if defined(MPQC_NEW_FEATURES)
# include <tiledarray.h>
#else
# error "sr_r12intermediates.h requires MPQC3 runtime, but it is not available"
#endif

#include <chemistry/qc/wfn/rdm.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <../bin/mpqc/mpqcinit.h>

namespace TA = TiledArray;

namespace sc {

  /// Tile of a 4-index tensor that's "evaluated" when needed by reading from DistArray4
  template <typename T>
  class DA4_Tile {
  private:
    TA::Array<T, 4, DA4_Tile>* owner_;
    std::array<std::size_t, 4> index_;
    Ref<DistArray4> darray4_;
    int te_type_;

  public:
    typedef T value_type;
    typedef TA::Tensor<T> eval_type;
    typedef typename eval_type::range_type range_type;

    DA4_Tile() { }

    DA4_Tile(TA::Array<T, 4, DA4_Tile>* owner,
             const std::array<std::size_t, 4>& index,
             const Ref<DistArray4>& darray4,
             int te_type) :
      owner_(owner), index_(index), darray4_(darray4), te_type_(te_type)
    { }

    operator TA::Tensor<T> () const;

    /// \tparam Archive The serialization archive type
    /// \param ar The serialization archive
    template <typename Archive>
    void serialize(Archive& ar) {
      MPQC_ASSERT(false);
    }

  };

  /// Tile of a <34> slice of <1234> that's "evaluated" when needed by reading from DistArray4 holding pqrs.
  template <typename T>
  class DA4_Tile34 {
  private:
    TA::Array<TA::Tensor<T>, 2, DA4_Tile34 >* owner_;
    std::array<std::size_t, 2> index_;
    Ref<DistArray4> darray4_;
    int te_type_;

  public:
    typedef T value_type;
    typedef TA::Tensor<TA::Tensor<T> > eval_type;
    typedef typename eval_type::range_type range_type;
    typedef T numeric_type;

    DA4_Tile34() { }

    DA4_Tile34(TA::Array<TA::Tensor<T>, 2, DA4_Tile34 >* owner,
               const std::array<std::size_t, 2>& index,
               const Ref<DistArray4>& darray4,
               int te_type) :
      owner_(owner), index_(index), darray4_(darray4), te_type_(te_type)
    { }

    /// conversion to eval_type
    operator TA::Tensor<TA::Tensor<T> > () const;

    /// \tparam Archive The serialization archive type
    /// \param ar The serialization archive
    template <typename Archive>
    void serialize(Archive& ar) {
      MPQC_ASSERT(false);
    }

  };

  /// Tile of a DIM-order tensor that's "evaluated" when needed by calling ElementGenerator({i0, i1, i2, .... i_DIM-1})
  template <typename T, unsigned int DIM, typename ElementGenerator>
  class LazyTensor {
  private:
    TA::Array<T, DIM, LazyTensor>* owner_;
    std::array<std::size_t, DIM> index_;
    ElementGenerator* element_generator_;

  public:
    typedef T value_type;
    typedef TA::Tensor<T> eval_type;
    typedef typename eval_type::range_type range_type;

    LazyTensor() { }

    LazyTensor(const LazyTensor& other) :
      owner_(other.owner_), index_(other.index_), element_generator_(other.element_generator_)
    {}

    LazyTensor& operator=(const LazyTensor& other) {
      owner_ = other.owner_;
      index_ = other.index_;
      element_generator_ = other.element_generator_;
      return *this;
    }

    LazyTensor(TA::Array<T, DIM, LazyTensor>* owner,
               const std::array<std::size_t, DIM>& index,
               ElementGenerator* gen) :
      owner_(owner), index_(index), element_generator_(gen)
    { }

    operator TA::Tensor<T> () const {

      eval_type tile(owner_->trange().make_tile_range(index_));

      auto* ptr = tile.data();
      for(auto i = tile.range().begin();
          i!=tile.range().end();
          ++i, ++ptr) {

        *ptr = (*element_generator_)(*i);

      }

      return tile;
    }

    /// \tparam Archive The serialization archive type
    /// \param ar The serialization archive
    template <typename Archive>
    void serialize(Archive& ar) {
      MPQC_ASSERT(false);
    }

  };

  namespace expressions {
    template <typename ArgType, bool Transpose> struct trace_tensor2_op;
    template <typename ArgType, bool Transpose> struct diag_tensor2_op;

    /// makes a geminal T tensor
    template <typename T>
    struct TGeminalGenerator {
        /// @param spin 0(singlet) or 1(triplet)
        TGeminalGenerator(unsigned int spin = 0) : s_(spin) {
          MPQC_ASSERT(s_ == 0 || s_ == 1);
        }

        template <typename Index> T operator()(const Index& i) {
          if (s_ == 0) {
            if (i[0] == i[1] && i[0] == i[2] && i[0] == i[3]) // iiii
              return 1.0/2.0;
            else if (i[0] == i[2] && i[1] == i[3]) // ijij
              return 3.0/8.0;
            else if (i[0] == i[3] && i[1] == i[2]) // ijji
              return 1.0/8.0;
          }
          if (s_ == 1) {
            if ((i[0] == i[2] && i[1] == i[3]) || (i[0] == i[3] && i[1] == i[2])) // ijij
              return 1.0/4.0;
          }
          return 0.0;
        }

      private:
        unsigned int s_;
    };

  };


  /// SingleReference_R12Intermediates computes R12/F12 intermediates using MPQC3 runtime
  /// @tparam T the numeric type supporting all tensors
  template <typename T>
  class SingleReference_R12Intermediates {
    public:

      /// standard 4-index tensor
      typedef TA::Array<T, 4 > TArray4; // Tile = Tensor<T>
      /// 4-index tensor with lazy tiles
      typedef TA::Array<T, 4, DA4_Tile<T> > TArray4d; // Tile = DA4_Tile<T>
      /// standard 2-index tensor
      typedef TA::Array<T, 2> TArray2; // Tile = Tensor<T>
      /// 2-index tensor with lazy tiles
      typedef TA::Array<T, 2, DA4_Tile<T> > TArray2d; // Tile = DA4_Tile<T>
      /// 2-index tensor with lazy tiles
      //typedef TA::Array<T, 2, LazyTensor<T, 2, ElementGenerator> > TArray2d; // Tile = LazyTensor<T, 2, ElementGenerator>
      /// 2-index tensor of 2-index tensors
      typedef TA::Array<TA::Tensor<T>, 2> TArray22; // Tile = Tensor<Tensor<T>>
      /// 2-index tensor of lazy 2-index tensors
      typedef TA::Array<TA::Tensor<T>, 2, DA4_Tile34<T> > TArray22d; // Tile = Tensor<Tensor<T>>

      /// geminal T tensor
      typedef TA::Array<T, 4, LazyTensor<T, 4, expressions::TGeminalGenerator<T> > > TArray4Tg;

      /**
       * Constructs an SingleReference_R12Intermediates object. There can be many.
       *
       * @param world MADNESS world in which the computed Arrays will live
       * @param r12world R12WavefunctionWorld that provides computation of integrals and OrbitalSpace objects
       */
      SingleReference_R12Intermediates(madness::World& world,
                                       const Ref<R12WavefunctionWorld>& r12world) :
                                         world_(world),
                                         r12world_(r12world)
      {
      }
      ~SingleReference_R12Intermediates() {
        r12world_ = 0;
      }

      const Ref<R12WavefunctionWorld>& r12world() const {
        return r12world_;
      }

      /** computes diagonal (spin-restricted, for now) V intermediate
      * @return \f$ V^{ij}_{ij} \f$ and \f$ V^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> V_diag();

      /** computes diagonal (spin-restricted, for now) X intermediate
      * @return \f$ X^{ij}_{ij} \f$ and \f$ X^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> X_diag();

      /** computes diagonal (spin-restricted, for now) B intermediate
      * @return \f$ B^{ij}_{ij} \f$ and \f$ B^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> B_diag();

      /**
       * Computes second-order Green's function IPs and EAs
       * \parame orbital the index of the orbital, -1 = HOMO, +1 = LUMO
       */
      void gf2_r12(int orbital);

      /** returns the 1-particle reduced density matrix
      * @return \f$ \gamma^{p}_{q} \f$, respectively
      */
      TArray2 rdm1();

      // compute multipole
      void compute_multipole();
      // compute B^P_Q which is summed over k
      TArray2 BPk_Qk(const char* p, const char* q,
                     const double C_0, const double C_1);
      TArray4 Bpr_qs(const char* p, const char* q);

      // compute V^Pq_Rs = 1/2 \bar{R}^Pq_AlphaBeta \bar{g}^AlphaBeta_Rs
      // P and R are in alpha or beta space
      TArray4 VPq_Rs(const char* p, const char* q,
                     const char* r, const char* s,
                     const double C_0, const double C_1);
      // V^Rk_Sk which is summed over k
      TArray2 VRk_Sk(const char* r, const char* s,
                     const double C_0, const double C_1);

      // compute Xam contribution from CABS Singles
      TArray2 Xam_CabsSingles(const TArray2& TmA, const TArray2& Tma);

      // compute Xam contribution from MP2
      TArray2 Xam_mp2(const TArray4& T2_ijab,
                      const TArray2& Dij, const TArray2& Dab);

      // compute Xam contribution from MP2 F12 coulping part
      TArray2 Xam_Cmp2f12(const double C_0, const double C_1,
                          const TArray4& T2_ijab, const TArray4& A_ijab,
                          const TArray2& Dij, const TArray2& Dab,
                          const TArray2& RT_apb);

      // compute Xam contribution from F12 V part
      TArray2 Xam_V(const double C_0, const double C_1);
      // compute Xam contribution from F12 X part
      TArray2 Xam_X(const double C_0, const double C_1);
      // compute Xam contribution from F12 B part
      TArray2 Xam_B(const double C_0, const double C_1);

      // compute Xii' contribution from F12 part for frozen-core systems
      TArray2 Xiip_VBX(const double C_0, const double C_1);
      // compute Xii' contribution from CCSD F12 coupling
      // for frozen-core systems
      TArray2 Xiip_CVT(const double C_0, const double C_1,
                       const TArray2& T1, const TArray4& T2);

      // Xam contribution of CT2 part resulted from CCSD F12 coupling
      TArray2 Xam_CT2_ccsd(const double C_0, const double C_1,
                           const TArray4& T2, const TArray2& RT2_aPb);
      // Xam contribution of VT(T1&T2) part resulted from CCSD F12 coupling
      TArray2 Xam_VT_ccsd(const double C_0, const double C_1,
                          const TArray2& T1, const TArray4& T2);

      // compute T1 & T2 amplitudes of CC2
      void compute_T_cc2(TArray2& T1, TArray4& T2);

      // compute Lambda_1 & Lambda_2 amplitudes of CC2
      // using formula from Schaefer III, JCP 87, 5361 (1987)
      void compute_lambda_cc2(const TArray2& t1, const TArray4& t2,
                         TArray2& L1, TArray4& L2);

      // use formula from Gauss and Stanton, JCP, 103 (1995)
      void compute_lambda_cc2_2(const TArray2& t1, const TArray4& t2,
                                TArray2& L1, TArray4& L2);

      // CCSD (Gauss and Stanton formula)
      // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
      //       & Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
      void compute_Delta_cc(const TArray2& fai, const TArray4d& g_abij,
                            TArray2& D_ai, TArray4& D_abij);
      // compute intermediates needed for T and lambda amplitudes
      void compute_intermediates_TFW_ccsd(
             const TArray2& t1, const TArray4& t2,
             const TArray4d& g_abij, const TArray4d& g_aikl,
             const TArray4d& g_aibj, const TArray4d& g_aijb,
             const TArray4d& g_abci,
             TArray2& TFkc, TArray2& TFac, TArray2& TFki,
             TArray4& TW_KbCj_ab,TArray4& TW_KbcJ_ba,
             TArray4& TW_AbCd_ab, TArray4& TW_KlIj_ab);
      void compute_intermediates_CW_ccsd(
             const TArray2& t1, const TArray4& t2,
             TArray2& CFkc, TArray2& CFac, TArray2& CFki,
             TArray4& CW_KbEj_ab, TArray4& CW_KbeJ_ba,
             TArray4& CW_AbCd_ab, TArray4& CW_AbCi_ab,
             TArray4& CW_KlIj_ab, TArray4& CW_KbIj_ab,
             TArray4& CW_KlIc_ab, TArray4& CW_KliC_ab,
             TArray4& CW_AkCd_ab);
      // compute CCSD T amplitudes
      void compute_T_ccsd(TArray2& t1, TArray4& t2);
      // compute CCSD lambda amplitudes
      void compute_lambda_ccsd(const TArray2& t1, const TArray4& t2,
                               TArray2& L1, TArray4& L2,
                               const bool include_F12contri);

      // compute CC2 one-electron density from amplitudes
      void compute_cc2_1rdm_amp(const TArray2& T1_cc2, const TArray4& T2_cc2,
                               const TArray2& L1_cc2, const TArray4& L2_cc2,
                               TArray2& Dij_cc2, TArray2& Dab_cc2,
                               TArray2& Dia_cc2, TArray2& Dai_cc2);

      // cpmute Gamma intermediates needed for CCSD Xam
      void compute_Gamma_ijab_ccsd(const TArray2& T1, const TArray4& T2,
                                   const TArray4& tau_ab,
                                   const TArray2& L1, const TArray4& L2,
                                   TArray4& Gamma_IjAb_ab);
      void compute_ccsd_1rdm_amp(const TArray2& T1, const TArray4& T2,
                                 const TArray2& L1, const TArray4& L2,
                                 TArray2& Dij, TArray2& Dab,
                                 TArray2& Dia, TArray2& Dai);

      void compute_Gamma2_ccsd(const TArray2& T1, const TArray4& T2,
                              const TArray2& L1, const TArray4& L2,
                              const TArray4& tau_aa, const TArray4& tau_ab,
                              TArray4& Gamma_IjKa_ab,
                              TArray4& Gamma_AbCi_ab,
                              TArray4& Gamma_iBjA_ab, TArray4& Gamma_iBJa_ba,
                              TArray4& Gamma_AbCd_ab, TArray4& Gamma_IjKl_ab,
                              TArray4& Gamma_IjAb_ab);

      // compute CCSD Xam (JCP, 103, 3561 (1995))
      void compute_Xam_ccsd(const TArray2& T1, const TArray4& T2,
                            const TArray2& L1, const TArray4& L2,
                            TArray2& Xam_tot, TArray2& Xiip);

      /** returns the 2-particle density matrix
      * @return \f$ \gamma^{pq}_{rs} \f$, respectively
      */
      TArray4 rdm2();

      /** computes spin-free V intermediate
        * @param symmetrize_p1_p2 if yes, make sure this holds for the result: V^{ij}_{mn} = V^{ji}_{nm}
        * @return \f$ V^{ij}_{mn} \f$, where ij refer to the geminal indices
      */
      TArray4 V_spinfree(bool symmetrize_p1_p2 = false);

      /** computes spin-free X intermediate
        * @param symmetrize_p1_p2 if yes, make sure this holds for the result: X^{ij}_{kl} = X^{ji}_{kl}
        * @return \f$ X^{ij}_{kl} \f$
      */
      TArray4 X_spinfree(bool symmetrize_p1_p2 = false);

      /** computes spin-free B intermediate
        * @param symmetrize_p1_p2 if yes, make sure this holds for the result: B^{ij}_{kl} = B^{ji}_{kl}
        * @return \f$ B^{ij}_{kl} \f$
      */
      TArray4 B_spinfree(bool symmetrize_p1_p2 = false);

      /**
       * provides T1 amplitude tensor
       * @param t1 act_occ by act_vir matrix
       */
      void set_T1(const RefSCMatrix& t1) { t1_ = t1; }
      /**
       * provides T1 CABS amplitude tensor
       * @param t1 act_occ by allvir (or CABS) matrix
       */
      void set_T1_cabs(const RefSCMatrix& t1_cabs) { t1_cabs_ = t1_cabs; }
      /**
       * provides Lambda1 amplitude tensor
       * @param l1 act_occ by act_vir matrix
       */
      void set_L1(const RefSCMatrix& l1) { l1_ = l1; }
      /**
       * provides T2 amplitudes
       * @param t2 array of T2 amplitudes, for AlphaBeta, AlphaAlpha, and (optionally) BetaBeta
       */
      void set_T2(const Ref<DistArray4> (&t2)[NSpinCases2]) { std::copy(t2, t2+NSpinCases2, t2_); }
      /**
       * provides Lambda2 amplitudes
       * @param l2 array of Lambda2 amplitudes, for AlphaBeta, AlphaAlpha, and (optionally) BetaBeta
       */
      void set_L2(const Ref<DistArray4> (&l2)[NSpinCases2]) { std::copy(l2, l2+NSpinCases2, l2_); }
      /**
       * provides (spin-free) RDM2
       * @param rdm2 a SpinFreeRDM<Two> object
       */
      void set_rdm2(const Ref<SpinFreeRDM<Two> >& rdm2) { rdm2_ = rdm2; }

      /// \sa _4()
      TArray4d& ijxy(const std::string& key);
      ///
      TArray22d& ij_xy(const std::string& key);
      /// \sa _2()
      TArray2& xy(const std::string& key);

      /// sieves x|o1|y -> x'|o1|y'
      /// does not throw only if each tile of the result depends only on 1 tile of the input
      TArray2& sieve(const TArray2& input, const std::string& output_annotation);

      /** Given a descriptive \c key, creates a rank-4 Array of integrals, or other related quantities
       *  The syntax of \c key is similar to that used by ParsedTwoBodyInt and TwoBodyMOIntsRuntime,
       *  but with te_type embedded into key.
        * The following te_types are understood:
        *   - g : \f$ r_{12}^{-1} \f$
        *   - r : \f$ f(r_{12}) \f$
        *   - gr : \f$ r_{12}^{-1} f(r_{12}) \f$
        *   - rTr : \f$ [f(r_{12}), [ \hat{T}_1 , f(r_{12})]] \f$
        *   - gamma : \f$ \Gamma \f$, 2-RDM, provided by \c set_rdm2()
        *   - T2 : 2-body excitation amplitudes, e.g. of the coupled-cluster method, provided by \c set_T2()
        *   - L2 : 2-body de-excitation amplitudes, e.g. of the coupled-cluster method, provided by \c set_L2()
        *
        *   Indices are interpreted using to_space().
        *
        * Example1: <i j|g|p a'> will create an array of <act_occ act_occ| r_{12}^{-1} |obs cabs> integrals
        * Example2: <i j_F(p')|g|p a'> will create an array of <act_occ cbs| r_{12}^{-1} |obs cabs> integrals
        * with the second index transformed using Fock matrix between act_occ and cbs spaces.
        *
        * Note that the indices are used to create the resulting tensor expression, hence these keys can be used
        * in composing expressions. For example, the full (non-diagonal) V intermediate without the CABS terms
        * can be computed as follows:
        * TArray4d V_ij_kl = _4("<i j|gr|k l>") - _4("<i j|g|p q>") * _4("<k l|r|p q>");
        *
        * \sa ijxy()
        */
      TA::expressions::TsrExpr<const TArray4d> _4(const std::string& key);

      /** Given a descriptive \c key, creates a rank-2 Array of integrals, or other related quantities
       *  The syntax of \c key is similar to that used by ParsedOneBodyInt,
       *  but with oe_type embedded into key.
        * The following oe_types are understood:
        *   - mu_i : electric dipole integral \f$ \bf{r}_i \f$ (i = x,y,z)
        *   - q_ij : electric quadrupole integral \f$ \bf{r}_i \bf{r}_j \f$ (ij = xx, xy, xz, yy, yz, zz)
        *   - gamma : \f$ \Gamma \f$, 1-RDM of the reference wave function
        *   - T1 : 1-body excitation amplitudes, e.g. of the coupled-cluster method, provided by \c set_T1() and/or \c set_T1_cabs()
        *   - L1 : 1-body de-excitation amplitudes, e.g. of the coupled-cluster method, provided by \c set_L1()
        *   - I : identity matrix
        *
        *   Indices are interpreted using to_space().
        *
        * Example1: <i|mu_x|p> will create an array of <act_occ| -x |obs> integrals
        * Example2: <j_F(p')|q_xy|a'> will create an array of <cbs| -xy |cabs> integrals
        * with the second index transformed using Fock matrix between act_occ and cbs spaces.
        *
        * Note that the indices are used to create the resulting tensor expression, hence these keys can be used
        * in composing expressions.
        *
        * \sa _4() \sa xy()
        */
      TA::expressions::TsrExpr<const TArray2> _2(const std::string& key);

      //TA::expressions::TensorExpression<TA::Tensor< TA::Tensor<T> > > _22(const std::string& key);

      /// like _4, produces geminal T tensor
      TA::expressions::TsrExpr<const TArray4Tg> _Tg(const std::string& key);

    private:
      madness::World& world_;
      Ref<R12WavefunctionWorld> r12world_;

      // extra data
      RefSCMatrix t1_;
      RefSCMatrix t1_cabs_;
      RefSCMatrix l1_;    // lambda1
      Ref<DistArray4> t2_[NSpinCases2];
      Ref<DistArray4> l2_[NSpinCases2]; // lambda2
      Ref<SpinFreeRDM<Two> > rdm2_;

      // utilities

      // deprecated
      /// computes a DistArray4 and converts into a TArray22
      std::shared_ptr<TArray22> ab_O_cd(const std::string& key,
                                        const int te_type);

      std::map<std::string, std::shared_ptr<TArray22> > tarray22_registry_;
      std::map<std::string, std::shared_ptr<TArray4d> > tarray4_registry_;
      std::map<std::string, std::shared_ptr<TArray2> > tarray2_registry_;
      std::map<std::string, std::shared_ptr<TArray22d> > tarray22d_registry_;
      std::map<std::string, std::shared_ptr<TArray4Tg> > tarray4tg_registry_;

      // helps to create TArray4Tg
      static expressions::TGeminalGenerator<T> tg_s0_gen;
      static expressions::TGeminalGenerator<T> tg_s1_gen;

      // deprecated
      TArray22& _(const std::string& key);

      /**
       *  Converts an index label to the label of the corresponding space
       *  Indices can be composite, e.g. also include a 1-particle operator transform
       *  Base indices are transformed using to_space_()
       *
       * converts index label to space label using the canonical dictionary documented in to_space_()
       * @param index label
       * @return space label
       */
      std::string to_space(const std::string& index);

      /**   Converts an index label to the label of the corresponding space.
       *    Composite indices are not allowed.
        *   Here's the key index dictionary that can be used (\sa to_space):
        *     - i,j,k,l -> i (active occupied, act_occ)
        *     - m,n -> m (occupied, occ)
        *     - a,b,c,d -> a (active virtual, act_vir)
        *     - e,f -> e (virtual, vir)
        *     - p,q,r,s -> p (orbital basis, obs)
        *     - p',q',r',s' -> p' (complete basis, cbs)
        *     - a', b', c', d' -> a' (cabs = complete basis - orbital basis)
        *     - A', B', C', D' -> A' (complete virtuals = a + a', cvir)
        *   Indices can have arbitrary numerical labels. Examples:
        *   to_space_("s'") -> "p'"
        *   to_space_("j17") -> "i"
        *
        * @param index
        * @return space label
        */
      static std::string to_space_(std::string index);

      /// returns the hashmarks for tiling a space corresponding to space_label.
      /// space_label should have been canonicalized using to_space().
      /// By default, occupied space chopped finely, and other spaces are chopped coarsely.
      std::vector<size_t> space_hashmarks(std::string space_label) const;

      template <size_t NDIM>
      TA::TiledRange
      make_trange(const std::array<std::string, NDIM>& spaces) const;

      /** returns the 1-particle reduced density matrix in spaces row and col. Will throw if the RDM of r12world()->refwfn()
       * is expressed in a space not supported by the same basis as \c row and \c col
      * @return \f$ \gamma^{r}_{c} \f$, respectively
      */
      RefSCMatrix rdm1(const Ref<OrbitalSpace>& row,
                       const Ref<OrbitalSpace>& col) const;

      typedef enum {
        ij=0, ji=1
      } ij_type;
      /// takes ij|o|ij or ij|o|ji
      template <typename Array22>
      TA::expressions::TsrExpr<const Array22> take(const Array22& ij_o_pq,
                                                             ij_type IJ) {
        //typedef typename Array22::value_type value_type;
          typedef TA::Tensor<TA::Tensor<T> > value_type;
            if (IJ == ij) {
              sc::expressions::diag_tensor2_op<value_type,false> diag_op;
              return make_unary_tensor(ij_o_pq("i,j"),
                                       diag_op);
            }
            // else: IJ == ji
            sc::expressions::diag_tensor2_op<value_type,true> diag_op;
            return make_unary_tensor(ij_o_pq("i,j"),
                                     diag_op);
          }

      /// computes ij|o1|pq . ij|o2|pq
      template <typename Array22>
      TA::expressions::TsrExpr<const TArray2> dotket(const Array22& ij_o1_pq,
                                                               const Array22& ij_o2_pq,
                                                               bool transpose_o2_ket = false) {
        //typedef typename Array22::value_type value_type;
        typedef TA::Tensor<TA::Tensor<T> > value_type;
        if (transpose_o2_ket == false) {
          sc::expressions::trace_tensor2_op<value_type,false> trace_op;
          return make_binary_tensor(ij_o1_pq("i,j"), ij_o2_pq("i,j"),
                                    trace_op);
        }
        // else: transpose_o2_ket = true
        sc::expressions::trace_tensor2_op<value_type,false> trace_op;
        return make_binary_tensor(ij_o1_pq("i,j"), ij_o2_pq("j,i"),
                                  trace_op);
      }

  };

}; // end of namespace sc

#include <chemistry/qc/mbptr12/sr_r12intermediates_util.h>
#include <chemistry/qc/mbptr12/sr_r12intermediates_VXB_diag.h>

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
