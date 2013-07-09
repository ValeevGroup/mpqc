//
// sr_r12intermediates_VXB_diag.h
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesVXBdiag_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesVXBdiag_h

#include <TiledArray/eigen.h>
#include <TiledArray/algebra/conjgrad.h>

namespace sc {

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::V_diag() {

    TArray2 V_ij_ij_cabs = dotket(ij_xy("<i j|g|m a'>"), ij_xy("<i j|r|m a'>"));
    // don't need this in closed shell!
    //TArray2 V_ij_ij_cabs1 = dotket(ij_xy("<i j|g|a' m>"), ij_xy("<i j|r|a' m>"));

    TArray2 V_ij_ij = take(ij_xy("<i j|gr|p q>"), ij) - dotket(ij_xy("<i j|g|p q>"), ij_xy("<i j|r|p q>"))
                      - V_ij_ij_cabs("i,j") - V_ij_ij_cabs("j,i");

    TArray2 V_ij_ji_cabs = dotket(ij_xy("<i j|g|m a'>"), ij_xy("<i j|r|m a'>"), true);
    TArray2 V_ij_ji = take(ij_xy("<i j|gr|p q>"), ji) - dotket(ij_xy("<i j|g|p q>"), ij_xy("<i j|r|p q>"), true)
                      - V_ij_ji_cabs("i,j") - V_ij_ji_cabs("j,i");

    return std::make_pair(V_ij_ij,V_ij_ji);

  }

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::X_diag() {

    TArray2 X_ij_ij_cabs = dotket(ij_xy("<i j|r|m a'>"), ij_xy("<i j|r|m a'>"));
    // don't need this in closed shell!
    //TArray2 X_ij_ij_cabs1 = dotket(ij_xy("<i j|r|a' m>"), ij_xy("<i j|r|a' m>"));

    TArray2 X_ij_ij = take(ij_xy("<i j|r2|p q>"), ij) - dotket(ij_xy("<i j|r|p q>"), ij_xy("<i j|r|p q>"))
                      - X_ij_ij_cabs("i,j") - X_ij_ij_cabs("j,i");

    TArray2 X_ij_ji_cabs = dotket(ij_xy("<i j|r|m a'>"), ij_xy("<i j|r|m a'>"), true);
    TArray2 X_ij_ji = take(ij_xy("<i j|r2|p q>"), ji) - dotket(ij_xy("<i j|r|p q>"), ij_xy("<i j|r|p q>"), true)
                      - X_ij_ji_cabs("i,j") - X_ij_ji_cabs("j,i");

    return std::make_pair(X_ij_ij,X_ij_ji);

  }

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::B_diag() {

    /// this is incomplete at the moment
    TArray2 B_ij_ij = take(ij_xy("<i j|rTr|p q>"), ij) + take(ij_xy("<i_hJ(p') j|r2|p q>"), ij)
                      - dotket(ij_xy("<i j|r|p' q'>"), ij_xy("<i j|r|p' q'_K(r')>"));

    TArray2 B_ij_ji = take(ij_xy("<i j|rTr|p q>"), ji) + take(ij_xy("<i_hJ(p') j|r2|p q>"), ji);

    return std::make_pair(B_ij_ij,B_ij_ji);

  }

  namespace detail {
    /** this functor helps to implement conjugate gradient CABS singles solver
     */
    template<typename T>
    struct _CABS_singles_h0t1 {

        typedef TiledArray::Array<T, 2> Array;

        /**
         * @param h0_AB allvirt/allvirt Fock operator
         * @param h0_ij occ/occ Fock operator
         */
        _CABS_singles_h0t1(const Array& h0_AB, const Array& h0_ij) :
            H0_AB(h0_AB), H0_IJ(h0_ij) {
        }

        const Array& H0_AB;
        const Array& H0_IJ;

        /**
         * @param[in] T1 t_i^A
         * @param[out] R1 R_i^A
         */
        void operator()(const Array& T1, Array& R1) {
          R1 = T1("i,b") * H0_AB("b,a") - H0_IJ("i,j") * T1("j,a");
        }
    };

    /** this functor helps to implement orbital response
     */
    template<typename T>
    struct _OrbResponse {

        typedef TiledArray::Array<T, 2> Array2;
        typedef TiledArray::Array<T, 4> Array4;

        /**
         * @param f_AB virt/virt Fock operator
         * @param f_ij occ/occ Fock operator
         * @param g_ij_ab <ij|ab>
         * @param g_ia_jb <ia|jb>
         */
        _OrbResponse(const Array2& f_AB, const Array2& f_ij,
                     const Array4& g_ij_ab, const Array4& g_ia_jb) :
            F_AB(f_AB), F_IJ(f_ij), G_IJ_AB(g_ij_ab), G_IA_JB(g_ia_jb) {
        }

        const Array2& F_AB;
        const Array2& F_IJ;
        const Array4& G_IJ_AB;
        const Array4& G_IA_JB;

        /**
         * @param[in] kappa kappa_i^a
         * @param[out] residual residual_i^a
         */
        void operator()(const Array2& kappa, Array2& residual) {
          residual =  kappa("i,b") * F_AB("b,a") - F_IJ("i,j") * kappa("j,a")
          + 4.0 * G_IJ_AB("i,j,a,b") *  kappa("j,b")
                   -       G_IJ_AB("i,j,b,a") *  kappa("j,b")
                   -       G_IA_JB("i,a,j,b") *  kappa("j,b");
        }
    };

    /// makes a diagonal 2-index preconditioner: pc_x^y = -1/ ( <x|O1|x> - <y|O2|y> )
    template <typename T>
    struct diag_precond2 {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        diag_precond2(const EigenMatrixX& O1_mat,
                      const EigenMatrixX& O2_mat) :
                          O1_mat_(O1_mat), O2_mat_(O2_mat) {
        }
        template <typename Index> T operator()(const Index& i) {
          return 1.0 / (- O1_mat_(i[0], i[0]) + O2_mat_(i[1], i[1]));
        }

      private:
        EigenMatrixX O1_mat_;
        EigenMatrixX O2_mat_;
    };

    /// makes a diagonal 4-index preconditioner: pc_xy^zw = -1/ ( <x|O1|x> + <y|O2|y> - <z|O3|z> - <w|O4|w> )
    template <typename T>
    struct diag_precond4 {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        diag_precond4(const EigenMatrixX& O1_mat,
                      const EigenMatrixX& O2_mat,
                      const EigenMatrixX& O3_mat,
                      const EigenMatrixX& O4_mat) :
                          O1_mat_(O1_mat), O2_mat_(O2_mat),
                          O3_mat_(O3_mat), O4_mat_(O4_mat) {
        }
        template <typename Index> T operator()(const Index& i) {
          return 1.0 / (- O1_mat_(i[0], i[0]) - O2_mat_(i[1], i[1]) + O3_mat_(i[2], i[2]) + O4_mat_(i[3], i[3]));
        }

      private:
        EigenMatrixX O1_mat_;
        EigenMatrixX O2_mat_;
        EigenMatrixX O3_mat_;
        EigenMatrixX O4_mat_;
    };

  } // namespace sc::detail

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::rdm1() {

    if (0) {
      {
      typedef TiledArray::Array<T,2> Array;
      Array FiA = _2("<i|F|A'>");
      Array FAi = _2("<A'|F|i>");
      std::cout << "<i|F|A'>:" << std::endl << FiA << std::endl;
      std::cout << "<A'|F|i>:" << std::endl << FAi << std::endl;
      Array FiA_2(FiA.get_world(), FiA.trange());
      FiA_2("i,A'") = FAi("A',i");
      std::cout << "<i|F|A'>=Perm(<A'|F|i>):" << std::endl << FiA_2 << std::endl;
      }
      {
      typedef TiledArray::Array<T,4> Array;
      Array g_ij_ab = _4("<i j|g|a b>");
      Array g_ab_ij = _4("<a b|g|i j>");
      std::cout << "<i j|g|a b>:" << std::endl << g_ij_ab << std::endl;
      std::cout << "<a b|g|i j>:" << std::endl << g_ab_ij << std::endl;
      Array g_ij_ab_2(g_ij_ab.get_world(), g_ij_ab.trange());
      g_ij_ab_2("i,j,a,b") = g_ab_ij("a,b,i,j");
      std::cout << "<i j|g|a b>=Perm(<a b|g|i j>):" << std::endl << g_ij_ab_2 << std::endl;
      Array should_be_zero = g_ij_ab("i,j,a,b") - g_ab_ij("a,b,i,j");
      std::cout << "<i j|g|a b> - Perm(<a b|g|i j>):" << std::endl << should_be_zero << std::endl;
      const double max_nonzero = norminf(should_be_zero("i,j,a,b"));
      std::cout << "|| <i j|g|a b> - Perm(<a b|g|i j>) ||_\infty = " << max_nonzero << std::endl;
      }
      {
      typedef TiledArray::Array<T,2> Array;
      Array mu_z_ij = _2("<i|mu_z|j>");
      Array gamma_ij = _2("<i|gamma|j>");
      const double mu_z_e = dot(mu_z_ij("i,j"), gamma_ij("i,j"));
      double mu_z_n = 0.0;
      Ref<Molecule> mol = r12world_->basis()->molecule();
      for(int a=0; a<mol->natom(); ++a) {
        mu_z_n += mol->Z(a) * mol->r(a, 2);
      }
      std::cout << "mu_z = " << -mu_z_e+mu_z_n << std::endl;
      }
    }

    TArray2 Fab = _2("<a|F|b>");
    TArray2 Iij = _2("<i|I|j>");
    //std::cout << "Fock(vir,vir)\n" << Fab << std::endl;
    //std::cout << "Idenity(occ,occ)\n" << Iij << std::endl;

    // can only ask for T1 with i in bra!
    // since we computed T1 CABS, they are expressed in terms of all virtuals = A'
    // if you turn off vir-CABS coupling, use a' (i.e. CABS only)
    TArray2 T1iA = _2("<i|T1|A'>");
    //t1_cabs_.print("T1(RefSCMatrix)");
    //std::cout << "T1(cabs)\n" << T1iA << std::endl;
    TArray2 T1ia = _2("<i|T1|a>");
    //std::cout << "T1(cabs) => i by a block\n" << T1ia << std::endl;

    // recompute E2(CABS) = T1_cabs . H1
    const double E2_cabs = 2.0 * dot(T1iA("i,A'"), _2("<i|F|A'>"));
    std::cout << "E2_cabs (recomputed) = " << E2_cabs << std::endl;

    // recompute T1_cabs and re-recompute E2_cabs
    {
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<A'|F|B'>");
      // this computes Z_i^A' = T_i^B' F_B'^A' - F_i^j T_j^A'
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|A'>");
      Array minus_FiA = -1.0 * FiA("i,A'");
      Array T1_recomp = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <A'|F|A'>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,A'");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_recomp,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      const double E2_cabs = 2.0 * dot(T1_recomp("i,A'"), _2("<i|F|A'>")); // 2 accounts for spin
      std::cout << "E2_cabs (re-recomputed) = " << E2_cabs << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|A'>"), T1_recomp("i,A'"))
                            - dot(_2("<i|mu_z|j>"), T1_recomp("i,A'") * T1_recomp("j,A'") )
                            + dot(_2("<A'|mu_z|B'>"), T1_recomp("i,A'") * T1_recomp("i,B'") );
        std::cout << "Mu_z (2)_S = " << 2*mu_z_e << std::endl; // 2 accounts for spin degeneracy
        std::cout << "Mu_z ref = " << (dot(_2("<i|mu_z|j>"),_2("<i|gamma|j>")))
                  << std::endl; // 2 is included in gamma!
      }

      // compute orbital rotation multipliers in the (2)_S Lagrangian
      if (1) {
        TArray2 Tia = T1_recomp("j,A'") * _2("<A'|I|a>");
        TArray2 Xia_1 = 2* (_2("<i|F|j>") * Tia("j,a") - T1_recomp("i,B'") * _2("<B'|F|a>"))
                        + 2 * _2("<i|F|C'>") * T1_recomp("j,C'") * Tia("j,a")
                        - 2 * (2 * _4("<j i|g|B' a>") - _4("<j i|g|a B'>") + 2 * _4("<j a|g|B' i>") - _4("<j a|g|i B'>")) * T1_recomp("j,B'")
                        - 2 * (2 * _4("<B' i|g|C' a>") - _4("<B' i|g|a C'>")) * T1_recomp("j,B'") * T1_recomp("j,C'")
                        + 2 * (2 * _4("<j i|g|k a>") - _4("<j i|g|a k>")) * T1_recomp("j,B'") * T1_recomp("k,B'");
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;

        // verify the solution:
        {
          TArray2 res = Xia_1;
          response_lhs_eval(kappa, res);
          std::cout << "should be zero = " << TA::expressions::norm2(res("i,a") - Xia_1("i,a"));
        }

        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << "mu_z_e (orb response) = " << 2*mu_z_e << std::endl;

      }
    }

    // now recompute T1_cabs and E2_cabs using only CABS for perturbation!
    {
      std::cout << "computing E2_cabs due to CABS only as the first-order space\n";
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<a'|F|b'>");
      // this computes Z_i^a' = T_i^b' F_b'^a' - F_i^j T_j^a'
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|a'>");
      Array minus_FiA = -1.0 * FiA("i,a'");
      Array T1_recomp = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <a'|F|a'>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,a'");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_recomp,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      const double E2_cabs = 2.0 * dot(T1_recomp("i,a'"), _2("<i|F|a'>")); // 2 accounts for spin
      std::cout << "E2_cabs (re-recomputed) = " << E2_cabs << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|a'>"), T1_recomp("i,a'"))
                            - dot(_2("<i|mu_z|j>"), T1_recomp("i,a'") * T1_recomp("j,a'") )
                            + dot(_2("<a'|mu_z|b'>"), T1_recomp("i,a'") * T1_recomp("i,b'") );
        std::cout << "Mu_z (2)_S = " << 2*mu_z_e << std::endl; // 2 accounts for spin degeneracy
      }

      // compute orbital rotation multipliers in the (2)_S Lagrangian
      if (1) {
        // only include the first-order terms in the Lagrangian derivative
        TArray2 Xia_1 = - 2 * T1_recomp("i,b'") * _2("<b'|F|a>")
                      - 2 * (2 * _4("<j i|g|b' a>") - _4("<j i|g|a b'>") + 2 * _4("<j a|g|b' i>") - _4("<j a|g|i b'>")) * T1_recomp("j,b'")
                      - 2 * (2 * _4("<b' i|g|c' a>") - _4("<b' i|g|a c'>")) * T1_recomp("j,b'") * T1_recomp("j,c'")
                      + 2 * (2 * _4("<j i|g|k a>") - _4("<j i|g|a k>")) * T1_recomp("j,b'") * T1_recomp("k,b'");
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;
        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << "mu_z_e (orb reponse) = " << 2*mu_z_e << std::endl;

      }
    }

    // compute T1_obs and E2_obs
    {
      std::cout << "computing E2_obs (virtuals is the first-order space)\n";
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<a|F|b>");
      // this computes Z_i^a = T_i^b F_b^a - F_i^j T_j^a
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|a>");
      Array minus_FiA = -1.0 * FiA("i,a");
      Array T1_obs = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <a|F|a>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,a");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_obs,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      // same as the lagrangian since solved for T1 exactly
      const double E2_obs = 2.0 * dot(T1_obs("i,a"), _2("<i|F|a>")); // 2 accounts for spin
      std::cout << "E2_obs = " << E2_obs << std::endl;
      const double E0_obs = dot(_2("<i|F|j>"), _2("<i|gamma|j>")); // 2 included in gamma
      std::cout << "E0_obs = " << scprintf("%15.10lf",E0_obs) << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|a>"), T1_obs("i,a"))
                            - dot(_2("<i|mu_z|j>"), T1_obs("i,a") * T1_obs("j,a") )
                            + dot(_2("<a|mu_z|b>"), T1_obs("i,a") * T1_obs("i,b") );
        std::cout << scprintf("Mu_z (2)_S OBS = %25.15lf", 2*mu_z_e) << std::endl; // 2 accounts for spin degeneracy
        std::cout << scprintf("Mu_z ref = %25.15lf", (dot(_2("<i|mu_z|j>"),_2("<i|gamma|j>")))) << std::endl;
      }

      // compute orbital rotation multipliers in the (2)_S OBS Lagrangian
      if (1) {
        // only include the first-order terms in the Lagrangian derivative
        TArray2 Xia_1 = _2("<i|F|j>") * T1_obs("j,a") - T1_obs("i,b") * _2("<b|F|a>")
                        + 2 * _2("<i|F|c>") * T1_obs("j,c") * T1_obs("j,a")
                        + 2 * T1_obs("i,c") * T1_obs("j,c") * _2("<j|F|a>")
                        ;
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;
        TArray2 Fia = _2("<i|F|a>");
        std::cout << "Fia, should be close to Xia_1:\n" << Fia << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;
        std::cout << scprintf("E_2 (orb reponse) = %25.15lf", 2*dot(kappa("i,a"), _2("<i|F|a>"))) << std::endl;
        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << scprintf("mu_z_e (orb reponse) = %25.15lf", 2*mu_z_e) << std::endl;

      }
    }

    if (1) {
      for(int i=0; i<10; ++i) {
        TArray4 should_be_zero = _4("<b j|g|a i>") - _4("<b i|g|a j>");
        std::cout << "should be 0: " << TA::expressions::norminf(should_be_zero("b,j,a,i")) << std::endl;
      }
    }

    // this now works ... thanks Justus!
#if 1
    //TArray4 A = _2("<i|I|j>") * _2("<a|F|b>") - _4("<i j|g|a b>");
    //std::cout << "A\n" << A << std::endl;
#endif

    /// this is just an example of how to compute the density
    TArray2 r2_i_j = _4("<i j|r|p q>") * _4("<k_F(p) j|r|p q>");
    //std::cout << "<ij|r|pq> . <kj|r|pq>\n" << r2_i_j << std::endl;

    // this is another random contraction, useful for non-diagonal X intermediate
    //TArray4 x = _4("<i j|r|p q>") * _4("<k l|r|p q>");

    return r2_i_j;
  }


  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::V_spinfree(bool symmetrize_p1_p2) {

    TArray4 V_ij_mn = _4("<i1 i2|gr|m1 m2>") - _4("<i1 i2|r|p1 p2>") * _4("<m1 m2|g|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<m1 m2|g|m3 a'>");

    if (symmetrize_p1_p2)
      return 0.5 * (V_ij_mn("i1,i2,m1,m2") + V_ij_mn("i2,i1,m2,m1"));
    else
      return V_ij_mn;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::X_spinfree(bool symmetrize_p1_p2) {

    TArray4 X_ij_kl = _4("<i1 i2|r2|j1 j2>") - _4("<i1 i2|r|p1 p2>") * _4("<j1 j2|r|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<j1 j2|r|m3 a'>");

    if (symmetrize_p1_p2)
      return 0.5 * (X_ij_kl("i1,i2,j1,j2") + X_ij_kl("i2,i1,j2,j1"));
    else
      return X_ij_kl;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::B_spinfree(bool symmetrize_p1_p2) {

    TArray4 B_ij_kl =

        // everything seems scaled up by factor of 2 relative to Eq.(12) in J. Chem. Phys. 135, 214105 (2011),
        // due to including particle 1 and particle 2 contributions?

        // diag                      Q
        _4("<i1 i2|rTr|j1 j2>") + 2.0 * _4("<i1 i2|r2|j1 j2_hJ(p')>")

        //           rKr
        -2.0 * _4("<i1 i2|r|r' s'>") * _4("<j1 j2|r|r' s'_K(p')>")

        //           rFr
        -2.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(p)>")

        //           rFr_2, extra 2 due to bra-ket symmetrization
        -4.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(a')>")

        //           rFGr
        - _4("<i1 i2|r|n_gamma(m) b'>") * _4("<j1 j2|r|n b'_F(a')>")

        //           rFGr_2
        - _4("<i1 i2|r|n_gamma(m) a'>") * _4("<j1 j2|r|n_F(p') a'>");

    B_ij_kl = 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));
    B_ij_kl = 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("j1,j2,i1,i2"));

    if (symmetrize_p1_p2)
      return 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));
    else
      return B_ij_kl;
  }


}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
