/*
 * ccsd_pno_impl.h
 *
 *  Created on: Jan 5, 2017
 *      Author: jinmei
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_

namespace mpqc {

  namespace lcao {

    namespace detail {
      // compute MP2 T2 amplitudes
      template <typename Tile, typename Policy>
      TA::DistArray<Tile, Policy> compute_mp2_t2(lcao::LCAOFactory<Tile, Policy> &lcao_factory,
                                                 const std::shared_ptr<Eigen::VectorXd> &orbital_energy,
                                                 const std::shared_ptr<mpqc::TRange1Engine> &trange1_engine,
                                                 bool df) {
        TA::DistArray<Tile, Policy> g_abij;
        ExEnv::out0() << indent << "Using density fitting: " << df << std::endl;
        g_abij = lcao_factory.compute(df ? L"<a b|G|i j>[df]" : L"<a b|G|i j>");

        auto n_occ = trange1_engine->get_occ();
        auto n_frozen = trange1_engine->get_nfrozen();

        // compute t2 amplitudes
        TA::DistArray<Tile, Policy> t2_abij = d_abij(g_abij, *orbital_energy, n_occ, n_frozen);

        // test mp2 energy:
        double mp2 =  TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), t2_abij("a,b,i,j"));
        ExEnv::out0() << "Test: MP2 Energy is " << mp2 << std::endl;

        return t2_abij;
      }

      // print out details of CCSD iterations
      inline void print_ccsd(int iter, double dE,
                             double error, double error_r1, double error_r2,
                             double E1, double time) {
        if (iter == 0) {
          std::printf("%3s \t %10s \t %10s \t %12s \t %12s \t %15s \t %10s \n", "iter", "deltaE",
                      "residual", "norm (r1)", "norm (r2)", "energy", "total time/s");
        }
        std::printf("%3i \t %10.5e \t %10.5e \t %10.5e \t %10.5e \t %15.12f \t %10.1f \n",
                    iter, dE, error, error_r1, error_r2, E1, time);
      }
    } // namespace detail

    /**
     * KeyVal constructor
     * @param kv
     *
     * keywords : all keywords for LCAOWavefunciton
     *
     * | KeyWord | Type | Default| Description |
     * |---------|------|--------|-------------|
     * | ref     | wfn  | none   | reference wavefunction, RHF for example |
     * | df      | bool | false  | choice of using density fitting
     */
    template<typename Tile, typename Policy>
    CCSD_PNO<Tile, Policy>::CCSD_PNO(const KeyVal &kv): CCSD<Tile, Policy>(kv) {
      tcut_ = kv.value<double>("tcut", 1e-6);
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::compute(PropertyBase *pb)  {
      throw std::runtime_error("Not Implemented!!");
    }

    // compute MP2 T2 amplitudes
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::compute_mp2_t2() {

      auto &lcao_factory = this->lcao_factory();
      t2_mp2_ = detail::compute_mp2_t2(lcao_factory, this->orbital_energy(),
                                       this->trange1_engine(), this->df_);
    }

    // compute converting matrices for reblocking MP2 T2
    // occ_convert(old_i, new_i), vir_convert(old_a, new_a)
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::compute_M_reblock(TA::DistArray<Tile, Policy> &occ_convert,
                                                   TA::DistArray<Tile, Policy> &vir_convert) {

      auto &lcao_factory = this->lcao_factory();
      auto &world = lcao_factory.world();

      std::size_t occ = this->trange1_engine()->get_occ();
      std::size_t vir = this->trange1_engine()->get_vir();
      std::size_t all = this->trange1_engine()->get_all();
      std::size_t n_frozen = this->trange1_engine()->get_nfrozen();

      std::size_t occ_blk_size = 1;
      std::size_t vir_blk_size = vir;

      TA::TiledRange1 old_occ = this->trange1_engine()->get_active_occ_tr1();
      TA::TiledRange1 old_vir = this->trange1_engine()->get_vir_tr1();

      auto new_tr1 =
          std::make_shared<TRange1Engine>(occ, all, occ_blk_size, vir_blk_size, n_frozen);

      TA::TiledRange1 new_occ = new_tr1->get_active_occ_tr1();
      TA::TiledRange1 new_vir = new_tr1->get_vir_tr1();

      mpqc::detail::parallel_print_range_info(world, new_occ, "PNO new Occ");
      mpqc::detail::parallel_print_range_info(world, new_vir, "PNO new Vir");

      occ_convert = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
                       world, old_occ, new_occ, 1.0);

      vir_convert = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
                       world, old_vir, new_vir, 1.0);
    }

    // decompose MP2 T2 amplitudes
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::decom_t2_mp2() {

      double threshold = tcut_;
      ExEnv::out0() << "tcut is " << tcut_ << std::endl;

      integer vir = this->trange1_engine()->get_vir();

      std::size_t orb_i = 1, orb_j = 1;

      // decompose and reconstruct T2
      auto transfrom = [&](Tile &in_tile) -> float {

        // copy in_tile which is ab matrix of T^ij
        Tile tab_tile = in_tile.clone();

        // get the dimensions of input tile
        const auto extent = in_tile.range().extent();
        const auto a = extent[0];
        const auto b = extent[1];

        std::stringstream ss;
        ss << "input tile: " << in_tile << std::endl;

#define PNO_DECOM 1
#define T2_SVD 0

#if PNO_DECOM
        auto Tab = TA::eigen_map(tab_tile, a, b);
        // compute virtual or PNO density
        TA::EigenMatrixXd Dab = Tab * (4.0 * Tab - 2.0 * Tab.transpose()).transpose()
                              + Tab.transpose() * (4.0 * Tab - 2.0 * Tab.transpose());
        // when i=j, Dab = Tab*Tab' + Tab'*Tab
        if ((Tab - Tab.transpose()).norm() < 1e-10) {
          ss << "i = j, (Tab - Tab.transpose()).norm(): "
             << (Tab - Tab.transpose()).norm() << std::endl;
          Dab = Dab/2.0;
        }

        // compute eigenvalues ans eigenvectors of Dab
        Eigen::SelfAdjointEigenSolver<TA::EigenMatrixXd> es(Dab);
        TA::EigenVectorXd S_es = es.eigenvalues();
        TA::EigenMatrixXd C_es = es.eigenvectors();

        // truncate eigenvectors with corresponding eigenvalues smaller than threshold
        int pos_trunc = a - 1;
        for(; pos_trunc >= 0; --pos_trunc) {
          if(std::abs(S_es(pos_trunc)) < threshold)
            break;
        }
        int rank = a - pos_trunc - 1;
        if (rank > 0) {
          auto tab_pno = TA::eigen_map(in_tile, a, b);
          tab_pno.noalias() = C_es.rightCols(rank) * C_es.rightCols(rank).transpose()
                            * Tab
                            * C_es.rightCols(rank) * C_es.rightCols(rank).transpose();
        } else {
            std::fill(in_tile.begin(), in_tile.end(), 0.0);
        }

        ss << "Dab: " << std::endl << Dab << std::endl;
        ss << "S_es: " << std::endl << S_es << std::endl;
        ss << "C_es: " << std::endl << C_es << std::endl;
        ss << "pair (" << orb_i << "," << orb_j << ")  rank: " << rank << std::endl;
        ss << "S_es (trunc): " << std::endl << S_es.tail(rank) << std::endl;
        ss << "C_es (trunc): " << std::endl << C_es.rightCols(rank) << std::endl;
#endif // PNO_DECOM

#if T2_SVD
        const auto x = std::min(a,b);

        // allocate memory for SVD
        std::unique_ptr<double[]> u_ax(new double[a * x]);
        std::unique_ptr<double[]> v_xb(new double[b * x]);
        std::unique_ptr<double[]> s(new double[x]);

        // SVD of tile: (t_ab)^T
        tensor::algebra::svd(in_tile.data(), s.get(), u_ax.get(), v_xb.get(), a, b, 'A');
//        // print out results for test purpose
//        ss << "s: " << std::endl;
//        for (std::size_t i = 0; i < a; ++i) {
//          ss << s.get()[i] << "  ";
//        }
//        ss << std::endl;
//        ss << "u: " << std::endl;
//        for (std::size_t i = 0; i < a; ++i) {
//          for (int j = 0; j < b; ++j) {
//            ss << u_ax.get()[i+x*j] << "  ";
//          }
//          ss << std::endl;
//        }
//        ss << std::endl;
//        ss << "v: " << std::endl;
//        for (std::size_t i = 0; i < a; ++i) {
//          for (int j = 0; j < b; ++j) {
//            ss << v_xb.get()[i+x*j] << "  ";
//          }
//          ss << std::endl;
//        }
//        ss << std::endl;

        std::size_t rank = 0;
        for(; rank < x; ++rank) {
          // check singular value for truncation threshold
          if(std::abs(s.get()[rank]) < threshold)
            break;

          // fold s into u_ab:
          // u_{a} = u_{a} * s_{a}
          madness::cblas::scal(a, s.get()[rank], u_ax.get() + (a * rank), 1);
        }

        // compute output tile: (t_ab)^T = u_a{c} * v_{c}b
        madness::cblas::gemm(madness::cblas::NoTrans, madness::cblas::NoTrans,
            a, b, rank, 1.0, u_ax.get(), a, v_xb.get(), b, 0.0, in_tile.data(), a);

        ss << "pair (" << orb_i << "," << orb_j << ")  rank: " << rank << std::endl;
#endif // T2_SVD

        ++orb_i;
        ++orb_j;

        ss << "output tile: " << in_tile << std::endl;

        // compute the norm of difference between input and output tile
        ss << "(t_new - t_old).norm(): " << (TA::subt(in_tile, tab_tile)).norm()
           << std::endl;

        std::printf("%s", ss.str().c_str());

        return in_tile.norm();
      };

      TA::foreach_inplace(t2_mp2_, transfrom);
    }

    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::compute_ccsdpno_df(TA::DistArray<Tile, Policy> &t1,
                                                      TA::DistArray<Tile, Policy> &t2) {

      using TArray = TA::DistArray<Tile, Policy>;

      auto &world = this->wfn_world()->world();
      bool accurate_time = this->lcao_factory().accurate_time();

      auto n_occ = this->trange1_engine()->get_occ();
      auto n_frozen = this->trange1_engine()->get_nfrozen();

      if (world.rank() == 0) {
          std::cout << "Use DF CCSD Compute" << std::endl;
      }

      auto tmp_time0 = mpqc::now(world, accurate_time);
      // get all two electron integrals
      TArray g_abij = this->get_abij();
      TArray g_ijkl = this->get_ijkl();
      TArray g_abcd = this->get_abcd();
      TArray X_ai = this->get_Xai();
      TArray g_iajb = this->get_iajb();
      TArray g_iabc = this->get_iabc();
      TArray g_aibc = this->get_aibc();
      TArray g_ijak = this->get_ijak();
      TArray g_ijka = this->get_ijka();
      auto tmp_time1 = mpqc::now(world, accurate_time);
      auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (this->print_detail_) {
          mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                                   "\n");
      }

      TArray f_ai = this->get_fock_ai();

      // store d1 to local
      TArray d1 = create_d_ai<Tile,Policy>(f_ai.world(), f_ai.trange(),
                                           *this->orbital_energy(), n_occ, n_frozen);

      // initial value for t1
      t1("a,i") = f_ai("a,i") * d1("a,i");
      t1.truncate();

      // intitial value for t2 is PNO constructed t2

      TArray tau;
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      double E0 = 0.0;
      double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                  TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
      double mp2 = E1;
      double dE = std::abs(E1 - E0);

      mpqc::utility::print_par(world, "MP2 Energy      ", mp2, "\n");

      // optimize t1 and t2
      std::size_t iter = 0ul;
      double error = 1.0;
      TArray r1;
      TArray r2;

//      bool less = this->kv_.value<bool>("less_memory", false);
      bool less = false;

      if (world.rank() == 0) {
          std::cout << "Start Iteration" << std::endl;
          std::cout << "Max Iteration: " << this->max_iter_ << std::endl;
          std::cout << "Convergence: " << this->converge_ << std::endl;
          std::cout << "AccurateTime: " << accurate_time << std::endl;
          std::cout << "PrintDetail: " << this->print_detail_ << std::endl;
          if (less) {
              std::cout << "Less Memory Approach: Yes" << std::endl;
          } else {
              std::cout << "Less Memory Approach: No" << std::endl;
          }
      }

      auto diis = this->get_diis(world);

      while (iter < this->max_iter_) {
          // start timer
          auto time0 = mpqc::fenced_now(world);
          TArray::wait_for_lazy_cleanup(world);

          auto t1_time0 = mpqc::now(world, accurate_time);
          TArray h_ki, h_ac;
          {
              // intermediates for t1
              // external index i and a
              // vir index a b c d
              // occ index i j k l
              TArray h_kc;

              // compute residual r1(n) = t1(n+1) - t1(n)
              // external index i and a
              tmp_time0 = mpqc::now(world, accurate_time);
              r1("a,i") = f_ai("a,i") - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k");

              {
                  h_ac("a,c") =
                  -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");
                  r1("a,i") += h_ac("a,c") * t1("c,i");
              }

              {
                  h_ki("k,i") =
                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");
                  r1("a,i") -= t1("a,k") * h_ki("k,i");
              }

              {
                  h_kc("k,c") =
                  f_ai("c,k") +
                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
                  r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                              t1("c,i") * t1("a,k"));
              }

              tmp_time1 = mpqc::now(world, accurate_time);
              tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
              if (this->print_detail_) {
                  mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
              }

              tmp_time0 = mpqc::now(world, accurate_time);
              r1("a,i") += (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k");

              r1("a,i") +=
              (2.0 * g_iabc("k,a,c,d") - g_iabc("k,a,d,c")) * tau("c,d,k,i");

              r1("a,i") -=
              (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l");

              r1("a,i") *= d1("a,i");

              r1("a,i") -= t1("a,i");

              tmp_time1 = mpqc::now(world, accurate_time);
              tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
              if (this->print_detail_) {
                  mpqc::utility::print_par(world, "t1 other time: ", tmp_time, "\n");
              }
          }
          auto t1_time1 = mpqc::now(world, accurate_time);
          auto t1_time = mpqc::duration_in_s(t1_time0, t1_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t1 total time: ", t1_time, "\n");
          }

          // intermediates for t2
          // external index i j a b

          auto t2_time0 = mpqc::now(world, accurate_time);

          // compute residual r2(n) = t2(n+1) - t2(n)

          // permutation part
          tmp_time0 = mpqc::now(world, accurate_time);

          {
              r2("a,b,i,j") =
              (g_iabc("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j");

              r2("a,b,i,j") -=
              (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k");
          }
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 other time: ", tmp_time, "\n");
          }

          tmp_time0 = mpqc::now(world, accurate_time);
          {
              // compute g intermediate
              TArray g_ki, g_ac;

              g_ki("k,i") = h_ki("k,i") + f_ai("c,k") * t1("c,i") +
              (2.0 * g_ijka("k,l,i,c") - g_ijka("l,k,i,c")) * t1("c,l");

              g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k") +
              (2.0 * g_aibc("a,k,c,d") - g_aibc("a,k,d,c")) * t1("d,k");

              r2("a,b,i,j") +=
              g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j");
          }
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 g term time: ", tmp_time, "\n");
          }

          tmp_time0 = mpqc::now(world, accurate_time);
          {
              TArray j_akic;
              TArray k_kaic;
              // compute j and k intermediate
              {
                  TArray T;

                  T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

                  j_akic("a,k,i,c") = g_abij("a,c,i,k");

                  j_akic("a,k,i,c") -= g_ijka("l,k,i,c") * t1("a,l");

                  j_akic("a,k,i,c") += g_aibc("a,k,d,c") * t1("d,i");

                  j_akic("a,k,i,c") -= X_ai("x,d,l") * T("d,a,i,l") * X_ai("x,c,k");

                  j_akic("a,k,i,c") += 0.5 *
                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) *
                  t2("a,d,i,l");

                  k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                  - g_ijka("k,l,i,c") * t1("a,l")

                  + g_iabc("k,a,d,c") * t1("d,i")

                  - g_abij("d,c,k,l") * T("d,a,i,l");
                  if (this->print_detail_) {
                      mpqc::detail::print_size_info(T, "T");
                      mpqc::detail::print_size_info(j_akic, "J_akic");
                      mpqc::detail::print_size_info(k_kaic, "K_kaic");
                  }
              }

              r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
              (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

              r2("a,b,i,j") += -0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j") -
              k_kaic("k,b,i,c") * t2("a,c,k,j");
          }
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
          }

          // perform the permutation
          r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

          r2("a,b,i,j") += g_abij("a,b,i,j");

          tmp_time0 = mpqc::now(world, accurate_time);
          {
              TArray a_klij;
              // compute a intermediate
              a_klij("k,l,i,j") = g_ijkl("k,l,i,j");

              a_klij("k,l,i,j") += g_ijka("k,l,i,c") * t1("c,j");

              a_klij("k,l,i,j") += g_ijak("k,l,c,j") * t1("c,i");

              a_klij("k,l,i,j") += g_abij("c,d,k,l") * tau("c,d,i,j");

              r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");

              if (this->print_detail_) {
                  mpqc::detail::print_size_info(a_klij, "A_klij");
              }
          }
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
          }

          tmp_time0 = mpqc::now(world, accurate_time);
          {
              // compute b intermediate
              if (less) {
                  // avoid store b_abcd
                  TArray b_abij;
                  b_abij("a,b,i,j") = g_abcd("a,b,c,d") * tau("c,d,i,j");

                  b_abij("a,b,i,j") -= g_aibc("a,k,c,d") * tau("c,d,i,j") * t1("b,k");

                  b_abij("a,b,i,j") -= g_iabc("k,b,c,d") * tau("c,d,i,j") * t1("a,k");

                  if (this->print_detail_) {
                      mpqc::detail::print_size_info(b_abij, "B_abij");
                  }

                  r2("a,b,i,j") += b_abij("a,b,i,j");
              } else {
                  TArray b_abcd;

                  b_abcd("a,b,c,d") = g_abcd("a,b,c,d") -
                  g_aibc("a,k,c,d") * t1("b,k") -
                  g_iabc("k,b,c,d") * t1("a,k");

                  if (this->print_detail_) {
                      mpqc::detail::print_size_info(b_abcd, "B_abcd");
                  }

                  r2("a,b,i,j") += b_abcd("a,b,c,d") * tau("c,d,i,j");
              }
          }
          tmp_time1 = mpqc::now(world, accurate_time);
          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
          }

          d_abij_inplace(r2, *this->orbital_energy(), n_occ, n_frozen);

          r2("a,b,i,j") -= t2("a,b,i,j");

          t1("a,i") = t1("a,i") + r1("a,i");
          t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
          t1.truncate();
          t2.truncate();

          auto t2_time1 = mpqc::now(world, accurate_time);
          auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
          }
          // recompute energy
          E0 = E1;
          E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
          TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                  tau("a,b,i,j"));
          dE = std::abs(E0 - E1);

          double error_r1 = r1("a,i").norm().get();
          double error_r2 = r2("a,b,i,j").norm().get();

          if (dE >= this->converge_ || error >= this->converge_) {
              tmp_time0 = mpqc::now(world, accurate_time);
              cc::T1T2<double, Tile, Policy> t(t1, t2);
              cc::T1T2<double, Tile, Policy> r(r1, r2);
              error = r.norm() / size(t);
              diis.extrapolate(t, r);

              // update t1 and t2
              t1("a,i") = t.first("a,i");
              t2("a,b,i,j") = t.second("a,b,i,j");

              if (this->print_detail_) {
                  mpqc::detail::print_size_info(r2, "R2");
                  mpqc::detail::print_size_info(t2, "T2");
              }

              tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
              tmp_time1 = mpqc::now(world, accurate_time);
              tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
              if (this->print_detail_) {
                  mpqc::utility::print_par(world, "diis time: ", tmp_time, "\n");
              }

              auto time1 = mpqc::fenced_now(world);
              auto duration = mpqc::duration_in_s(time0, time1);

              if (world.rank() == 0) {
                  detail::print_ccsd(iter, dE, error, error_r1, error_r2, E1, duration);
              }

              iter += 1ul;
          } else {
              auto time1 = mpqc::fenced_now(world);
              auto duration = mpqc::duration_in_s(time0, time1);

              if (world.rank() == 0) {
                  detail::print_ccsd(iter, dE, error, error_r1, error_r2, E1, duration);
              }

              break;
          }
      }
      if (iter >= this->max_iter_) {
          utility::print_par(this->wfn_world()->world(),
                             "\n Warning!! Exceed Max Iteration! \n");
      }
      if (world.rank() == 0) {
          std::cout << "CCSD Energy  " << E1 << std::endl;
      }
      return E1;
    }

    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::value()  {

      auto &world = this->wfn_world()->world();

      double time;
      auto time0 = mpqc::fenced_now(world);

      double ref_energy = this->ref_wfn_->value();

      auto time1 = mpqc::fenced_now(world);
      time = mpqc::duration_in_s(time0, time1);
      utility::print_par(world, "Total Ref Time: ", time, " S \n");

      // initialize
      this->init();

      // compute MP2 T2 amplitudes
      ExEnv::out0() << std::endl << "Computing MP2 amplitudes" << std::endl;
      compute_mp2_t2();
      //ExEnv::out0() << "MP2 amplitudes: " << t2_mp2_ << std::endl;
      // copy MP2 amplitudes for test purpose
      TA::DistArray<Tile, Policy> t2_mp2_orig = t2_mp2_.clone();

      // reblock MP2 T2 amplitudes
      // in each block: n_i=n_j=1, n_a=n_b=nvir
      ExEnv::out0() << std::endl << "Reblocking MP2 amplitudes" << std::endl;
      // compute converting matrices
      TA::DistArray<Tile, Policy> occ_convert, vir_convert;
      compute_M_reblock(occ_convert, vir_convert);
      // obtain MP2 T2 with new blocking structure
      t2_mp2_("a,b,i,j") = t2_mp2_("c,d,k,l")
                         * vir_convert("c,a") * vir_convert("d,b")
                         * occ_convert("k,i") * occ_convert("l,j");
      //ExEnv::out0() << "MP2 amplitudes after reblocking: " << t2_mp2_ << std::endl;

      // decompose MP2 T2 amplitudes
      // using either eigen or SVD decomposition
      ExEnv::out0() << std::endl << "Decomposing MP2 amplitudes" << std::endl;
      decom_t2_mp2();

      // transform MP2 T2 amplitudes back into its original blocking structure
      t2_mp2_("a,b,i,j") = t2_mp2_("c,d,k,l")
                         * vir_convert("a,c") * vir_convert("b,d")
                         * occ_convert("i,k") * occ_convert("j,l");

      // compute the difference between original and decomposed MP2 T2
      ExEnv::out0() << std::endl << "Test: t2 - t2 (decomposed): "
                                 << (t2_mp2_orig("a,b,i,j") - t2_mp2_("a,b,i,j")).norm().get()
                                 << std::endl << std::endl;

      // compute CCSD with decomposed MP2 T2 as initial guess
      TA::DistArray<Tile, Policy> t1, t2;
      t2 = t2_mp2_.clone();
      double ccsd_energy = compute_ccsdpno_df(t1,t2);

      double energy = ref_energy + ccsd_energy;

      return energy;
    }

} // namespace lcao
} // namespace mpqc




#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_
