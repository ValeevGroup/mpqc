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

      template <typename Tile, typename Policy>
      double get_eccd_pno(const int nocc,
                          const std::vector<TA::DistArray<Tile, Policy>>& vec_t2_pno,
                          const std::vector<TA::DistArray<Tile, Policy>>& vec_gabij_pno) {

        double E_ccd = 0.0;
        for(int i = 0; i < nocc; ++i) {
          for (int j = 0; j < i; ++j) {
            const int idx = i * (i + 1) / 2 + j;
            TA::DistArray<Tile, Policy> tab = vec_t2_pno[idx];
            TA::DistArray<Tile, Policy> gab = vec_gabij_pno[idx];
            E_ccd += TA::dot(2.0 * gab("a,b") - gab("b,a"), tab("a,b")) * 2.0;
          }

          // add the diagonal contribution
          const int idx_ii = i * (i + 1) / 2 + i;
          TA::DistArray<Tile, Policy> tab_ii = vec_t2_pno[idx_ii];
          TA::DistArray<Tile, Policy> gab_ii = vec_gabij_pno[idx_ii];
          E_ccd += TA::dot(2.0 * gab_ii("a,b") - gab_ii("b,a"), tab_ii("a,b"));
        }

        return E_ccd;
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

    // compute MP2 T2 amplitudes
    template<typename Tile, typename Policy>
    TA::DistArray<Tile, Policy> CCSD_PNO<Tile, Policy>::compute_mp2_t2() {

      auto &lcao_factory = this->lcao_factory();
      return detail::compute_mp2_t2(lcao_factory, this->orbital_energy(),
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

    // obtain PNO coefficients
    template<typename Tile, typename Policy>
    TA::DistArray<Tile, Policy>
    CCSD_PNO<Tile, Policy>::compute_pno_coef(const TA::DistArray<Tile, Policy> &t2_mp2) {

      ExEnv::out0() << "tcut is " << tcut_ << std::endl;

      // decompose and reconstruct T2
      auto decom = [&](Tile &dab_tile, const Tile &tab_tile) -> float {

        // compute ij pair index
        const auto orb_i = tab_tile.range().upbound()[2];
        const auto orb_j = tab_tile.range().upbound()[3];

        // get the dimensions of input tile
        const auto extent = tab_tile.range().extent();
        const auto a = extent[0];
        const auto b = extent[1];

        auto Tab = TA::eigen_map(tab_tile, a, b);
        // compute virtual or PNO density
        TA::EigenMatrixXd Dab = Tab * (4.0 * Tab - 2.0 * Tab.transpose()).transpose()
                              + Tab.transpose() * (4.0 * Tab - 2.0 * Tab.transpose());
        // when i=j, Dab = Tab*Tab' + Tab'*Tab
        if (orb_i == orb_j)
          Dab = Dab/2.0;

        // compute eigenvalues ans eigenvectors of Dab
        Eigen::SelfAdjointEigenSolver<TA::EigenMatrixXd> es(Dab);
        TA::EigenVectorXd S_es = es.eigenvalues();
        TA::EigenMatrixXd C_es = es.eigenvectors();

        // truncate eigenvectors with corresponding eigenvalues smaller than threshold
        int pos_trunc = a - 1;
        for(; pos_trunc >= 0; --pos_trunc) {
          if(std::abs(S_es(pos_trunc)) < tcut_)
            break;
        }

        const size_t rank = a - pos_trunc - 1;
        if (rank > 0) {
          const size_t a_lo = tab_tile.range().lobound()[0];
          const size_t b_lo = tab_tile.range().lobound()[1];
          const size_t i_lo = tab_tile.range().lobound()[2];
          const size_t j_lo = tab_tile.range().lobound()[3];
          const size_t a_up = tab_tile.range().upbound()[0];
          const size_t b_up = rank;
          std::array<const size_t, 4> dab_lo = {a_lo, b_lo, i_lo, j_lo};
          std::array<const size_t, 4> dab_up = {a_up, b_up, orb_i, orb_j};

          dab_tile = Tile(TA::Range(dab_lo, dab_up));
          auto dab_pno = TA::eigen_map(dab_tile, a, rank);

          dab_pno.noalias() = C_es.rightCols(rank);

          // print out pair and rank information
          std::stringstream ss;
          ss << "("<< orb_i << "," << orb_j << ") pair  rank: " << rank
             << std::endl;
          //ss << "dab tile: " << dab_tile << std::endl;
          std::printf("%s", ss.str().c_str());
        }

        return dab_tile.norm();
      }; // end of decom function

      auto result = TA::foreach (t2_mp2, decom);
      return result;
    }

    // obtain the PNO coefficients in a vector
    // each vector entry contain ab matrix
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::comput_vec_pnocoef(const TA::DistArray<Tile, Policy>& t2_mp2,
                                                    const TA::TiledRange1& trange1_a) {

      TA::DistArray<Tile, Policy> dab_ij = compute_pno_coef(t2_mp2);

      // obtain the size of the vector
      const int nocc = t2_mp2.trange().elements_range().extent()[2];
      const int size_vecij = nocc * (nocc + 1) / 2;

      vecij_dab_.resize(size_vecij);
      comput_vecij_ab(dab_ij, trange1_a, vecij_dab_, true);
    }

    // construct vector<TArray> for a tensor abij,
    // for abij, tile is indexed by ij and has size a*b
    // each vector entry contain ab matrix
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::comput_vecij_ab(const TA::DistArray<Tile, Policy>& ab_ij,
                                                 const TA::TiledRange1& trange1_a,
                                                 std::vector<TA::DistArray<Tile, Policy>>& vecij_dab,
                                                 const bool compute_pno_coef) {

       auto &world = this->wfn_world()->world();
       const int ni = ab_ij.trange().elements_range().extent()[2];
       const int nij = ni * ni;

       // loop over tiles
       for(int index = 0; index < nij; ++index) {
         if(ab_ij.is_zero(index))
           continue;

         const Tile ab_tile = ab_ij.find(index);

         // only loop over i>=j pairs
         // skip i < j pairs
         const int orb_i = ab_tile.range().lobound()[2]; // index i
         const int orb_j = ab_tile.range().lobound()[3]; // index j
         if (orb_i < orb_j)
           continue;

         const int b = ab_tile.range().upbound()[1];

         // create TiledRange for result tensor
         TA::TiledRange trange_ab;
         if (compute_pno_coef) {
           // create TileRange1 for PNO dimension
           const int size_pno_blk = 4;
           const int n_pno_blk = b / size_pno_blk;
           const int n_last_blk = b % size_pno_blk;

           std::vector<int> tr_pno;
           tr_pno.push_back(0);
           for (int i = 1; i < n_pno_blk; ++i) {
             tr_pno.push_back(i*size_pno_blk);
           }
           tr_pno.push_back(b);

           TA::TiledRange1 trange1_pno(tr_pno.begin(), tr_pno.end());


           trange_ab = {trange1_a, trange1_pno};
         } else {
           trange_ab = {trange1_a, trange1_a};
         }

         const TA::DistArray<Tile, Policy>
         array_from_tile = TA::make_array<TA::DistArray<Tile, Policy>>
                               (world, trange_ab,
                                [=] (Tile& tile, const TA::Range& range) -> double {
                                  std::array<std::size_t, 2> lobound{ ab_tile.range().lobound_data()[0],
                                                                      ab_tile.range().lobound_data()[1] };
                                  std::array<std::size_t, 2> upbound{ ab_tile.range().upbound_data()[0],
                                                                      ab_tile.range().upbound_data()[1] };
                                  TA::Range dab_range(lobound, upbound);
                                  TA::detail::TensorInterface<const double, TA::BlockRange>
                                  t_ab(TA::BlockRange(dab_range, range.lobound(), range.upbound()),
                                      ab_tile.data());

                                  tile = Tile(range);
                                  tile = t_ab;
//                                tile = TA::Tensor<double>(range);
//                                for(std::size_t i = 0; i < tile.size(); ++i)
//                                  tile[i] = 1.0;
                                return tile.norm();}
                               );

         // compute the index for the vector:
         const int idx = orb_i * (orb_i + 1) / 2 + orb_j;
         vecij_dab[idx] = array_from_tile;

         // test: print out information
//         std::stringstream ss;
//         ss << std::endl;
//         ss << "ab_tile: " << ab_tile << std::endl;
//         ss << "i = " <<  orb_i << "  j = " << orb_j << "  idx (vec) = " << idx
//            << "  b = " << b << std::endl;
//         ss << "ab range: " << ab_tile.range() << std::endl
//            << "  trange_ab: " << trange_ab << std::endl;
//         ss << "array_from_tile: " << vecij_dab[idx] << std::endl;
//         std::printf("%s", ss.str().c_str());
       } // end of loop over tile

     }

    // compute PNO transformed F^a_b
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::get_fab_pno(std::vector<TA::DistArray<Tile, Policy>>& vecij_fab) {

      if (!vecij_dab_.empty()) {
        // get regular fock integrals
        TA::DistArray<Tile, Policy> fab;
        if (this->df_) {
          fab = this->lcao_factory().compute(L"<a|F|b>[df]");
        } else {
          fab = this->lcao_factory().compute(L"<a|F|b>");
        }

        int index = 0;
        for (const auto& ele_ij: vecij_dab_) {
          const int rank_vir = ele_ij.trange().elements_range().extent()[0];
          const int rank_pno = ele_ij.trange().elements_range().extent()[1];
          TA::DistArray<Tile, Policy> fab_pno;
//          if (rank_vir == rank_pno) {
//            fab_pno = fab;
//          } else {
            fab_pno("a,b") = fab("e,f") * ele_ij("e,a") * ele_ij("f,b");
//          }
          vecij_fab[index] = fab_pno;
          ++index;
        }
      } else {
        throw std::runtime_error("PNO coefficients have not been computed yet!");
      }
    }

    // compute PNO transformed T2^ab_ij or G^ab_ij
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::get_abij_pno(const TA::DistArray<Tile, Policy>& ab_ij,
                                              const TA::TiledRange1& trange1_a,
                                              std::vector<TA::DistArray<Tile, Policy>>& vecij_pno) {

      // reblock abij array
      const int size_vecij = vecij_pno.size();
      std::vector<TA::DistArray<Tile, Policy>> vecij_ab(size_vecij);
      comput_vecij_ab(ab_ij, trange1_a, vecij_ab, false);

      if (!vecij_dab_.empty()) {
        int index = 0;
        for (const auto& ele_ij: vecij_dab_) {
          const int rank_vir = ele_ij.trange().elements_range().extent()[0];
          const int rank_pno = ele_ij.trange().elements_range().extent()[1];
//          ExEnv::out0() << std::endl << "rank_vir: " << rank_vir
//                        << "  rank_vir: " << rank_pno << std::endl;

          // obtain T2_ab
          const TA::DistArray<Tile, Policy> ab = vecij_ab[index];
          TA::DistArray<Tile, Policy> ab_pno;
//          if (rank_vir == rank_pno) {
//            ab_pno = ab;
//          } else {
            ab_pno("a,b") = ab("c,d") * ele_ij("c,a") * ele_ij("d,b");
//          }
          vecij_pno[index] = ab_pno;
          ++index;
        }
      } else {
        throw std::runtime_error("PNO coefficients have not been computed yet!");
      }
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::get_abcd_pno(std::vector<TA::DistArray<Tile, Policy>>& vecij_gabcd) {

      if (!vecij_dab_.empty()) {
        // get regular two electron integrals
        TA::DistArray<Tile, Policy> g_abcd = this->get_abcd();

        int index = 0;
        for (const auto& ele_ij: vecij_dab_) {
          const int rank_vir = ele_ij.trange().elements_range().extent()[0];
          const int rank_pno = ele_ij.trange().elements_range().extent()[1];
          TA::DistArray<Tile, Policy> gabcd_pno;
//          if (rank_vir == rank_pno) {
//            gabcd_pno = g_abcd;
//          } else {
            gabcd_pno("a,b,c,d") = g_abcd("e,f,g,h") * ele_ij("e,a") * ele_ij("f,b")
                                   * ele_ij("g,c") * ele_ij("h,d");
//          }
          vecij_gabcd[index] = gabcd_pno;
          ++index;
        }
      } else {
        throw std::runtime_error("PNO coefficients have not been computed yet!");
      }
    }

    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::compute_ccd_pno() {

      using TArray = TA::DistArray<Tile, Policy>;

      auto &world = this->wfn_world()->world();
      bool accurate_time = this->lcao_factory().accurate_time();

      // obtain MP2 T2 amplitudes
      ExEnv::out0() << std::endl << "Computing MP2 amplitudes" << std::endl;
      TArray t2_mp2 = compute_mp2_t2();

      // compute converting matrices for reblocking T2 amplitudes
      // in each block: n_i=n_j=1, n_a=n_b=nvir
      TArray occ_convert, vir_convert;
      compute_M_reblock(occ_convert, vir_convert);

//       // test function: decompose PNO with different receipts
//       double corr_energy = test_pno_decom(false, t2_mp2, occ_convert, vir_convert);

      ExEnv::out0() << std::endl << "Reblocking T2 amplitudes" << std::endl;
      // obtain MP2 T2 with new blocking structure
      t2_mp2("a,b,i,j") = t2_mp2("c,d,k,l")
                        * vir_convert("c,a") * vir_convert("d,b")
                        * occ_convert("k,i") * occ_convert("l,j");

      // obtain pno coefficients in a vector with each element as dab matrix
      // dab (b: PNO index) is indexed by ij (i>=j) and b is reblocked
      // obtain tile range for virtual index: a
      const TA::TiledRange1 trange1_a = vir_convert.trange().data()[0];
      comput_vec_pnocoef(t2_mp2, trange1_a);

      // obtain the size of vectors for PNO transformed tensors
      const int nocc = occ_convert.trange().elements_range().extent()[1];
      const int size_vecij = nocc * (nocc + 1) / 2;
//      ExEnv::out0() << std::endl << "nocc: " << nocc << "  size_vecij: " << size_vecij
//                    << "  trange1_a: " << trange1_a << std::endl;

      // compute PNO integrals: g^ab_cd
      ExEnv::out0() << std::endl << "Computing PNO transformed g^ab_cd" << std::endl;
      std::vector<TArray> vec_gabcd_pno(size_vecij);
      get_abcd_pno(vec_gabcd_pno);

      // compute PNO t2^ab_ij
      ExEnv::out0() << "Computing PNO transformed t2^ab_ij" << std::endl;
      std::vector<TArray> vec_t2_pno(size_vecij);
      get_abij_pno(t2_mp2, trange1_a, vec_t2_pno);

      // compute PNO integrals: g^ab_ij
      ExEnv::out0() << "Computing PNO transformed g^ab_ij" << std::endl;
      TA::DistArray<Tile, Policy> g_abij = this->get_abij();
      // obtain g_abij with new blocking structure
      g_abij("a,b,i,j") = g_abij("c,d,k,l")
                        * vir_convert("c,a") * vir_convert("d,b")
                        * occ_convert("k,i") * occ_convert("l,j");
      std::vector<TArray> vec_gabij_pno(size_vecij);
      get_abij_pno(g_abij, trange1_a, vec_gabij_pno);

      // compute PNO transformed f^a_b
      ExEnv::out0() << "Computing PNO transformed f^a_b" << std::endl;
      std::vector<TArray> vec_fab_pno(size_vecij);
      get_fab_pno(vec_fab_pno);
//      ExEnv::out0() << "vec_fab_pno: " << std::endl;
//      for (int i = 0; i < vec_fab_pno.size(); ++i) {
//      ExEnv::out0() << i << " " << vec_fab_pno[i] << std::endl;
//      }

      // compute CCD energy
      ExEnv::out0() << std::endl << "Computing initial PNO-CCD energy" << std::endl;
      double E0 = 0.0;
      double E1 = detail::get_eccd_pno(nocc, vec_t2_pno, vec_gabij_pno);
      double dE = std::abs(E0 - E1);
      mpqc::utility::print_par(world, "Initial CCD (PNO) Energy      ", E1, "\n");

      // obtain orbital energy
      const Eigen::VectorXd e_orb = *this->orbital_energy();

      const int nij = nocc * nocc;
      std::vector<TArray> vecij_r2(size_vecij);
      TArray r2;

      if (world.rank() == 0)
        std::printf("%3s \t %10s \t %10s \t %15s \t %10s \n", "iter", "deltaE",
                  "residual", "energy", "total time/s");
      int iter = 0;
      while (iter < this->max_iter_) {
//      while (iter < 2) {
        auto time0 = mpqc::fenced_now(world);
        TArray::wait_for_lazy_cleanup(world);

        for(int i = 0; i < nocc; ++i) {
          for (int j = 0; j <=i; ++j) {
            const int idx = i * (i + 1) / 2 + j;

            auto t2_time0 = mpqc::now(world, accurate_time);
            // obtain t2_ab
            TArray tab = vec_t2_pno[idx];
            TArray gab = vec_gabij_pno[idx];
            TArray gabcd = vec_gabcd_pno[idx];
            // r2("a,b") = gab("a,b") + gabcd("a,b,c,d") * tab("c,d");
            r2 = gab.clone();

            TArray fab_ta = vec_fab_pno[idx];
            auto fab = array_ops::array_to_eigen(fab_ta);
            const Eigen::VectorXd faa = fab.diagonal();
            const double e_i = e_orb[i];
            const double e_j = e_orb[j];
            d_ab_inplace(r2, faa, e_i, e_j);

            r2("a,b") -= tab("a,b");
            tab("a,b") = tab("a,b") + r2("a,b");
            tab.truncate();

            // update amplitudes and residual
            vec_t2_pno[idx] = tab;
            vecij_r2[idx] = r2;

            auto t2_time1 = mpqc::now(world, accurate_time);
            auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
            if (this->print_detail_) {
                mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
            }
          }
        }

        // update CCD energy
        E0 = E1;
        E1 = detail::get_eccd_pno(nocc, vec_t2_pno, vec_gabij_pno);
        dE = std::abs(E0 - E1);

        // compute the norm of all the residuals
        const double norm_r2 = norm_vec_tensor(vecij_r2);
        if (dE > this->converge_ || norm_r2 > this->converge_) {
          if (world.rank() == 0)
            std::printf("%3i \t %10.5e \t %10.5e \t %15.12f \t %10.1f \n", iter, dE, norm_r2, E1,
                      time);
            ++iter;
        } else {
            break;
        }

      } // end of CCD iteration

      if (iter >= this->max_iter_) {
          utility::print_par(this->wfn_world()->world(),
                             "\n Warning!! Exceed Max Iteration! \n");
      }

      if (world.rank() == 0) {
          std::cout << "CCD Energy  " << E1 << std::endl;
      }

      return E1;
    }

    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::compute_ccd_df() {

      using TArray = TA::DistArray<Tile, Policy>;

      auto &world = this->wfn_world()->world();
      bool accurate_time = this->lcao_factory().accurate_time();

      auto n_occ = this->trange1_engine()->get_occ();
      auto n_frozen = this->trange1_engine()->get_nfrozen();

      if (world.rank() == 0) {
          std::cout << "Start DF CCD Computation" << std::endl;
      }

      auto tmp_time0 = mpqc::now(world, accurate_time);
      // get all two electron integrals
      TArray g_abij = this->get_abij();
//      TArray g_ijkl = this->get_ijkl();
      TArray g_abcd = this->get_abcd();
//      TArray X_ai = this->get_Xai();
//      TArray g_iajb = this->get_iajb();
//      TArray g_iabc = this->get_iabc();
//      TArray g_aibc = this->get_aibc();
//      TArray g_ijak = this->get_ijak();
//      TArray g_ijka = this->get_ijka();
      TArray f_ai = this->get_fock_ai();
      auto tmp_time1 = mpqc::now(world, accurate_time);
      auto tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
      if (this->print_detail_) {
          mpqc::utility::print_par(world, "Integral Prepare Time: ", tmp_time,
                                   "\n");
      }

      TArray t2 = compute_mp2_t2();
      double E0 = 0.0;
      double E1 = TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), t2("a,b,i,j"));
      double dE = std::abs(E1 - E0);

      mpqc::utility::print_par(world, "Initial CCD Energy      ", E1, "\n");

      // optimize t1 and t2
      std::size_t iter = 0ul;
      double error = 1.0;
      TArray r2;

      if (world.rank() == 0) {
          std::cout << "Start Iteration" << std::endl;
          std::cout << "Max Iteration: " << this->max_iter_ << std::endl;
          std::cout << "Convergence: " << this->converge_ << std::endl;
          std::cout << "AccurateTime: " << accurate_time << std::endl;
          std::cout << "PrintDetail: " << this->print_detail_ << std::endl;
      }

      TArray t1(world, f_ai.trange(), f_ai.get_shape());
      TArray r1(world, f_ai.trange(), f_ai.get_shape());
      t1.fill(0.0);
      r1.fill(0.0);

      auto diis = this->get_diis(world);
      while (iter < this->max_iter_) {
          // start timer
          auto time0 = mpqc::fenced_now(world);
          TArray::wait_for_lazy_cleanup(world);

//          auto t1_time0 = mpqc::now(world, accurate_time);
//          TArray h_ki, h_ac;
//          {
//              // intermediates for t1
//              // external index i and a
//              // vir index a b c d
//              // occ index i j k l
//
//              // compute residual r1(n) = t1(n+1) - t1(n)
//              // external index i and a
//              tmp_time0 = mpqc::now(world, accurate_time);
//
//              {
//                  h_ac("a,c") =
//                  -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * t2("a,d,k,l");
//              }
//
//              {
//                  h_ki("k,i") =
//                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t2("c,d,i,l");
//              }
//
//              tmp_time1 = mpqc::now(world, accurate_time);
//              tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
//              if (this->print_detail_) {
//                  mpqc::utility::print_par(world, "t1 h term time: ", tmp_time, "\n");
//              }
//          }

          // intermediates for t2
          // external index i j a b

          auto t2_time0 = mpqc::now(world, accurate_time);

          // compute residual r2(n) = t2(n+1) - t2(n)

//          tmp_time0 = mpqc::now(world, accurate_time);
//          {
//              // compute g intermediate
//              TArray g_ki, g_ac;
//
//              r2("a,b,i,j") = h_ac("a,c") * t2("c,b,i,j") - h_ki("k,i") * t2("a,b,k,j");
//          }
//          tmp_time1 = mpqc::now(world, accurate_time);
//          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
//          if (this->print_detail_) {
//              mpqc::utility::print_par(world, "t2 h term time: ", tmp_time, "\n");
//          }

//          tmp_time0 = mpqc::now(world, accurate_time);
//          {
//              TArray j_akic;
//              TArray k_kaic;
//              // compute j and k intermediate
//              {
//                  TArray T;
//
//                  T("d,b,i,l") = 0.5 * t2("d,b,i,l");
//
//                  j_akic("a,k,i,c") = g_abij("a,c,i,k");
//
//                  j_akic("a,k,i,c") -= X_ai("x,d,l") * T("d,a,i,l") * X_ai("x,c,k");
//
//                  j_akic("a,k,i,c") += 0.5 *
//                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) *
//                  t2("a,d,i,l");
//
//                  k_kaic("k,a,i,c") = g_iajb("k,a,i,c") - g_abij("d,c,k,l") * T("d,a,i,l");
//                  if (this->print_detail_) {
//                      mpqc::detail::print_size_info(T, "T");
//                      mpqc::detail::print_size_info(j_akic, "J_akic");
//                      mpqc::detail::print_size_info(k_kaic, "K_kaic");
//                  }
//              }
//
//              r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
//              (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));
//
//              r2("a,b,i,j") += -0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j") -
//              k_kaic("k,b,i,c") * t2("a,c,k,j");
//          }
//          tmp_time1 = mpqc::now(world, accurate_time);
//          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
//          if (this->print_detail_) {
//              mpqc::utility::print_par(world, "t2 j,k term time: ", tmp_time, "\n");
//          }

//          // perform the permutation
//          r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");
//
//          r2("a,b,i,j") += g_abij("a,b,i,j");
          r2 = g_abij.clone();

//          tmp_time0 = mpqc::now(world, accurate_time);
//          {
//              TArray a_klij;
//              // compute a intermediate
//              a_klij("k,l,i,j") = g_ijkl("k,l,i,j");
//
//              a_klij("k,l,i,j") += g_abij("c,d,k,l") * t2("c,d,i,j");
//
//              r2("a,b,i,j") += a_klij("k,l,i,j") * t2("a,b,k,l");
//
//              if (this->print_detail_) {
//                  mpqc::detail::print_size_info(a_klij, "A_klij");
//              }
//          }
//          tmp_time1 = mpqc::now(world, accurate_time);
//          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
//          if (this->print_detail_) {
//              mpqc::utility::print_par(world, "t2 a term time: ", tmp_time, "\n");
//          }

//          tmp_time0 = mpqc::now(world, accurate_time);
//          {
//            // compute b intermediate
//            TArray b_abcd;
//
//            b_abcd("a,b,c,d") = g_abcd("a,b,c,d");
//
//            if (this->print_detail_) {
//                mpqc::detail::print_size_info(b_abcd, "B_abcd");
//            }
//
//            r2("a,b,i,j") += b_abcd("a,b,c,d") * t2("c,d,i,j");
//          }
//          tmp_time1 = mpqc::now(world, accurate_time);
//          tmp_time = mpqc::duration_in_s(tmp_time0, tmp_time1);
//          if (this->print_detail_) {
//              mpqc::utility::print_par(world, "t2 b term time: ", tmp_time, "\n");
//          }

          d_abij_inplace(r2, *this->orbital_energy(), n_occ, n_frozen);

          r2("a,b,i,j") -= t2("a,b,i,j");
          t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");
          t2.truncate();

          auto t2_time1 = mpqc::now(world, accurate_time);
          auto t2_time = mpqc::duration_in_s(t2_time0, t2_time1);
          if (this->print_detail_) {
              mpqc::utility::print_par(world, "t2 total time: ", t2_time, "\n");
          }
          // recompute energy
          E0 = E1;
          E1 = TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), t2("a,b,i,j"));
          dE = std::abs(E0 - E1);

          double error_r1 = 0.0;
          double error_r2 = r2("a,b,i,j").norm().get();

          if (dE >= this->converge_ || error >= this->converge_) {
              tmp_time0 = mpqc::now(world, accurate_time);
              cc::T1T2<TArray, TArray> t(t1, t2);
              cc::T1T2<TArray, TArray> r(r1, r2);
              error = r.norm() / size(t2);  // error = residual norm per element
              diis.extrapolate(t, r);

              // update t1 and t2
              t2("a,b,i,j") = t.t2("a,b,i,j");

              if (this->print_detail_) {
                  mpqc::detail::print_size_info(r2, "R2");
                  mpqc::detail::print_size_info(t2, "T2");
              }

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
          std::cout << "CCD Energy  " << E1 << std::endl;
      }
      return E1;
    }

    // decompose T2 amplitudes
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::decom_t2(TA::DistArray<Tile, Policy> &t2_mp2,
                                          const TA::DistArray<Tile, Policy> &t2_ccsd,
                                          bool use_diff_t2) {

      double threshold = tcut_;
      ExEnv::out0() << "tcut is " << tcut_ << std::endl;

      // decompose and reconstruct T2
      auto decom_1 = [&](Tile &in_tile) -> float {

        // copy in_tile which is ab matrix of T^ij
        Tile tab_tile = in_tile.clone();

        // compute ij pair index
        const auto orb_i = in_tile.range().upbound()[2];
        const auto orb_j = in_tile.range().upbound()[3];

        // get the dimensions of input tile
        const auto extent = in_tile.range().extent();
        const auto a = extent[0];
        const auto b = extent[1];

        int rank = 0;

        std::stringstream ss;
//        ss << "input tile: " << in_tile << std::endl;

#define PNO_DECOM 1
#define T2_SVD 0

#if PNO_DECOM
        auto Tab = TA::eigen_map(tab_tile, a, b);
        // compute virtual or PNO density
        TA::EigenMatrixXd Dab = Tab * (4.0 * Tab - 2.0 * Tab.transpose()).transpose()
                              + Tab.transpose() * (4.0 * Tab - 2.0 * Tab.transpose());
        // when i=j, Dab = Tab*Tab' + Tab'*Tab
        if (orb_i == orb_j)
          Dab = Dab/2.0;

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

        rank = a - pos_trunc - 1;
        if (rank > 0) {
          auto tab_pno = TA::eigen_map(in_tile, a, b);
          tab_pno.noalias() = C_es.rightCols(rank) * C_es.rightCols(rank).transpose()
                            * Tab
                            * C_es.rightCols(rank) * C_es.rightCols(rank).transpose();
        } else {
            std::fill(in_tile.begin(), in_tile.end(), 0.0);
        }

//        ss << "Dab: " << std::endl << Dab << std::endl;
//        ss << "S_es: " << std::endl << S_es << std::endl;
//        ss << "C_es: " << std::endl << C_es << std::endl;
//        ss << "S_es (trunc): " << std::endl << S_es.tail(rank) << std::endl;
//        ss << "C_es (trunc): " << std::endl << C_es.rightCols(rank) << std::endl;
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
#endif // T2_SVD

//        ss << "output tile: " << in_tile << std::endl;

        // compute the norm of difference between input and output tile
        ss << "("<< orb_i << "," << orb_j << ") pair  rank: " << rank
           << "  (t_new - t_old).norm(): " << (TA::subt(in_tile, tab_tile)).norm()
           << std::endl;

        std::printf("%s", ss.str().c_str());

        return in_tile.norm();
      };

      // decompose and reconstruct T2
      auto decom_2 = [&](Tile &tile_mp2, const Tile &tile_ccsd) -> float {

        // copy in_tile which is ab matrix of T^ij
        Tile tab_mp2 = tile_mp2.clone();

        // compute ij pair index
        const auto orb_i = tile_mp2.range().upbound_data()[2];
        const auto orb_j = tile_mp2.range().upbound_data()[3];

        // get the dimensions of input tile
        const auto extent = tile_mp2.range().extent();
        const auto a = extent[0];
        const auto b = extent[1];

        int rank = 0;

        // PNO decomposition using MP2 T2 amplitudes
        auto Tab_mp2 = TA::eigen_map(tab_mp2, a, b);
        // compute virtual or PNO density
        TA::EigenMatrixXd Dab_mp2 = Tab_mp2 * (4.0 * Tab_mp2 - 2.0 * Tab_mp2.transpose()).transpose()
                                  + Tab_mp2.transpose() * (4.0 * Tab_mp2 - 2.0 * Tab_mp2.transpose());
        // when i=j, Dab = Tab*Tab' + Tab'*Tab
        if (orb_i == orb_j)
          Dab_mp2 = Dab_mp2/2.0;

        // compute eigenvalues ans eigenvectors of Dab
        Eigen::SelfAdjointEigenSolver<TA::EigenMatrixXd> es_mp2(Dab_mp2);
        TA::EigenVectorXd S_es_mp2 = es_mp2.eigenvalues();
        TA::EigenMatrixXd C_es_mp2 = es_mp2.eigenvectors();

        // PNO decomposition using CCSD T2 amplitudes
        auto Tab_ccsd = TA::eigen_map(tile_ccsd, a, b);
        // compute virtual or PNO density
        TA::EigenMatrixXd Dab_ccsd = Tab_ccsd * (4.0 * Tab_ccsd - 2.0 * Tab_ccsd.transpose()).transpose()
                                   + Tab_ccsd.transpose() * (4.0 * Tab_ccsd - 2.0 * Tab_ccsd.transpose());
        // when i=j, Dab = Tab*Tab' + Tab'*Tab
        if (orb_i == orb_j)
          Dab_ccsd = Dab_ccsd/2.0;

        // compute eigenvalues ans eigenvectors of Dab
        Eigen::SelfAdjointEigenSolver<TA::EigenMatrixXd> es_ccsd(Dab_ccsd);
        TA::EigenVectorXd S_es_ccsd = es_ccsd.eigenvalues();
        TA::EigenMatrixXd C_es_ccsd = es_ccsd.eigenvectors();

        // truncate eigenvectors with corresponding eigenvalues smaller than threshold
        int pos_trunc = a - 1;
        for(; pos_trunc >= 0; --pos_trunc) {
          if(std::abs(S_es_mp2(pos_trunc)) < threshold)
            break;
        }

        rank = a - pos_trunc - 1;
        // use different C and Tab to construct output tile
        if (rank > 0) {
          auto tab_pno = TA::eigen_map(tile_mp2, a, b);
          tab_pno.noalias() = C_es_mp2.rightCols(rank) * C_es_ccsd.rightCols(rank).transpose()
                            * Tab_ccsd
                            * C_es_ccsd.rightCols(rank) * C_es_mp2.rightCols(rank).transpose();
        } else {
            std::fill(tile_mp2.begin(), tile_mp2.end(), 0.0);
        }

        std::stringstream ss;
//        // test: print out all the data
//        if (orb_i == 2 && orb_j == 4) {
//          ss << "input tile: " << tab_mp2 << std::endl;
//          ss << "Tab_mp2: " << std::endl << Tab_mp2 << std::endl;
//          ss << "Dab_mp2: " << std::endl << Dab_mp2 << std::endl;
//          ss << "S_es_mp2 (trunc): " << std::endl << S_es_mp2.tail(rank) << std::endl;
//          ss << "C_es_mp2 (trunc): " << std::endl << C_es_mp2.rightCols(rank) << std::endl;
//
//          ss << "Tab_ccsd: " << std::endl << Tab_ccsd << std::endl;
//          ss << "Dab_ccsd: " << std::endl << Dab_ccsd << std::endl;
//          ss << "S_es_ccsd (trunc): " << std::endl << S_es_ccsd.tail(rank) << std::endl;
//          ss << "C_es_ccsd (trunc): " << std::endl << C_es_ccsd.rightCols(rank) << std::endl;
//
//          ss << "output tile: " << tile_mp2 << std::endl;
//        }

        // compute the norm of difference between input and output tile
        ss << "("<< orb_i << "," << orb_j << ") pair  rank: " << rank
           << "  (t_new - t_old).norm(): " << (TA::subt(tile_mp2, tab_mp2)).norm()
           << std::endl;

        std::printf("%s", ss.str().c_str());

        return tile_mp2.norm();
      };

      if (use_diff_t2) {
        ExEnv::out0() << std::endl << "Using MP2 and CCSD T2 amplitudes for PNO decomposition" << std::endl;
        TA::foreach_inplace(t2_mp2, t2_ccsd, decom_2);
      } else {
        ExEnv::out0() << std::endl << "Using MP2 T2 amplitudes for PNO decomposition" << std::endl;
        TA::foreach_inplace(t2_mp2, decom_1);
      }
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
              cc::T1T2<TArray, TArray> t(t1, t2);
              cc::T1T2<TArray, TArray> r(r1, r2);
              error = r.norm() / (size(t1) + size(t2));  // error = residual norm per element
              diis.extrapolate(t, r);

              // update t1 and t2
              t1("a,i") = t.t1("a,i");
              t2("a,b,i,j") = t.t2("a,b,i,j");

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

    // test PNO decomposition
    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::test_pno_decom(const bool use_diff_t2,
                                                  const TA::DistArray<Tile, Policy> &t2_mp2_orig,
                                                  const TA::DistArray<Tile, Policy> &occ_convert,
                                                  const TA::DistArray<Tile, Policy> &vir_convert) {

      // obtain MP2 T2 with new blocking structure
      TA::DistArray<Tile, Policy> t2_mp2 = t2_mp2_orig.clone();
      t2_mp2("a,b,i,j") = t2_mp2("c,d,k,l")
                        * vir_convert("c,a") * vir_convert("d,b")
                        * occ_convert("k,i") * occ_convert("l,j");

      // compute CCSD T2 amplitudes
      TA::DistArray<Tile, Policy> t1_ccsd, t2_ccsd;
      if (use_diff_t2) {
        t2_ccsd = t2_mp2_orig.clone();
        double ccsd_energy = compute_ccsdpno_df(t1_ccsd, t2_ccsd);

        // obtain CCSD T2 with new blocking structure
        t2_ccsd("a,b,i,j") = t2_ccsd("c,d,k,l")
                           * vir_convert("c,a") * vir_convert("d,b")
                           * occ_convert("k,i") * occ_convert("l,j");
      }

      // decompose T2 amplitudes
      // using either eigen or SVD decomposition
      ExEnv::out0() << "Decomposing T2 amplitudes" << std::endl;
//      // test
//      TA::DistArray<Tile, Policy> t2_mp2_orig_blk = t2_mp2.clone();
//      decom_t2(t2_mp2, t2_mp2_orig_blk, use_diff_t2);
      decom_t2(t2_mp2, t2_ccsd, use_diff_t2);

      // transform MP2 T2 amplitudes back into its original blocking structure
      t2_mp2("a,b,i,j") = t2_mp2("c,d,k,l")
                         * vir_convert("a,c") * vir_convert("b,d")
                         * occ_convert("i,k") * occ_convert("j,l");

      // compute the difference between original and decomposed MP2 T2
      ExEnv::out0() << std::endl << "Test: t2 - t2 (decomposed): "
                                 << (t2_mp2_orig("a,b,i,j") - t2_mp2("a,b,i,j")).norm().get()
                                 << std::endl << std::endl;

      // compute CCSD with decomposed MP2 T2 as initial guess
      TA::DistArray<Tile, Policy> t1, t2;
      t2 = t2_mp2.clone();

      return compute_ccsdpno_df(t1,t2);
    }

    // computing function
    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::evaluate(Energy* result) {

      if (!this->computed()) {
        /// cast ref_wfn to Energy::Evaluator
        auto ref_evaluator = std::dynamic_pointer_cast<typename Energy::Evaluator>(this->ref_wfn_);
        if(ref_evaluator == nullptr) {
          std::ostringstream oss;
          oss << "RefWavefunction in CCSD-PNO" << this->ref_wfn_->class_key()
              << " cannot compute Energy" << std::endl;
          throw InputError(oss.str().c_str(), __FILE__, __LINE__);
        }

        ref_evaluator->evaluate(result);

        double ref_energy = this->get_value(result).derivs(0)[0];

        // initialize
        this->init();

        // test CCD equation
        double ccd_energy = compute_ccd_df();

        double corr_energy = compute_ccd_pno();

        this->computed_ = true;
        this->set_value(result, ref_energy + corr_energy);
      }
    }

} // namespace lcao
} // namespace mpqc

#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_
