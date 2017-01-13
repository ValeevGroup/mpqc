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
    }

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
    CCSD_PNO<Tile, Policy>::CCSD_PNO(const KeyVal &kv) : LCAOWavefunction<Tile, Policy>(kv) {
      if (kv.exists("ref")) {
        ref_wfn_ = kv.keyval("ref").class_ptr<Wavefunction>();
      } else {
        throw std::invalid_argument("Default Ref Wfn in CCSD is not support! \n");
      }

      df_ = kv.value<bool>("df", false);
      tcut_ = kv.value<double>("tcut", 1e-6);
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::compute(PropertyBase *pb)  {
      throw std::runtime_error("Not Implemented!!");
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::init() {
      if (this->trange1_engine_ == nullptr || this->orbital_energy_ == nullptr) {
        auto mol = this->lcao_factory().ao_factory().molecule();
        Eigen::VectorXd orbital_energy;
        this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(this->lcao_factory(),
                                                                      orbital_energy, mol,
                                                                      this->is_frozen_core(),
                                                                      this->occ_block(),
                                                                      this->unocc_block());
        this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
      }
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::compute_mp2_t2() {

      auto &lcao_factory = this->lcao_factory();
      t2_mp2_ = detail::compute_mp2_t2(lcao_factory, this->orbital_energy(),
                                    this->trange1_engine(), df_);
    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::reblock() {
      auto &lcao_factory = this->lcao_factory();
      auto &world = lcao_factory.world();

      std::size_t occ = this->trange1_engine_->get_occ();
      std::size_t vir = this->trange1_engine_->get_vir();
      std::size_t all = this->trange1_engine_->get_all();
      std::size_t n_frozen = this->trange1_engine_->get_nfrozen();

      std::size_t b_occ = 1;
      std::size_t b_vir = vir;

      TA::TiledRange1 old_occ = this->trange1_engine_->get_active_occ_tr1();
      TA::TiledRange1 old_vir = this->trange1_engine_->get_vir_tr1();

      auto new_tr1 =
          std::make_shared<TRange1Engine>(occ, all, b_occ, b_vir, n_frozen);

      TA::TiledRange1 new_occ = new_tr1->get_active_occ_tr1();
      TA::TiledRange1 new_vir = new_tr1->get_vir_tr1();

      mpqc::detail::parallel_print_range_info(world, new_occ, "CCSD-PNO Occ");
      mpqc::detail::parallel_print_range_info(world, new_vir, "CCSD-PNO Vir");

      this->trange1_engine_ = new_tr1;

      TA::DistArray<Tile, Policy> occ_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_occ, new_occ, 1.0);

      TA::DistArray<Tile, Policy> vir_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_vir, new_vir, 1.0);

      // get occupied and virtual orbitals
      auto occ_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"i"));
      auto vir_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"a"));

      auto new_occ_space = occ_space;
      new_occ_space("k,i") = occ_space("k,j") * occ_convert("j,i");

      auto new_vir_space = vir_space;
      new_vir_space("k,a") = vir_space("k,b") * vir_convert("b,a");

      lcao_factory.orbital_space().clear();
      lcao_factory.orbital_space().add(new_occ_space);
      lcao_factory.orbital_space().add(new_vir_space);

      // obtain t2 with new blocking structure
      t2_mp2_("a,b,i,j") = t2_mp2_("c,d,k,l") * vir_convert("c,a") * vir_convert("d,b") *
                      occ_convert("k,i") * occ_convert("l,j");

    }

    template<typename Tile, typename Policy>
    void CCSD_PNO<Tile, Policy>::pno_decom() {

      integer vir = this->trange1_engine_->get_vir();
      double threshold = tcut_;
      ExEnv::out0() << "tcut is " << tcut_ << std::endl;

      std::size_t orb_i = 1, orb_j = 1;


      // compute t_ab = d_a{c} t_{c}{d} d_{d}{b} ({}: pno indices)
      auto transfrom = [&](Tile &in_tile) -> float {

        std::stringstream ss;
//        ss << "input tile: " << in_tile << std::endl;

        // copy in_tile for test purpose
        Tile test_tile = in_tile.clone();

        // get the dimensions of the input tile
        const auto extent = in_tile.range().extent();
        const auto a = extent[0];
        const auto b = extent[1];
        const auto x = std::min(a,b);

        // allocate memory for SVD
        std::unique_ptr<double[]> u_ax(new double[a * x]);
        std::unique_ptr<double[]> v_xb(new double[b * x]);
        std::unique_ptr<double[]> s(new double[x]);

        // SVD of tile: (t_ab)^T
        tensor::algebra::svd(in_tile.data(), s.get(), u_ax.get(), v_xb.get(), a, b, 'A');
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
          // check s for the truncation threshold
          if(s.get()[rank] < threshold)
            break;

          // fold s into u_ab:
          // u_{a} = u_{a} * s_{a}
          madness::cblas::scal(a, s.get()[rank], u_ax.get() + (a * rank), 1);
        }
        ss << "pair (" << orb_i << "," << orb_j << ")  rank: " << rank << std::endl;
        ++orb_i;
        ++orb_j;

        // compute output tile: (t_ab)^T = u_a{c} * v_{c}b
        madness::cblas::gemm(madness::cblas::NoTrans, madness::cblas::NoTrans,
            a, b, rank, 1.0, u_ax.get(), a, v_xb.get(), b, 0.0, in_tile.data(), a);
//        ss << "output tile: " << in_tile << std::endl;

        // test:
        // when threshold is small enough, input and output tiles should be the same
        madness::cblas::gemm(madness::cblas::NoTrans, madness::cblas::NoTrans,
            a, b, rank, 1.0, u_ax.get(), a, v_xb.get(), b, -1.0, test_tile.data(), a);
        ss << "(t_new - t_old).norm(): " << test_tile.norm() << std::endl;

        std::printf("%s", ss.str().c_str());

        return in_tile.norm();
      };

      TA::foreach_inplace(t2_mp2_, transfrom);
    }

    template<typename Tile, typename Policy>
    double CCSD_PNO<Tile, Policy>::value()  {

      auto &world = this->wfn_world()->world();

      double time;
      auto time0 = mpqc::fenced_now(world);

      double ref_energy = ref_wfn_->value();

      auto time1 = mpqc::fenced_now(world);
      time = mpqc::duration_in_s(time0, time1);
      utility::print_par(world, "Total Ref Time: ", time, " S \n");

      // initialize
      init();

      // compute MP2 amplitudes
      ExEnv::out0() << std::endl << "Computing MP2 amplitudes" << std::endl;
      compute_mp2_t2();
      //ExEnv::out0() << indent << "MP2 amplitudes: " << t2_mp2_ << std::endl;

      // reblock MP2 amplitudes
      ExEnv::out0() << std::endl << "Reblocking MP2 amplitudes" << std::endl;
      reblock();
      //ExEnv::out0() << indent << "MP2 amplitudes after reblocking: " << t2_mp2_ << std::endl;

      // pno decomposition
      ExEnv::out0() << std::endl << "Doing PNO decomposition (SVD of MP2 amplitudes)" << std::endl;
      pno_decom();
      double energy = ref_energy;

      return energy;
    }

} // namespace lcao
} // namespace mpqc




#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_
