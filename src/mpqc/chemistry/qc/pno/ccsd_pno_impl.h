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
        g_abij = lcao_factory.compute(df ? L"<a b|G|i j>[df]" : L"<a b|G|i j>");

        auto n_occ = trange1_engine->get_occ();
        auto n_frozen = trange1_engine->get_nfrozen();

        // compute t2 amplitudes
        TA::DistArray<Tile, Policy> t2_ijab = d_abij(g_abij, *orbital_energy, n_occ, n_frozen);

        // test mp2 energy:
        double mp2 =  TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), t2_ijab("a,b,i,j"));
        ExEnv::out0() << indent << "Test: MP2 Energy: " << mp2 << std::endl;

        return t2_ijab;
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
     * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
     * | method | string | standard | method to compute ccsd (standard, df,
     * direct) |
     * | converge | double | 1.0e-07 | converge limit |
     * | max_iter | int | 20 | maxmium iteration in CCSD |
     * | print_detail | bool | false | if print more information in CCSD iteration
     * |
     * | less_memory | bool | false | avoid store another abcd term in standard
     * and df method |
     */
    template<typename Tile, typename Policy>
    CCSD_PNO<Tile, Policy>::CCSD_PNO(const KeyVal &kv) : LCAOWavefunction<Tile, Policy>(kv) {
      if (kv.exists("ref")) {
        ref_wfn_ = kv.keyval("ref").class_ptr<Wavefunction>();
      } else {
        throw std::invalid_argument("Default Ref Wfn in CCSD is not support! \n");
      }

      df_ = kv.value<bool>("df", false);
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
    TA::DistArray<Tile, Policy> CCSD_PNO<Tile, Policy>::compute_mp2_t2() {

      auto &lcao_factory = this->lcao_factory();
      return detail::compute_mp2_t2(lcao_factory, this->orbital_energy(),
                                    this->trange1_engine(), false);
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

      TA::DistArray<Tile, Policy> t2_ijab = compute_mp2_t2();

    }

} // namespace lcao
} // namespace mpqc




#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_IMPL_H_
