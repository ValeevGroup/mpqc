//
// Created by Chong Peng on 12/6/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_

namespace mpqc {
namespace lcao {

namespace detail {
template <typename Tile, typename Policy>
double compute_mp2(lcao::LCAOFactory<Tile, Policy> &lcao_factory,
                   std::shared_ptr<Eigen::VectorXd> orbital_energy,
                   std::shared_ptr<mpqc::TRange1Engine> tr1_engine, bool df) {
  auto& world = lcao_factory.world();
  TA::DistArray<Tile, Policy> g_ijab;
  g_ijab = lcao_factory.compute(df ? L"<i j|G|a b>[df]" : L"<i j|G|a b>");
  // compute mp2 energy
  double energy_mp2 =
      (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
          .reduce(mbpt::detail::Mp2Energy<Tile>(
              orbital_energy, tr1_engine->get_occ(),
              tr1_engine->get_nfrozen()));

  utility::print_par(world, (df ? "RI-" : ""), "MP2 Energy: ", energy_mp2, "\n");
  return energy_mp2;
}

} // end of namespace detail

template<typename Tile, typename Policy>
RMP2<Tile,Policy>::RMP2(const KeyVal &kv) : LCAOWavefunction<Tile,Policy>(kv) {
  if (kv.exists("ref")) {
    ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
  } else {
    throw InputError("Default RefWavefunction in RMP2 is not support! \n", __FILE__,__LINE__,"ref");
  }
}

template<typename Tile, typename Policy>
void RMP2<Tile,Policy>::obsolete() {
  mp2_corr_energy_ = 0.0;
  LCAOWavefunction<Tile, Policy>::obsolete();
  ref_wfn_->obsolete();
}

template<typename Tile, typename Policy>
bool RMP2<Tile,Policy>::can_evaluate(Energy* energy) {
  // can only evaluate the energy
  return energy->order() == 0;
};

template<typename Tile, typename Policy>
void RMP2<Tile,Policy>::evaluate(Energy* result) {
  if(!this->computed()){
    /// cast ref_wfn to Energy::Evaluator
    auto ref_evaluator = std::dynamic_pointer_cast<typename Energy::Evaluator>(ref_wfn_);
    if(ref_evaluator == nullptr) {
      std::ostringstream oss;
      oss << "RefWavefunction in CCSD" << ref_wfn_->class_key()
          << " cannot compute Energy" << std::endl;
      throw InputError(oss.str().c_str(), __FILE__, __LINE__);
    }

    ref_evaluator->evaluate(result);

    double ref_energy = this->get_value(result).derivs(0)[0];

    // initialize
    init();

    // compute
    double mp2_corr_energy_ = compute();

    this->computed_ = true;
    this->set_value(result, ref_energy + mp2_corr_energy_);
  }
};


template<typename Tile, typename Policy>
void RMP2<Tile,Policy>::init() {
  if (this->trange1_engine_ == nullptr || this->orbital_energy_ == nullptr) {
    auto mol = this->lcao_factory().ao_factory().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
        this->lcao_factory(), orbital_energy, mol, this->is_frozen_core(),
        this->occ_block(), this->unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
  }
}

template<typename Tile, typename Policy>
double RMP2<Tile,Policy>::compute() {
  auto &lcao_factory = this->lcao_factory();
  return detail::compute_mp2(lcao_factory, this->orbital_energy(),
                             this->trange1_engine(), false);
}


template<typename Tile, typename Policy>
const std::shared_ptr<Wavefunction> RMP2<Tile,Policy>::refwfn() const {
  return ref_wfn_;
}

//
// Member function of RIRMP2 class
//

template<typename Tile, typename Policy>
RIRMP2<Tile,Policy>::RIRMP2(const KeyVal &kv) : RMP2<Tile,Policy>(kv) {}

template<typename Tile, typename Policy>
double RIRMP2<Tile,Policy>::compute() {
  auto &lcao_factory = this->lcao_factory();
  return detail::compute_mp2(lcao_factory, this->orbital_energy(),
                             this->trange1_engine(), true);
}
}  // end of namespace lcao
}  // end of namespace mpqc

#endif //SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_
