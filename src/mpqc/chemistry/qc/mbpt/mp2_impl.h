//
// Created by Chong Peng on 12/6/16.
//

#ifndef MPQC_MP2_IMPL_H
#define MPQC_MP2_IMPL_H

namespace mpqc {
namespace mbpt {

namespace detail {
template <typename Tile, typename Policy>
double compute_mp2(integrals::LCAOFactory<Tile, Policy> &lcao_factory,
                   std::shared_ptr<Eigen::VectorXd> orbital_energy,
                   std::shared_ptr<mpqc::TRange1Engine> tr1_engine, bool df) {
  auto& world = lcao_factory.world();
  TA::DistArray<Tile, Policy> g_ijab;
  g_ijab = lcao_factory.compute(df ? L"<i j|G|a b>[df]" : L"<i j|G|a b>");
  // compute mp2 energy
  double energy_mp2 =
      (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
          .reduce(detail::Mp2Energy<Tile>(
              orbital_energy, tr1_engine->get_occ(),
              tr1_engine->get_nfrozen()));

  utility::print_par(world, (df ? "RI-" : ""), "MP2 Energy: ", energy_mp2, "\n");
  return energy_mp2;
}

} // end of namespace detail

template<typename Tile, typename Policy>
RMP2<Tile,Policy>::RMP2(const KeyVal &kv) : qc::LCAOWavefunction<Tile,Policy>(kv) {
  if (kv.exists("ref")) {
    ref_wfn_ = kv.keyval("ref").class_ptr<qc::Wavefunction>();
  } else {
    throw std::invalid_argument("Default Ref Wfn in RMP2 is not support! \n");
  }
}

template<typename Tile, typename Policy>
RMP2<Tile,Policy>::~RMP2() = default;

template<typename Tile, typename Policy>
void RMP2<Tile,Policy>::obsolete() {
  this->energy_ = 0.0;
  qc::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>::obsolete();
  ref_wfn_->obsolete();
}

template<typename Tile, typename Policy>
double RMP2<Tile,Policy>::value() {
  if (this->energy_ == 0.0) {
    auto &world = this->wfn_world()->world();

    double time;
    auto time0 = mpqc::fenced_now(world);

    double ref_energy = ref_wfn_->value();

    auto time1 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "Total Ref Time: ", time, " S \n");

    // initialize
    init();

    double mp2_energy = compute();

    this->energy_ = mp2_energy + ref_energy;

    auto time2 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time1, time2);
    utility::print_par(world, "Total MP2 Correlation Time: ", time, " S \n");

    time = mpqc::duration_in_s(time0, time2);
    utility::print_par(world, "Total MP2 Time: ", time, " S \n");
  }
  return this->energy_;
}

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
void RMP2<Tile,Policy>::compute(qc::PropertyBase *pb) {}

template<typename Tile, typename Policy>
const std::shared_ptr<qc::Wavefunction> RMP2<Tile,Policy>::refwfn() const {
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
}  // end of namespace mbpt
}  // end of namespace mpqc

#endif //MPQC_MP2_IMPL_H
