//
// Created by Chong Peng on 12/6/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_

namespace mpqc {
namespace lcao {

namespace detail {
template <typename Tile, typename Policy>
double compute_mp2(
    lcao::LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const std::shared_ptr<const Eigen::VectorXd>& orbital_energy,
    const std::shared_ptr<const ::mpqc::utility::TRange1Engine>& tr1_engine,
    bool df) {
  auto& world = lcao_factory.world();
  TA::DistArray<Tile, Policy> g_ijab;
  g_ijab = lcao_factory.compute(df ? L"<i j|G|a b>[df]" : L"<i j|G|a b>");
  // compute mp2 energy
  double energy_mp2 =
      (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
          .reduce(mbpt::detail::Mp2Energy<Tile>(orbital_energy,
                                                tr1_engine->get_occ(),
                                                tr1_engine->get_nfrozen()));

  utility::print_par(world, (df ? "RI-" : ""), "MP2 Energy: ", energy_mp2,
                     "\n");
  return energy_mp2;
}

}  // namespace detail

template <typename Tile, typename Policy>
RMP2<Tile, Policy>::RMP2(const KeyVal& kv)
    : LCAOWavefunction<Tile, Policy>(kv) {
  if (kv.exists("ref")) {
    ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
  } else {
    throw InputError("Default RefWavefunction in RMP2 is not support! \n",
                     __FILE__, __LINE__, "ref");
  }
}

template <typename Tile, typename Policy>
void RMP2<Tile, Policy>::obsolete() {
  mp2_corr_energy_ = 0.0;
  LCAOWavefunction<Tile, Policy>::obsolete();
  ref_wfn_->obsolete();
}

template <typename Tile, typename Policy>
bool RMP2<Tile, Policy>::can_evaluate(Energy* energy) {
  // can only evaluate the energy
  return energy->order() == 0;
};

template <typename Tile, typename Policy>
void RMP2<Tile, Policy>::evaluate(Energy* result) {
  auto target_precision = result->target_precision(0);
  // compute only if never computed, or requested with higher precision than
  // before
  if (!this->computed() || computed_precision_ > target_precision) {
    // compute reference to higher precision than required of correlation energy
    auto target_ref_precision = target_precision / 100.;
    auto ref_energy = std::make_shared<Energy>(ref_wfn_, target_ref_precision);
    ::mpqc::evaluate(*ref_energy, ref_wfn_);

    this->init_sdref(ref_wfn_, target_ref_precision);

    // compute
    mp2_corr_energy_ = compute();

    this->computed_ = true;
    this->set_value(result, ref_energy->energy() + mp2_corr_energy_);
    computed_precision_ = target_precision;
  }
};

template <typename Tile, typename Policy>
double RMP2<Tile, Policy>::compute() {
  return detail::compute_mp2(
      this->lcao_factory(),
      make_diagonal_fpq(this->lcao_factory(), this->ao_factory(),false),
      this->trange1_engine(), false);
}

template <typename Tile, typename Policy>
const std::shared_ptr<Wavefunction> RMP2<Tile, Policy>::refwfn() const {
  return ref_wfn_;
}

//
// Member function of RIRMP2 class
//

template <typename Tile, typename Policy>
RIRMP2<Tile, Policy>::RIRMP2(const KeyVal& kv) : RMP2<Tile, Policy>(kv) {}

template <typename Tile, typename Policy>
double RIRMP2<Tile, Policy>::compute() {
  return detail::compute_mp2(
      this->lcao_factory(),
      make_diagonal_fpq(this->lcao_factory(), this->ao_factory(), true),
      this->trange1_engine(), true);
}
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_IMPL_H_
