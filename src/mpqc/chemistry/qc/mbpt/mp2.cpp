//
// Created by Chong Peng on 10/6/16.
//

#include "mp2.h"

MPQC_CLASS_EXPORT_KEY2(mpqc::mbpt::RMP2, "RMP2");
MPQC_CLASS_EXPORT_KEY2(mpqc::mbpt::RIRMP2, "RI-RMP2");

namespace mpqc{
namespace mbpt {

RMP2::RMP2(const KeyVal &kv) : LCAOWavefunction(kv) {

  rmp2_energy_ = 0.0;

  if (kv.exists("ref")) {
    ref_wfn_ = kv.keyval("ref").class_ptr<qc::Wavefunction>();
  } else {
    throw std::invalid_argument("Default Ref Wfn in RMP2 is not support! \n");
  }

}

double RMP2::value() {
  if (rmp2_energy_ == 0.0) {
    double ref_energy = ref_wfn_->value();

    auto mol = this->lcao_factory().atomic_integral().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(), occ_block(), unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    double mp2_energy = compute();

    rmp2_energy_ = mp2_energy + ref_energy;
  }
  return rmp2_energy_;
}

double RMP2::compute() {
  auto &lcao_factory = this->lcao_factory();
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>");
  // compute mp2 energy
  double energy_mp2 = (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
      .reduce(detail::Mp2Energy<TA::TensorD>(orbital_energy_, trange1_engine_->get_occ(),
                                             trange1_engine_->get_nfrozen()));

  if (g_ijab.get_world().rank() == 0) {
    std::cout << "MP2 Energy  " << energy_mp2 << std::endl;
  }
  return energy_mp2;
}

void RMP2::compute(qc::PropertyBase *pb) {}

RIRMP2::RIRMP2(const KeyVal &kv) : RMP2(kv) {}

double RIRMP2::compute() {

  auto &lcao_factory = this->lcao_factory();

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>[df]");
  // compute mp2 energy
  double energy_mp2 =
      (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
          .reduce(detail::Mp2Energy<TA::TensorD>(orbital_energy_, trange1_engine_->get_occ(),
                                                 trange1_engine_->get_nfrozen()));

  if (g_ijab.get_world().rank() == 0) {
    std::cout << "MP2 Energy With DF: " << energy_mp2 << std::endl;
  }

  return energy_mp2;
}
} // end of namespace mbpt
} // end of namespace mpqc
