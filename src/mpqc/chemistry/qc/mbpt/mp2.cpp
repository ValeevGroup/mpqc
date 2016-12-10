//
// Created by Chong Peng on 10/6/16.
//

#include "mp2.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("RMP2", mpqc::lcao::RMP2);
MPQC_CLASS_EXPORT2("RI-RMP2", mpqc::lcao::RIRMP2);

namespace mpqc {
namespace lcao {

RMP2::RMP2(const KeyVal &kv) : LCAOWavefunction(kv) {
  if (kv.exists("ref")) {
    ref_wfn_ = kv.keyval("ref").class_ptr<lcao::Wavefunction>();
  } else {
    throw std::invalid_argument("Default Ref Wfn in RMP2 is not support! \n");
  }
}

RMP2::~RMP2() = default;

void RMP2::obsolete() {
  this->energy_ = 0.0;
  lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>::obsolete();
  ref_wfn_->obsolete();
}

double RMP2::value() {
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

void RMP2::init() {
  if (this->trange1_engine_ == nullptr || this->orbital_energy_ == nullptr) {
    auto mol = this->lcao_factory().ao_factory().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(),
        occ_block(), unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
  }
}

double RMP2::compute() {
  auto &lcao_factory = this->lcao_factory();
  return mbpt::detail::compute_mp2(lcao_factory, this->orbital_energy(),
                                   this->trange1_engine(), false);
}

void RMP2::compute(lcao::PropertyBase *pb) {}

const std::shared_ptr<lcao::Wavefunction> RMP2::refwfn() const {
  return ref_wfn_;
}

//
// Member function of RIRMP2 class
//

RIRMP2::RIRMP2(const KeyVal &kv) : RMP2(kv) {}

double RIRMP2::compute() {
  auto &lcao_factory = this->lcao_factory();
  return mbpt::detail::compute_mp2(lcao_factory, this->orbital_energy(),
                                   this->trange1_engine(), true);
}
}  // namespace lcao
}  // namespace mpqc
