#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/integrals/periodic_lcao_factory.h"
#include "mpqc/chemistry/qc/mbpt/mp2.h"
#include "mpqc/chemistry/qc/scf/zrhf.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class GammaPointMP2 : public PeriodicLCAOWavefunction<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

  GammaPointMP2() = default;

  GammaPointMP2(const KeyVal &kv) : PeriodicLCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.keyval("ref").class_ptr<PeriodicAOWavefunction<Tile,Policy>>();
    } else {
      throw std::invalid_argument(
          "Default ref wfn in GammaPointMP2 is not supported!");
    }
  }

  ~GammaPointMP2() = default;

  void compute(PropertyBase *pb) override {
    throw std::runtime_error("Not Implemented!!");
  }

  void obsolete() override {
    e_mp2_ = 0.0;
    PeriodicLCAOWavefunction<Tile, Policy>::obsolete();
    ref_wfn_->obsolete();
  }

  double value() override {
    if (this->energy_ == 0.0) {
      double ref_energy = ref_wfn_->value();

      init();

      e_mp2_ = compute_gamma_point_mp2();

      this->energy_ = ref_energy + e_mp2_;
    }
    return this->energy_;
  }

 private:
  std::shared_ptr<PeriodicAOWavefunction<Tile,Policy>> ref_wfn_;
  double e_ref_;
  double e_mp2_;

 private:
  /*!
   * \brief This initializes gamma-point MP2
   */
  void init() {
    auto unitcell = this->lcao_factory().pao_factory().molecule();
    auto co_coeff = ref_wfn_->co_coeff();

    mo_insert_gamma_point(this->lcao_factory(), co_coeff,
                          ref_wfn_->nk(), unitcell, this->occ_block(),
                          this->unocc_block());
  }

  /*!
   * \brief This computes gamma-point MP2 energy
   * \return gamma-point MP2 energy
   */
  double compute_gamma_point_mp2() {
    ExEnv::out0() << "Computing conventional gamma-point MP2 ..." << std::endl;

    auto g_abij = this->lcao_factory().compute(L"<a b|G|i j>");

    ExEnv::out0() << "\n<a b |G|i j> :\n" << g_abij << std::endl;

    return 0.0;
  }
};

}  // namespace lcao

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
