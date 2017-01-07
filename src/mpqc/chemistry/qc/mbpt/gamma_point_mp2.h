#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/integrals/periodic_lcao_factory.h"
#include "mpqc/chemistry/qc/mbpt/mp2.h"
#include "mpqc/chemistry/qc/scf/zrhf.h"

namespace mpqc {
namespace lcao {

namespace detail {

/*!
 * \brief This computes D_abij = g_abij / (e_i + e_j - e_a - a_b)
 */
template <typename Tile, typename Policy>
TA::Array<std::complex<double>, 4, Tile, Policy> d_abij(
    TA::Array<std::complex<double>, 4, Tile, Policy> &abij, const Vectorz &ens,
    std::size_t n_occ, std::size_t n_frozen) {
  auto convert = [&ens, n_occ, n_frozen](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    float norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens[j + n_frozen];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto old = arg_tile[tile_idx];
            const auto result_abij = old / (e_iajb);
            norm += std::abs(result_abij) * std::abs(result_abij);
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach (abij, convert);
  abij.world().gop.fence();
  return result;
}

}  // namespace detail

template <typename Tile, typename Policy>
class GammaPointMP2 : public PeriodicLCAOWavefunction<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

  GammaPointMP2() = default;

  GammaPointMP2(const KeyVal &kv) : PeriodicLCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ =
          kv.keyval("ref").class_ptr<PeriodicAOWavefunction<Tile, Policy>>();
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

      ExEnv::out0() << "\nTotal gamma-point MP2 energy = " << e_mp2_
                    << std::endl;

      this->energy_ = ref_energy + e_mp2_;
    }
    return this->energy_;
  }

 private:
  std::shared_ptr<PeriodicAOWavefunction<Tile, Policy>> ref_wfn_;
  double e_ref_;
  double e_mp2_;
  Matrixz C_;
  Vectorz eps_;

 private:
  /*!
   * \brief This initializes gamma-point MP2
   */
  void init() {
    auto nk = ref_wfn_->nk();
    auto k_size = 1 + detail::k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
    auto unitcell = this->lcao_factory().pao_factory().molecule();

    int64_t gamma_point;
    if (k_size % 2 == 1)
      gamma_point = (k_size - 1) / 2;
    else
      throw std::invalid_argument(
          "# of k points must be odd in order to run gamma-point methods");

    C_ = ref_wfn_->co_coeff()[gamma_point];
    eps_ = ref_wfn_->co_energy()[gamma_point];

    mo_insert_gamma_point(this->lcao_factory(), C_, unitcell, this->occ_block(),
                          this->unocc_block());
  }

  /*!
   * \brief This computes gamma-point MP2 energy
   * \return gamma-point MP2 energy
   */
  double compute_gamma_point_mp2() {
    ExEnv::out0() << "Computing conventional gamma-point MP2 ..." << std::endl;

    auto g_abij = this->lcao_factory().compute(L"<a b|G|i j>");

    auto n_occ = this->lcao_factory().pao_factory().molecule().occupation() / 2;
    std::size_t n_frozen = 0;
    auto t2 = detail::d_abij<Tile, Policy>(g_abij, eps_, n_occ, n_frozen);

    std::complex<double> e_complex = TA::dot(
        2.0 * g_abij("a, b, i, j") - g_abij("b, a, i, j"), t2("a, b, i, j"));
    double e_mp2 = e_complex.real();

    return e_mp2;
  }
};

}  // namespace lcao

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
