//
// Created by Chong Peng on 6/24/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_H_

#include "mpqc/mpqc_config.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/chemistry/qc/lcao/scf/mo_build.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/util/external/madworld/parallel_print.h"

using namespace mpqc;

namespace mpqc {
namespace lcao {

namespace mbpt {
namespace detail {

template <typename Tile>
struct Mp2Energy {
  using result_type = double;
  using argument_type = Tile;

  std::shared_ptr<Eigen::VectorXd> vec_;
  std::size_t n_occ_;
  std::size_t n_frozen_;

  Mp2Energy(std::shared_ptr<Eigen::VectorXd> vec, std::size_t n_occ,
            std::size_t n_frozen)
      : vec_(vec), n_occ_(n_occ), n_frozen_(n_frozen) {}

  Mp2Energy(Mp2Energy const &) = default;

  result_type operator()() const { return 0.0; }

  result_type operator()(result_type const &t) const { return t; }

  void operator()(result_type &me, result_type const &other) const {
    me += other;
  }

  void operator()(result_type &me, argument_type const &tile) const {
    auto const &range = tile.range();
    auto const &vec = *vec_;
    auto const st = range.lobound_data();
    auto const fn = range.upbound_data();
    auto tile_idx = 0;

    auto sti = st[0];
    auto fni = fn[0];
    auto stj = st[1];
    auto fnj = fn[1];
    auto sta = st[2];
    auto fna = fn[2];
    auto stb = st[3];
    auto fnb = fn[3];

    for (auto i = sti; i < fni; ++i) {
      const auto e_i = vec[i + n_frozen_];
      for (auto j = stj; j < fnj; ++j) {
        const auto e_ij = e_i + vec[j + n_frozen_];
        for (auto a = sta; a < fna; ++a) {
          const auto e_ija = e_ij - vec[a + n_occ_];
          for (auto b = stb; b < fnb; ++b, ++tile_idx) {
            const auto e_iajb = e_ija - vec[b + n_occ_];
            me += 1.0 / (e_iajb)*tile.data()[tile_idx];
          }
        }
      }
    }
  }
};

template <typename Tile, typename Policy>
double compute_mp2(lcao::LCAOFactory<Tile, Policy> &lcao_factory,
                   std::shared_ptr<Eigen::VectorXd> orbital_energy,
                   std::shared_ptr<mpqc::TRange1Engine> tr1_engine, bool df);

}  // end of namespce detail
}  // end of namespce mbpt

/**
 *  \brief MP2 class for closed-shell system
 *
 *  KeyVal type of this class RMP2
 */

template<typename Tile, typename Policy>
class RMP2 : public lcao::LCAOWavefunction<Tile,Policy>, public CanEvaluate<Energy> {
 public:

  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords: takes all keywords from LCAOWavefunction
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
   */
  // clang-format on
  RMP2(const KeyVal &kv);
  virtual ~RMP2() = default;

  void obsolete() override;

  const std::shared_ptr<lcao::Wavefunction> refwfn() const;

 protected:

  bool can_evaluate(Energy* energy) override;

  void evaluate(Energy* result) override;

  /// function to compute mp2 correlation energy
  virtual double compute();

  /// initialize orbitals
  virtual void init();
  std::shared_ptr<lcao::Wavefunction> ref_wfn_;

 private:
   double mp2_corr_energy_;
};

/**
 *  \brief MP2 class for closed-shell system with density fitting
 *
 *  KeyVal type of this class RI-RMP2
 */

template <typename Tile, typename Policy>
class RIRMP2 : public RMP2<Tile, Policy> {
 public:
  /**
  * KeyVal constructor
  * @param kv
  *
  * keywords: inherit all keywords from RMP2
  */
  RIRMP2(const KeyVal &kv);
  ~RIRMP2() = default;

 protected:
  /// override the compute function from RMP2
  double compute() override;
};

#if TA_DEFAULT_POLICY == 0
extern template class RMP2<TA::TensorD, TA::DensePolicy>;
extern template class RIRMP2<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class RMP2<TA::TensorD, TA::SparsePolicy>;
extern template class RIRMP2<TA::TensorD, TA::SparsePolicy>;
#endif
}  // namespace lcao
}  // namespace mpqc

#include "mp2_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_MP2_H_
