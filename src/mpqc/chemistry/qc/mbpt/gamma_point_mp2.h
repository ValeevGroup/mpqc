#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_

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

/*!
 * \brief This computes D_ai = 1 / (e_i - e_a)
 */
template <typename Tile, typename Policy>
TA::DistArray<
    Tile, typename std::enable_if<std::is_same<Policy, TA::SparsePolicy>::value,
                                  TA::SparsePolicy>::type>
d_ai(madness::World &world, const TA::TiledRange &trange, const Vectorz &ens,
     std::size_t n_occ, std::size_t n_frozen) {
  typedef typename TA::DistArray<Tile, Policy>::range_type range_type;

  auto make_tile = [&ens, n_occ, n_frozen](range_type &range, std::size_t ord, Tile *out_tile, TA::TensorF *norms) {

    auto result_tile = Tile(range);

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];

    auto tile_idx = 0;
    typename Tile::value_type tmp = 1.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto i = i0; i < in; ++i, ++tile_idx) {
        const auto e_i = ens[i + n_frozen];
        const auto e_ia = e_i - e_a;
        const auto result_ai = tmp / e_ia;
        result_tile[tile_idx] = result_ai;
      }
    }

    const auto tile_volume = result_tile.range().volume();
    const auto tile_norm = result_tile.norm();
    bool save_norm =
            tile_norm >= tile_volume * TA::SparseShape<float>::threshold();
    if (save_norm) {
        *out_tile = result_tile;
        (*norms)[ord] = tile_norm;
    }
  };

  const auto tvolume = trange.tiles_range().volume();
  std::vector<Tile> tiles(tvolume);
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);
  auto pmap = TA::SparsePolicy::default_pmap(world, tvolume);

  for (auto const ord : *pmap) {
    world.taskq.add(make_tile, trange.make_tile_range(ord), ord, &tiles[ord],
                    &tile_norms);
  }

  world.gop.fence();

  TA::SparseShape<float> shape(world, tile_norms, trange);
  TA::DistArray<Tile, Policy> result(world, trange, shape, pmap);

  for (auto const ord : *pmap) {
    if (result.is_local(ord) && !result.is_zero(ord)) {
      auto &tile = tiles[ord];
      assert(!tile.empty());
      result.set(ord, tile);
    }
  }

  world.gop.fence();
  result.truncate();
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

    print_detail_ = kv.value<bool>("print_detail", false);
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

      auto time0 = mpqc::now(this->wfn_world()->world(), false);
      init();

      e_mp2_ = compute_gamma_point_mp2();

      ExEnv::out0() << "\nGamma-Point MP2 Energy = " << e_mp2_
                    << std::endl;

      this->energy_ = ref_energy + e_mp2_;
      auto time1 = mpqc::now(this->wfn_world()->world(), false);
      auto duration = mpqc::duration_in_s(time0, time1);
      if (print_detail_) {
          ExEnv::out0() << " Total time for gamma-point mp2: " << duration << " s\n";
      }
    }
    return this->energy_;
  }

 private:
  std::shared_ptr<PeriodicAOWavefunction<Tile, Policy>> ref_wfn_;
  double e_ref_;
  double e_mp2_;
  Matrixz C_;
  Vectorz eps_;

  bool print_detail_;

 private:
  /*!
   * \brief This initializes gamma-point MP2
   */
  void init() {
    auto nk = ref_wfn_->nk();
    auto k_size = 1 + detail::k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
    auto unitcell = this->lcao_factory().pao_factory().unitcell();

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
    auto &world = this->wfn_world()->world();

    ExEnv::out0() << "Computing conventional gamma-point MP2 ..." << std::endl;

    auto g_abij = this->lcao_factory().compute(L"<a b|G|i j>");

    auto n_occ = this->lcao_factory().pao_factory().molecule().occupation() / 2;
    std::size_t n_frozen = 0;
    auto t2 = detail::d_abij<Tile, Policy>(g_abij, eps_, n_occ, n_frozen);

    auto time0 = mpqc::now(world, false);

    std::complex<double> e_complex = TA::dot(
        2.0 * g_abij("a, b, i, j") - g_abij("b, a, i, j"), t2("a, b, i, j"));
    double e_mp2 = e_complex.real();

    auto time1 = mpqc::now(world, false);
    auto duration = mpqc::duration_in_s(time0, time1);

    if (print_detail_) {
        ExEnv::out0() << " Time for energy computation (CO transformation not included): " << duration << " s\n";
    }
    return e_mp2;
  }
};

#if TA_DEFAULT_POLICY == 0
extern template class GammaPointMP2<TA::TensorZ, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class GammaPointMP2<TA::TensorZ, TA::SparsePolicy>;
#endif


}  // namespace lcao

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
