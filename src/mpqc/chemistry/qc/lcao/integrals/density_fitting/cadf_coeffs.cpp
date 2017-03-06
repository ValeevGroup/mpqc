/*
 * Created by Drew Lewis 02/23/2017
 */

#include "cadf_coeffs.h"

namespace mpqc {
namespace lcao {
namespace detail {
gaussian::Basis by_center_basis(gaussian::Basis const &in) {
  auto const shells = in.flattened_shells();
  auto it = shells.begin();
  auto end = shells.end();

  std::vector<gaussian::ShellVec> out;
  while (it != end) {
    auto center = it->O;
    auto next_center = center;

    gaussian::ShellVec atom;
    while (it != end && next_center == center) {
      atom.push_back(*it);
      ++it;
      next_center = it->O;
    }

    out.emplace_back(std::move(atom));
  }

  return gaussian::Basis(std::move(out));
}

TA::TiledRange cadf_trange(gaussian::Basis const &obs_by_atom,
                           gaussian::Basis const &dfbs_by_atom) {
  std::vector<TA::TiledRange1> trange1s;
  trange1s.emplace_back(dfbs_by_atom.create_trange1());
  trange1s.emplace_back(obs_by_atom.create_trange1());
  trange1s.emplace_back(obs_by_atom.create_trange1());

  return TA::TiledRange(trange1s.begin(), trange1s.end());
}

TA::SparseShape<float> cadf_shape(madness::World &world,
                                  TA::TiledRange const &trange) {
  auto &tiles = trange.tiles_range();
  TA::Tensor<float> norms(tiles, 0.0);
  auto natoms = tiles.extent_data()[0];  // trange is tiled by atom

  for (auto i = 0ul; i < natoms; ++i) {
    for (auto j = 0ul; j < natoms; ++j) {
      norms(i, i, j) = std::numeric_limits<float>::max();
      norms(j, i, j) = std::numeric_limits<float>::max();
    }
  }

  return TA::SparseShape<float>(world, norms, trange);
}

TA::SparseShape<float> cadf_shape_cluster(
    TA::DistArray<TA::Tensor<double>, TA::SparsePolicy> const &C_atom,
    TA::TiledRange const &trange,
    std::unordered_map<int64_t, std::vector<int64_t>> &c2a  // cluster to atom
    ) {
  auto const &tiles = trange.tiles_range();
  auto norms = TA::Tensor<float>(tiles, 0.0);

  auto const &trange_atom = C_atom.trange();
  const auto volume_atom = trange_atom.tiles_range().volume();

  // Loop over all tiles of C_atom to figure which tiles of C_cluster they got
  // to
  for (auto ord_atom = 0ul; ord_atom < volume_atom; ++ord_atom) {
    if (C_atom.is_zero(ord_atom)) {
      continue;
    }
    TA::Range range_atom = trange_atom.make_tile_range(ord_atom);
    auto lo_atom = range_atom.lobound_data();
    auto ord_cluster = tiles.ordinal(trange.element_to_tile(lo_atom));

    norms[ord_cluster] = std::numeric_limits<float>::max();

    // Write the atom to the cluster
    auto it = c2a.find(int64_t(ord_cluster));
    if (it != c2a.end()) {
      it->second.push_back(ord_atom);
    } else {
      c2a[int64_t(ord_cluster)] = std::vector<int64_t>{int64_t(ord_atom)};
    }
  }

  // No world needed every node computed entire shape
  return TA::SparseShape<float>(norms, trange);
}
//
// Function to compute the by atom screener for cadf eri3 integrals
std::shared_ptr<gaussian::SchwarzScreen> cadf_by_atom_screener(
    madness::World &world, gaussian::Basis const &obs,
    gaussian::Basis const &dfbs, double threshold) {
  const auto three_ref_array = utility::make_array_of_refs(dfbs, obs, obs);

  const auto three_array = gaussian::BasisVector{{dfbs, obs, obs}};

  auto eng4 = make_engine_pool(libint2::Operator::coulomb, three_ref_array,
                               libint2::BraKet::xx_xx);

  return std::make_shared<gaussian::SchwarzScreen>(
      gaussian::create_schwarz_screener(world, eng4, three_array, threshold));
}

TA::DistArray<TA::Tensor<double>, TA::SparsePolicy> cadf_by_atom_coeffs(
    madness::World &world, gaussian::Basis const &by_cluster_obs,
    gaussian::Basis const &by_cluster_dfbs) {
  auto obs = detail::by_center_basis(by_cluster_obs);
  auto dfbs = detail::by_center_basis(by_cluster_dfbs);

  const auto dfbs_ref_array = utility::make_array_of_refs(dfbs, dfbs);
  const auto dfbs_array = gaussian::BasisVector{{dfbs, dfbs}};
  auto eng2 = make_engine_pool(libint2::Operator::coulomb, dfbs_ref_array,
                               libint2::BraKet::xs_xs);
  auto M = gaussian::sparse_integrals(world, eng2, dfbs_array);

  auto screener = detail::cadf_by_atom_screener(world, obs, dfbs, 1e-12);
  auto eng3 = make_engine_pool(libint2::Operator::coulomb,
                               utility::make_array_of_refs(dfbs, obs, obs),
                               libint2::BraKet::xs_xx);

  auto eri3_norms = [&](TA::Tensor<float> const &in) {
    const auto thresh = screener->skip_threshold();

    auto const &range = in.range();
    auto ext = range.extent_data();

    std::vector<std::pair<std::array<int,3>, float>> norms;
    norms.reserve(ext[1] * ext[2]);

    using idx_type = std::array<int, 3>;
    auto val_max = std::numeric_limits<float>::max();
    for (auto a = 0; a < ext[1]; ++a) {
      for (auto b = 0; b < ext[2]; ++b) {
        auto in_val = std::max(in(a, a, b), in(b, a, b));

        if (in_val >= thresh) {
          norms.emplace_back(std::make_pair(idx_type{a,a,b}, val_max));
          norms.emplace_back(std::make_pair(idx_type{b,a,b}, val_max));
        }
      }
    }

    return norms;
  };

  auto three_array = std::vector<gaussian::Basis>{dfbs, obs, obs};

  auto trange = gaussian::detail::create_trange(three_array);
  auto pmap =
      TA::SparsePolicy::default_pmap(world, trange.tiles_range().volume());

  bool replicate = true;
  auto norms =
      eri3_norms(screener->norm_estimate(world, three_array, *pmap, replicate));

  auto eri3 = direct_sparse_integrals(world, eng3, three_array, norms,
                                      std::move(screener));

  return cadf_by_atom_array(M, eri3, detail::cadf_trange(obs, dfbs));
}

}  // namespace detail
}  // namespace lcao
}  // namespace mpqc
