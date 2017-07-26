/*
 * Created by Drew Lewis 02/23/2017
 */

#include "cadf_coeffs.h"

#include <fstream>

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
    gaussian::ShellVec atom;
    while (it != end && it->O == center) {
      atom.push_back(*it);
      ++it;
    }

    out.emplace_back(std::move(atom));
  }

  return gaussian::Basis(std::move(out));
}

void print_shape(TA::Tensor<float> const &t, std::string const &file_name) {
  auto rank = t.range().rank();
  auto ext = t.range().extent_data();

  std::ofstream out_file(file_name);

  if (rank == 3) {
    for (auto X = 0ul; X < ext[0]; ++X) {
      auto ord = 0;
      for (auto a = 0ul; a < ext[1]; ++a) {
        for (auto b = 0ul; b < ext[2]; ++b, ++ord) {
          out_file << t(X, a, b);
          if (ord < (ext[1] * ext[2] - 1)) {
            out_file << ", ";
          }
        }
      }
      out_file << "\n";
    }
  } else if (rank == 4) {
    for (auto a = 0ul; a < ext[0]; ++a) {
      for (auto b = 0ul; b < ext[1]; ++b) {
        auto ord = 0;
        for (auto c = 0ul; c < ext[2]; ++c) {
          for (auto d = 0ul; d < ext[3]; ++d, ++ord) {
            out_file << t(a, b, c, d);
            if (ord < (ext[2] * ext[3] - 1)) {
              out_file << ", ";
            }
          }
        }
        out_file << "\n";
      }
    }
  }
}

TA::TiledRange cadf_trange(gaussian::Basis const &obs_by_atom,
                           gaussian::Basis const &dfbs_by_atom) {
  std::vector<TA::TiledRange1> trange1s;
  trange1s.emplace_back(dfbs_by_atom.create_trange1());
  trange1s.emplace_back(obs_by_atom.create_trange1());
  trange1s.emplace_back(obs_by_atom.create_trange1());

  return TA::TiledRange(trange1s.begin(), trange1s.end());
}

TA::TiledRange cadf_trange(gaussian::Basis const &bs0_by_atom,
                           gaussian::Basis const &bs1_by_atom,
                           gaussian::Basis const &dfbs_by_atom) {
  std::vector<TA::TiledRange1> trange1s;
  trange1s.emplace_back(dfbs_by_atom.create_trange1());
  trange1s.emplace_back(bs0_by_atom.create_trange1());
  trange1s.emplace_back(bs1_by_atom.create_trange1());

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

std::shared_ptr<gaussian::SchwarzScreen> cadf_by_atom_screener(
    madness::World &world, gaussian::Basis const &bs0,
    gaussian::Basis const &bs1, gaussian::Basis const &dfbs, double threshold) {
  const auto three_ref_array = utility::make_array_of_refs(dfbs, bs0, bs1);

  const auto three_array = gaussian::BasisVector{{dfbs, bs0, bs1}};

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

    std::vector<std::pair<std::array<int, 3>, float>> norms;
    norms.reserve(ext[1] * ext[2]);

    using idx_type = std::array<int, 3>;
    auto val_max = std::numeric_limits<float>::max();
    for (auto a = 0; a < ext[1]; ++a) {
      for (auto b = 0; b < ext[2]; ++b) {
        auto in_val = std::max(in(a, a, b), in(b, a, b));

        if (in_val >= thresh) {
          norms.emplace_back(std::make_pair(idx_type{{a, a, b}}, val_max));
          norms.emplace_back(std::make_pair(idx_type{{b, a, b}}, val_max));
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

TA::DistArray<TA::Tensor<double>, TA::SparsePolicy> cadf_by_atom_coeffs(
    const TA::DistArray<TA::Tensor<double>, TA::SparsePolicy> &M,
    gaussian::Basis const &by_cluster_bs0,
    gaussian::Basis const &by_cluster_bs1,
    gaussian::Basis const &by_cluster_dfbs, size_t const &natoms_per_uc,
    Vector3i const &lattice_range0,
    Vector3i const &lattice_range1,
    Vector3i const &lattice_range_df,
    Vector3i const &lattice_center0,
    Vector3i const &lattice_center1,
    Vector3i const &lattice_center_df) {
  auto &world = M.world();
//  mpqc::time_point t0, t1;

//  t0 = mpqc::fenced_now(world);
  auto bs0 = detail::by_center_basis(by_cluster_bs0);
  auto bs1 = detail::by_center_basis(by_cluster_bs1);
  auto dfbs = detail::by_center_basis(by_cluster_dfbs);
//  t1 = mpqc::fenced_now(world);
//  double t_by_center_basis = mpqc::duration_in_s(t0, t1);

//  t0 = mpqc::fenced_now(world);
  auto screener = detail::cadf_by_atom_screener(world, bs0, bs1, dfbs, 1e-12);
//  t1 = mpqc::fenced_now(world);
//  double t_cadf_by_atom_screener = mpqc::duration_in_s(t0, t1);

//  t0 = mpqc::fenced_now(world);
  auto eng3 = make_engine_pool(libint2::Operator::coulomb,
                               utility::make_array_of_refs(dfbs, bs0, bs1),
                               libint2::BraKet::xs_xx);

  using ::mpqc::lcao::detail::direct_3D_idx;
  using ::mpqc::lcao::detail::direct_ord_idx;

  auto eri3_norms = [&](TA::Tensor<float> const &in) {
    const auto thresh = screener->skip_threshold();

    auto const &range = in.range();
    auto ext = range.extent_data();

    std::vector<std::pair<std::array<int, 3>, float>> norms;
    norms.reserve(ext[1] * ext[2]);

    using idx_type = std::array<int, 3>;
    auto val_max = std::numeric_limits<float>::max();
    for (auto a = 0; a < ext[1]; ++a) {
      const auto uc0_ord = a / natoms_per_uc;
      const auto uc0_3D = direct_3D_idx(uc0_ord, lattice_range0) + lattice_center0;
      const auto uc0_ord_in_df = direct_ord_idx(uc0_3D - lattice_center_df, lattice_range_df);
      const auto a_in_df = uc0_ord_in_df * natoms_per_uc + a % natoms_per_uc;

      for (auto b = 0; b < ext[2]; ++b) {
        auto uc1_ord = b / natoms_per_uc;
        auto uc1_3D = direct_3D_idx(uc1_ord, lattice_range1) + lattice_center1;
        auto uc1_ord_in_df = direct_ord_idx(uc1_3D - lattice_center_df, lattice_range_df);
        auto b_in_df = uc1_ord_in_df * natoms_per_uc + b % natoms_per_uc;

        auto in_val = std::max(in(a_in_df, a, b), in(b_in_df, a, b));

        if (in_val >= thresh) {
          norms.emplace_back(
              std::make_pair(idx_type{{int(a_in_df), a, b}}, val_max));
          norms.emplace_back(
              std::make_pair(idx_type{{int(b_in_df), a, b}}, val_max));
        }
      }
    }

    return norms;
  };

  auto three_array = std::vector<gaussian::Basis>{dfbs, bs0, bs1};

  auto trange = gaussian::detail::create_trange(three_array);
  auto pmap =
      TA::SparsePolicy::default_pmap(world, trange.tiles_range().volume());

  bool replicate = true;
  auto norms =
      eri3_norms(screener->norm_estimate(world, three_array, *pmap, replicate));
//  t1 = mpqc::fenced_now(world);
//  double t_norms = mpqc::duration_in_s(t0, t1);

//  t0 = mpqc::fenced_now(world);
  auto eri3 = direct_sparse_integrals(world, eng3, three_array, norms,
                                      std::move(screener));
//  t1 = mpqc::fenced_now(world);
//  double t_direct_eri3 = mpqc::duration_in_s(t0, t1);

//  ExEnv::out0() << "  by_center_basis time: " << t_by_center_basis << " s\n"
//                << "  by_center screener:   " << t_cadf_by_atom_screener << " s\n"
//                << "  eri3 norms:           " << t_norms << " s\n"
//                << "  eri3 direct array:    " << t_direct_eri3 << "s\n";


  return cadf_by_atom_array(M, eri3, detail::cadf_trange(bs0, bs1, dfbs), natoms_per_uc,
                            lattice_range0, lattice_range1, lattice_range_df,
                            lattice_center0, lattice_center1, lattice_center_df);
}

}  // namespace detail
}  // namespace lcao
}  // namespace mpqc
