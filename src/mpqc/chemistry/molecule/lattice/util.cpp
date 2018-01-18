#include "mpqc/chemistry/molecule/lattice/util.h"

namespace mpqc {
namespace detail {

Vector3d direct_vector(int64_t ord_idx, Vector3i const &lattice_max,
                       Vector3d const &dcell) {
  auto z = ord_idx % (2 * lattice_max(2) + 1);
  auto y = (ord_idx / (2 * lattice_max(2) + 1)) % (2 * lattice_max(1) + 1);
  auto x = ord_idx / (2 * lattice_max(2) + 1) / (2 * lattice_max(1) + 1);
  Vector3d result((x - lattice_max(0)) * dcell(0),
                  (y - lattice_max(1)) * dcell(1),
                  (z - lattice_max(2)) * dcell(2));
  return result;
}

Vector3d k_vector(int64_t ord_idx, Vector3i const &nk, Vector3d const &dcell) {
  Vector3d result;
  auto x = ord_idx / nk(2) / nk(1);
  auto y = (ord_idx / nk(2)) % nk(1);
  auto z = ord_idx % nk(2);
  result(0) = (dcell(0) == 0.0)
                  ? 0.0
                  : (-1.0 + (2.0 * (x + 1) - 1.0) / nk(0)) * (M_PI / dcell(0));
  result(1) = (dcell(1) == 0.0)
                  ? 0.0
                  : (-1.0 + (2.0 * (y + 1) - 1.0) / nk(1)) * (M_PI / dcell(1));
  result(2) = (dcell(2) == 0.0)
                  ? 0.0
                  : (-1.0 + (2.0 * (z + 1) - 1.0) / nk(2)) * (M_PI / dcell(2));
  return result;
}

int64_t direct_ord_idx(Vector3i const &in_3D_idx, Vector3i const &lattice_max) {
  return direct_ord_idx(in_3D_idx(0), in_3D_idx(1), in_3D_idx(2), lattice_max);
}

int64_t direct_ord_idx(int64_t x, int64_t y, int64_t z,
                       Vector3i const &lattice_max) {
  if (lattice_max(0) >= 0 && lattice_max(1) >= 0 && lattice_max(2) >= 0 &&
      std::abs(x) <= lattice_max(0) && std::abs(y) <= lattice_max(1) &&
      std::abs(z) <= lattice_max(2)) {
    int64_t idx = (x + lattice_max(0)) * (2 * lattice_max(1) + 1) *
                      (2 * lattice_max(2) + 1) +
                  (y + lattice_max(1)) * (2 * lattice_max(2) + 1) +
                  (z + lattice_max(2));
    return idx;
  } else {
    throw "invalid lattice sum index/boundaries";
  }
}

int64_t k_ord_idx(Vector3i const &in_3D_idx, Vector3i const &nk) {
  return k_ord_idx(in_3D_idx(0), in_3D_idx(1), in_3D_idx(2), nk);
}

int64_t k_ord_idx(int64_t x, int64_t y, int64_t z, Vector3i const &nk) {
  if (nk(0) >= 1 && nk(1) >= 1 && nk(2) >= 1 && x >= 0 && y >= 0 && z >= 0 &&
      x < nk(0) && y < nk(1) && z < nk(2)) {
    int64_t idx = x * nk(0) * nk(1) + y * nk(1) + z;
    return idx;
  } else {
    throw "invalid k-space index/boundaries";
  }
}

Vector3i direct_3D_idx(const int64_t ord_idx, Vector3i const &lattice_max) {
  if (lattice_max(0) >= 0 && lattice_max(1) >= 0 && lattice_max(2) >= 0) {
    auto z = ord_idx % (2 * lattice_max(2) + 1);
    auto y = (ord_idx / (2 * lattice_max(2) + 1)) % (2 * lattice_max(1) + 1);
    auto x = ord_idx / (2 * lattice_max(2) + 1) / (2 * lattice_max(1) + 1);

    Vector3i result(x - lattice_max(0), y - lattice_max(1), z - lattice_max(2));
    return result;
  } else {
    throw "invalid lattice boundaries";
  }
}

std::shared_ptr<Molecule> shift_mol_origin(Molecule const &mol,
                                           Vector3d const &shift) {
  std::vector<AtomBasedClusterable> vec_of_clusters;
  for (auto &cluster : mol) {
    AtomBasedCluster shifted_cluster;
    for (auto &atom : collapse_to_atoms(cluster)) {
      Atom shifted_atom(atom.center() + shift, atom.mass(), atom.charge());
      shifted_cluster.add_clusterable(shifted_atom);
    }
    shifted_cluster.update_cluster();
    vec_of_clusters.emplace_back(shifted_cluster);
  }

  Molecule result(vec_of_clusters);

  auto result_ptr = std::make_shared<Molecule>(result);
  return result_ptr;
}

}  // namespace detail
}  // namespace mpqc
