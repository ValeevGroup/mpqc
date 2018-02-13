#include "mpqc/chemistry/molecule/lattice/util.h"

#include "mpqc/util/core/exception.h"
#include "mpqc/util/misc/assert.h"

namespace mpqc {
namespace detail {

Vector3d direct_vector(int64_t ord_idx, Vector3i const &lattice_max,
                       Vector3d const &dcell) {
  if ((lattice_max.array() < 0).any()) {
    throw ProgrammingError("invalid lattice range", __FILE__, __LINE__);
  }
  auto z = ord_idx % (2 * lattice_max(2) + 1);
  auto y = (ord_idx / (2 * lattice_max(2) + 1)) % (2 * lattice_max(1) + 1);
  auto x = ord_idx / (2 * lattice_max(2) + 1) / (2 * lattice_max(1) + 1);
  Vector3d result((x - lattice_max(0)) * dcell(0),
                  (y - lattice_max(1)) * dcell(1),
                  (z - lattice_max(2)) * dcell(2));
  return result;
}

Vector3d k_vector(int64_t ord_idx, Vector3i const &nk, Vector3d const &dcell) {
  if ((nk.array() < 1).any()) {
    throw ProgrammingError("invalid number of k points", __FILE__, __LINE__);
  }

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
  if ((lattice_max.array() >= 0).all() && std::abs(x) <= lattice_max(0) &&
      std::abs(y) <= lattice_max(1) && std::abs(z) <= lattice_max(2)) {
    int64_t idx = (x + lattice_max(0)) * (2 * lattice_max(1) + 1) *
                      (2 * lattice_max(2) + 1) +
                  (y + lattice_max(1)) * (2 * lattice_max(2) + 1) +
                  (z + lattice_max(2));
    return idx;
  } else {
    throw ProgrammingError("invalid lattice sum index/boundaries", __FILE__,
                           __LINE__);
  }
}

int64_t k_ord_idx(Vector3i const &in_3D_idx, Vector3i const &nk) {
  return k_ord_idx(in_3D_idx(0), in_3D_idx(1), in_3D_idx(2), nk);
}

int64_t k_ord_idx(int64_t x, int64_t y, int64_t z, Vector3i const &nk) {
  if ((nk.array() >= 1).all() && x >= 0 && y >= 0 && z >= 0 && x < nk(0) &&
      y < nk(1) && z < nk(2)) {
    int64_t idx = x * nk(1) * nk(2) + y * nk(2) + z;
    return idx;
  } else {
    throw ProgrammingError("invalid k-space index/boundaries", __FILE__,
                           __LINE__);
  }
}

Vector3i direct_3D_idx(const int64_t ord_idx, Vector3i const &lattice_max) {
  if ((lattice_max.array() >= 0).all()) {
    auto z = ord_idx % (2 * lattice_max(2) + 1);
    auto y = (ord_idx / (2 * lattice_max(2) + 1)) % (2 * lattice_max(1) + 1);
    auto x = ord_idx / (2 * lattice_max(2) + 1) / (2 * lattice_max(1) + 1);

    Vector3i result(x - lattice_max(0), y - lattice_max(1), z - lattice_max(2));
    return result;
  } else {
    throw ProgrammingError("invalid lattice boundaries", __FILE__, __LINE__);
  }
}

Vector3i k_3D_idx(const int64_t ord_idx, const Vector3i &nk) {
  MPQC_ASSERT((nk.array() > 0).all());

  auto z = ord_idx % nk(2);
  auto xy = ord_idx / nk(2);
  auto y = xy % nk(1);
  auto x = xy / nk(1);

  // assume odd number of k points in each direction (e.g. uniformly spaced k)
  MPQC_ASSERT((nk(0) % 2 == 1) && (nk(1) % 2 == 1) && (nk(2) % 2 == 1));

  Vector3i ref((nk(0) - 1) / 2, (nk(1) - 1) / 2, (nk(2) - 1) / 2);
  return Vector3i(x - ref(0), y - ref(1), z - ref(2));
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

/*!
 * \brief This determines if a unit cell is included by the give lattice range
 * \param in_idx 3D index of a unit cell
 * \param range lattice range
 * \param center center of \range
 * \return
 */
bool is_in_lattice_range(Vector3i const &in_idx, Vector3i const &range,
                         Vector3i const &center) {
  if (in_idx(0) <= center(0) + range(0) && in_idx(0) >= center(0) - range(0) &&
      in_idx(1) <= center(1) + range(1) && in_idx(1) >= center(1) - range(1) &&
      in_idx(2) <= center(2) + range(2) && in_idx(2) >= center(2) - range(2))
    return true;
  else
    return false;
}

}  // namespace detail
}  // namespace mpqc
