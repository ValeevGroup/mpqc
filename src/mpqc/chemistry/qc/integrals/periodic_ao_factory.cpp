#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

namespace mpqc {
namespace integrals {
namespace pbc {

TA::TiledRange1 extend_trange1(TiledArray::TiledRange1 tr0, int64_t size) {
    auto blocking = std::vector<int64_t>{0};
    for (auto idx = 0; idx < size; ++idx) {
      for (auto u = 0; u < tr0.tile_extent(); ++u) {
        auto next = blocking.back() + tr0.tile(u).second - tr0.tile(u).first;
        blocking.emplace_back(next);
      }
    }
    TA::TiledRange1 tr1(blocking.begin(), blocking.end());
    return tr1;
}

void sort_eigen(Vectorc &eigVal, Matrixc &eigVec) {
    auto val = eigVal.real();

    // Sort by ascending eigenvalues
    std::vector<std::pair<double, int>> sortedVal;
    sortedVal.reserve(val.size());
    for (auto i = 0; i != val.size(); ++i) {
      auto pair = std::make_pair(val(i), i);
      sortedVal.push_back(pair);
    };
    std::sort(sortedVal.begin(), sortedVal.end());

    // Build sorted eigenvalues and eigenvectors
    Vectorc sortedEigVal(eigVal);
    Matrixc sortedEigVec(eigVec);
    for (auto i = 0; i != val.size(); ++i) {
      sortedEigVal(i) = eigVal(sortedVal[i].second);
      sortedEigVec.col(i) = eigVec.col(sortedVal[i].second);
    }

    eigVal = sortedEigVal;
    eigVec = sortedEigVec;
}

Vector3d R_vector(int64_t idx_lattice, Vector3i vec, Vector3d dcell) {
    auto z = idx_lattice % (2 * vec(2) + 1);
    auto y = (idx_lattice / (2 * vec(2) + 1)) % (2 * vec(1) + 1);
    auto x = idx_lattice / (2 * vec(2) + 1) / (2 * vec(1) + 1);
    Vector3d result((x - vec(0)) * dcell(0), (y - vec(1)) * dcell(1),
                    (z - vec(2)) * dcell(2));
    return result;

}

Vector3d k_vector(int64_t idx_k, Vector3i nk, Vector3d dcell) {
    Vector3d result;
    auto x = idx_k / nk(2) / nk(1);
    auto y = (idx_k / nk(2)) % nk(1);
    auto z = idx_k % nk(2);
    result(0) = (dcell(0) == 0.0) ? 0.0
                                   : (-1.0 + (2.0 * (x + 1) - 1.0) / nk(0)) *
                                         (M_PI / dcell(0));
    result(1) = (dcell(1) == 0.0) ? 0.0
                                   : (-1.0 + (2.0 * (y + 1) - 1.0) / nk(1)) *
                                         (M_PI / dcell(1));
    result(2) = (dcell(2) == 0.0) ? 0.0
                                   : (-1.0 + (2.0 * (z + 1) - 1.0) / nk(2)) *
                                         (M_PI / dcell(2));
    return result;
}

int64_t idx_lattice(int x, int y, int z, Vector3i vec) {
    if (vec(0) >= 0 && vec(1) >= 0 && vec(2) >= 0 && abs(x) <= vec(0) &&
        abs(y) <= vec(1) && abs(z) <= vec(2)) {
      int64_t idx = (x + vec(0)) * (2 * vec(0) + 1) * (2 * vec(1) + 1) +
                    (y + vec(1)) * (2 * vec(1) + 1) + (z + vec(2));
      return idx;
    } else {
      throw "invalid lattice sum index/boundaries";
    }

}

int64_t idx_k(int x, int y, int z, Vector3i nk) {
    if (nk(0) >= 1 && nk(1) >= 1 && nk(2) >= 1 && x >= 0 && y >= 0 && z >= 0 &&
        x < nk(0) && y < nk(1) && z < nk(2)) {
      int64_t idx = x * nk(0) * nk(1) + y * nk(1) + z;
      return idx;
    } else {
      throw "invalid k-space index/boundaries";
    }
}

std::shared_ptr<basis::Basis> shift_basis_origin(basis::Basis &basis,
                                                 Vector3d shift) {
    std::vector<ShellVec> vec_of_shells;
    for (auto shell_vec : basis.cluster_shells()) {
      ShellVec shells;
      for (auto shell : shell_vec) {
        std::array<double, 3> new_origin = {
            shell.O[0] + shift(0), shell.O[1] + shift(1), shell.O[2] + shift(2)};
        shell.move(new_origin);
        shells.push_back(shell);
      }
      vec_of_shells.push_back(shells);
    }
    basis::Basis result(vec_of_shells);
    auto result_ptr = std::make_shared<basis::Basis>(result);
    return result_ptr;
}

std::shared_ptr<basis::Basis> shift_basis_origin(basis::Basis &basis,
                                                 Vector3d shift_base,
                                                 Vector3i nshift,
                                                 Vector3d dcell,
                                                 bool is_real_space) {
    std::vector<ShellVec> vec_of_shells;

    int64_t shift_size =
        is_real_space
            ? (1 + idx_lattice(nshift(0), nshift(1), nshift(2), nshift))
            : (1 + idx_k(nshift(0) - 1, nshift(1) - 1, nshift(2) - 1, nshift));

    for (auto idx_shift = 0; idx_shift < shift_size; ++idx_shift) {
      Vector3d shift = is_real_space ? (R_vector(idx_shift, nshift, dcell) + shift_base)
                                     : (k_vector(idx_shift, nshift, dcell) + shift_base);

      for (auto shell_vec : basis.cluster_shells()) {
        ShellVec shells;
        for (auto shell : shell_vec) {
          std::array<double, 3> new_origin = {shell.O[0] + shift(0),
                                              shell.O[1] + shift(1),
                                              shell.O[2] + shift(2)};
          shell.move(new_origin);
          shells.push_back(shell);
        }
        vec_of_shells.push_back(shells);
      }
    }

    basis::Basis result(vec_of_shells);
    auto result_ptr = std::make_shared<basis::Basis>(result);
    return result_ptr;

}

std::shared_ptr<Molecule> shift_mol_origin(const Molecule &mol, Vector3d shift) {
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

}  // namespace pbc
}  // namespace integrals
}  // namespace mpqc
