#include "mpqc/chemistry/qc/basis/basis_set.h"
#include "mpqc/chemistry/qc/basis/basis_set_maps.h"
#include "mpqc/chemistry/qc/basis/cluster_shells.h"

#include "mpqc/chemistry/molecule/common.h"
#include "mpqc/chemistry/molecule/molecule.h"

#include <libint2/basis.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

BasisSet::BasisSet(std::string const &s) : basis_set_name_{s} {}

std::vector<ShellVec> BasisSet::get_cluster_shells(Molecule const &mol) const {
  std::vector<ShellVec> cs;
  for (auto const &cluster : mol) {
    const auto libint_atoms = to_libint_atom(collapse_to_atoms(cluster));

    std::streambuf *cout_sbuf = std::cout.rdbuf();  // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    libint2::BasisSet libint_basis(basis_set_name_, libint_atoms);
    std::cout.rdbuf(cout_sbuf);

    // Shells that go with this cluster
    ShellVec cluster_shells;
    cluster_shells.reserve(libint_basis.size());

    for (auto &&shell : libint_basis) {
      cluster_shells.emplace_back(std::move(shell));
    }

    cs.emplace_back(std::move(cluster_shells));
  }

  return cs;
}

ShellVec BasisSet::get_flat_shells(Molecule const &mol) const {
  ShellVec cs;
  for (auto const &cluster : mol) {
    const auto libint_atoms = to_libint_atom(collapse_to_atoms(cluster));

    std::streambuf *cout_sbuf = std::cout.rdbuf();  // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    libint2::BasisSet libint_basis(basis_set_name_, libint_atoms);
    std::cout.rdbuf(cout_sbuf);

    for (auto &&shell : libint_basis) {
      cs.emplace_back(std::move(shell));
    }
  }

  return cs;
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
