#include "basis_set.h"
#include "basis_set_maps.h"
#include "atom_basisset.h"
#include "cluster_shells.h"

#include "../molecule/cluster_collapse.h"
#include "../molecule/common.h"
#include "../molecule/cluster.h"
#include "../molecule/atom.h"

#include "../include/libint.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace tcc {
namespace basis {

BasisSet::BasisSet(std::string const &s) : basis_set_name_{s} {}

std::vector<ClusterShells> BasisSet::create_basis(
    std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const {

    using namespace molecule;

    std::vector<ClusterShells> cs;
    for (auto const &c : clusters) {

        const auto libint_atoms = to_libint_atom(collapse_to_atoms(*c));

        // Sneaky Libint2::Basis inherits from std::vector<libint2::Shell> !
        libint2::BasisSet libint_basis(basis_set_name_, libint_atoms); 

        // Create vector of shells binned by angular momentum need 1 because 
        // s shells have number 0.
        const auto n_am = libint_basis.max_l()+1; 
        std::vector<std::vector<libint2::Shell>> binned_shells(n_am);

        for (auto const &shell : libint_basis) {
            const auto ang_mo = shell.contr[0].l;
            binned_shells[ang_mo].push_back(shell);
        }

        cs.emplace_back(std::move(binned_shells), c);
    }

    return cs;
}

std::vector<ClusterShells> BasisSet::create_soad_basis(
    std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const {

    using namespace molecule;

    std::vector<ClusterShells> cs;
    for (auto const &c : clusters) {

        const auto libint_atoms = to_libint_atom(collapse_to_atoms(*c));

        // Sneaky Libint2::Basis inherits from std::vector<libint2::Shell> !
        libint2::BasisSet libint_basis(basis_set_name_, libint_atoms); 

        // For soad just wrap shells once
        std::vector<std::vector<libint2::Shell>> binned_shells(1);

        for (auto const &shell : libint_basis) {
            binned_shells[0].push_back(shell);
        }

        cs.emplace_back(std::move(binned_shells), c);
    }

    return cs;
}

} // namespace basis
} // namespace tcc
