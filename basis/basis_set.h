#pragma once
#ifndef TILECLUSTERCHEM_BASIS_BASISSET_H
#define TILECLUSTERCHEM_BASIS_BASISSET_H

#include "../molecule/molecule_fwd.h"
#include "basis_fwd.h"

#include <string>
#include <vector>
#include <iosfwd>

// FWD Decl
namespace libint2 {
struct Shell;
}

namespace tcc {
namespace basis {

class BasisSet {
  public:
    BasisSet();
    BasisSet(BasisSet const &b);
    BasisSet(BasisSet &&b);
    BasisSet &operator=(BasisSet const &b);
    BasisSet &operator=(BasisSet &&b);

    BasisSet(std::string const &s);

    std::vector<AtomBasisSet> const & atom_basis_set() const;
    std::vector<libint2::Shell> atom_basis(molecule::Atom const &a) const;
    
    std::vector<ClusterShells> create_basis(
        std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const;

  private:
    void read_basis(std::string const &s);

    std::vector<AtomBasisSet> atom_bases_;
};

std::ostream &operator<<(std::ostream &os, BasisSet const &bs);

} // namespace basis
} // namespace tcc

#endif // TILECLUSTERCHEM_BASIS_BASISSET_H
