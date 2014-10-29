#pragma once
#ifndef TILECLUSTERCHEM_BASIS_BASISSET_H
#define TILECLUSTERCHEM_BASIS_BASISSET_H

#include "../molecule/molecule_fwd.h"
#include "basis_fwd.h"

#include <string>
#include <vector>
#include <iosfwd>

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
    std::vector<ShellCluster> create_basis(molecule::Molecule const &) const;

  private:
    void read_basis(std::string const &s);

    std::vector<AtomBasisSet> atom_bases_;
};

std::ostream &operator<<(std::ostream &os, BasisSet const &bs);

} // namespace basis
} // namespace tcc

#endif // TILECLUSTERCHEM_BASIS_BASISSET_H
