#ifndef TILECLUSTERCHEM_BASIS_BASISSET_H
#define TILECLUSTERCHEM_BASIS_BASISSET_H

#include "atom_basisset.h"
#include <string>

class BasisSet {
  BasisSet() = default;
  BasisSet(BasisSet const &b) = default;
  BasisSet(BasisSet &&b) = default;
  BasisSet& operator=(BasisSet const &b) = default;
  BasisSet& operator=(BasisSet &&b) = default;

  BasisSet(std::string const&s) : atom_bases() {
    read_basis(s);
  }

private:
  void read_basis(std::string const &s) const;

  std::vector<AtomBasisSet> atom_bases;
};

#endif // TILECLUSTERCHEM_BASIS_BASISSET_H

