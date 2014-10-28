#ifndef TILECLUSTERCHEM_BASIS_BASISSET_H
#define TILECLUSTERCHEM_BASIS_BASISSET_H

#include "atom_basisset.h"
#include <string>
#include <iosfwd>

namespace tcc {
namespace basis {

class BasisSet {
  public:
    BasisSet() = default;
    BasisSet(BasisSet const &b) = default;
    BasisSet(BasisSet &&b) = default;
    BasisSet &operator=(BasisSet const &b) = default;
    BasisSet &operator=(BasisSet &&b) = default;

    BasisSet(std::string const &s) : atom_bases_() { read_basis(s); }

    std::vector<AtomBasisSet> const & basis() const {
        return atom_bases_;
    }

  private:
    void read_basis(std::string const &s);

    std::vector<AtomBasisSet> atom_bases_;
};

std::ostream &operator<<(std::ostream &os, BasisSet const &bs);

} // namespace basis
} // namespace tcc

#endif // TILECLUSTERCHEM_BASIS_BASISSET_H
