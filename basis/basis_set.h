#pragma once
#ifndef MPQC_BASIS_BASISSET_H
#define MPQC_BASIS_BASISSET_H

#include "../molecule/molecule_fwd.h"
#include "basis_fwd.h"

#include <string>
#include <vector>
#include <memory>
#include <iosfwd>

// FWD Decl
namespace libint2 {
struct Shell;
}

namespace mpqc {
namespace basis {

class BasisSet {
  public:
    BasisSet() = delete; // Can't init a basis without name.
    BasisSet(BasisSet const &b) = default;
    BasisSet(BasisSet &&b) = default;
    BasisSet &operator=(BasisSet const &b) = default;
    BasisSet &operator=(BasisSet &&b) = default;

    /// BasisSet takes a string which specifies the name of the basis set to
    /// use.
    BasisSet(std::string const &s);

    std::vector<ClusterShells> create_basis(
        std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const;

    std::vector<ClusterShells> create_soad_basis(
        std::vector<std::shared_ptr<molecule::Cluster>> const &clusters) const;

  private:
    std::string basis_set_name_;
};

} // namespace basis
} // namespace mpqc

#endif // MPQC_BASIS_BASISSET_H
