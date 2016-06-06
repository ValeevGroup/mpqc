#pragma once
#ifndef MPQC_BASIS_BASISSET_H
#define MPQC_BASIS_BASISSET_H

#include <mpqc/chemistry/molecule/molecule_fwd.h>
#include <mpqc/chemistry/qc/basis/basis_fwd.h>

#include "../../../../../common/namespaces.h"
#include "../../../../../common/typedefs.h"

#include <string>
#include <vector>
#include <memory>
#include <iosfwd>


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

    /*! \brief Constructs a vector of ShellVecs
     *
     * Each ShellVec represents the shells for a cluster.
     */
    std::vector<ShellVec> get_cluster_shells(mol::Molecule const &) const;


    /*! \brief returns a single vector of all shells in the molecule */
    ShellVec get_flat_shells(mol::Molecule const &) const;

  private:
    std::string basis_set_name_;
};

} // namespace basis
} // namespace mpqc

#endif // MPQC_BASIS_BASISSET_H
