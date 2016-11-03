
#ifndef MPQC_BASIS_BASISSET_H
#define MPQC_BASIS_BASISSET_H

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include <libint2/shell.h>

#include "mpqc/chemistry/molecule/molecule_fwd.h"
#include "mpqc/chemistry/qc/basis/basis_fwd.h"



namespace mpqc {
namespace basis {

using Shell = libint2::Shell;
using ShellVec = std::vector<Shell>;

class BasisSet {
 public:
  typedef mpqc::basis::Shell Shell;

  BasisSet() = delete;  // Can't init a basis without name.
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
  std::vector<ShellVec> get_cluster_shells(Molecule const &) const;

  /*! \brief returns a single vector of all shells in the molecule */
  ShellVec get_flat_shells(Molecule const &) const;

 private:
  std::string basis_set_name_;
};

}  // namespace basis
}  // namespace mpqc

#endif  // MPQC_BASIS_BASISSET_H
