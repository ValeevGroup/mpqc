
#ifndef MPQC_BASIS_CLUSTERSHELLS_H
#define MPQC_BASIS_CLUSTERSHELLS_H

#include "mpqc/chemistry/molecule/molecule_fwd.h"

#include <libint2/shell.h>

#include <memory>
#include <vector>

namespace mpqc {
namespace basis {

class ClusterShells {
 public:
  ClusterShells();
  ClusterShells(ClusterShells const &);
  ClusterShells(ClusterShells &&);
  ClusterShells &operator=(ClusterShells const &);
  ClusterShells &operator=(ClusterShells &&);

  enum class ang_mo { s = 0, p = 1, d = 2, f = 3, g = 4, h = 5, i = 6 };

  ClusterShells(std::vector<std::vector<libint2::Shell>> shell,
                std::shared_ptr<Cluster> c);

  std::vector<libint2::Shell> const &shells(unsigned int) const;
  std::vector<libint2::Shell> flattened_shells() const;

  unsigned int nfunctions(unsigned int) const;
  /// gives number of shells in am provided
  unsigned int nshells(unsigned int) const;

  /// Gives total number of total shells in this ClusterShells
  unsigned int nshells() const;

  /// total number of basis functions in this ClusterShells
  unsigned int flattened_nfunctions() const;

  unsigned int max_am() const;
  bool has_am(unsigned int) const;

  Cluster const &cluster() const;

 private:
  std::vector<std::vector<libint2::Shell>> shells_;
  std::shared_ptr<Cluster> cluster_;

};  // class ClusterShells

}  // namespace basis
}  // namespace mpqc

#endif  // MPQC_BASIS_CLUSTERSHELLS_H
