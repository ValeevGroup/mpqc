#pragma once
#ifndef MPQC_BASIS_CLUSTERSHELLS_H
#define MPQC_BASIS_CLUSTERSHELLS_H

#include "../molecule/molecule_fwd.h"

#include <memory>
#include <vector>

namespace libint2 {
struct Shell;
}

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
                  std::shared_ptr<molecule::Cluster> c);

    std::vector<libint2::Shell> const &shells(unsigned int) const;
    std::vector<libint2::Shell> flattened_shells() const;

    unsigned int nfunctions(unsigned int) const;
    /// gives number of shells in am provided
    unsigned int nshells(unsigned int) const;
    /// Gives number of total shells
    unsigned int nshells() const;
    unsigned int flattened_nfunctions() const;

    unsigned int max_am() const;
    bool has_am(unsigned int) const;

    molecule::Cluster const &cluster() const;

  private:
    std::vector<std::vector<libint2::Shell>> shells_;
    std::shared_ptr<molecule::Cluster> cluster_;

}; // class ClusterShells

} // namespace basis
} // namespace mpqc

#endif // MPQC_BASIS_CLUSTERSHELLS_H
