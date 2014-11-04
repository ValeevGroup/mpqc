#pragma once
#ifndef TCC_BASIS_CLUSTERSHELLS_H
#define TCC_BASIS_CLUSTERSHELLS_H

#include "../include/libint.h"
#include "../molecule/cluster.h"

#include <memory>
#include <vector>

namespace tcc {
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

    std::vector<libint2::Shell> const &shells(ang_mo am) const;
    molecule::Cluster const &cluster() const;

  private:
    std::vector<std::vector<libint2::Shell>> shells_;
    std::shared_ptr<molecule::Cluster> cluster_;

}; // class ClusterShells

} // namespace basis
} // namespace tcc

#endif // TCC_BASIS_CLUSTERSHELLS_H
