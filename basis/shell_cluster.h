#pragma once
#ifndef TCC_BASIS_SHELLCLUSTER_H
#define TCC_BASIS_SHELLCLUSTER_H

#include "../include/libint.h"
#include "molecule/cluster.h"

#include <memory>
#include <vector>

namespace tcc {
namespace basis {

class ShellCluster {
  public:
    ShellCluster();
    ShellCluster(ShellCluster const &);
    ShellCluster(ShellCluster &&);
    ShellCluster &operator=(ShellCluster const &);
    ShellCluster &operator=(ShellCluster &&);

    ShellCluster(std::vector<libint2::Shell> shell, molecule::Cluster c);

    std::vector<libint2::Shell> const &shells() const; 
    molecule::Cluster const &cluster() const; 

  private:
    std::vector<libint2::Shell> shells_;
    std::shared_ptr<molecule::Cluster> cluster_;
}; // class ShellCluster

} // namespace basis
} // namespace tcc

#endif // TCC_BASIS_SHELLCLUSTER_H
