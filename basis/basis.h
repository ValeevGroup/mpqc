#pragma once
#ifndef TCC_BASIS_BASIS_H
#define TCC_BASIS_BASIS_H

#include "basis_fwd.h"
#include "../molecule/molecule_fwd.h"

#include <vector>

namespace tcc {
namespace basis {

class Basis {
  public:
    Basis();
    ~Basis();
    Basis(Basis const &);
    Basis(Basis &&);
    Basis &operator=(Basis const &);
    Basis &operator=(Basis &&);

    Basis(std::vector<ClusterShells> cs);

  private:
    std::vector<ClusterShells> cluster_shells_;
};

} // namespace basis
} // namespace tcc
#endif /* end of include guard: TCC_BASIS_BASIS_H */
