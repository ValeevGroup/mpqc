#pragma once
#ifndef TCC_BASIS_BASIS_H
#define TCC_BASIS_BASIS_H

#include "basis_fwd.h"
#include "../molecule/molecule_fwd.h"

#include <vector>
#include <iosfwd>

namespace TiledArray { class TiledRange1; }

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

    std::vector<ClusterShells> const & cluster_shells() const;

    TiledArray::TiledRange1 create_flattend_trange1() const;
    TiledArray::TiledRange1 create_trange1() const;

  private:
    std::vector<ClusterShells> cluster_shells_;

    std::vector<unsigned int> am_blocking_generator() const;
    std::vector<unsigned int> flattened_blocking_generator() const;
};

std::ostream & operator<<(std::ostream &, Basis const &);

} // namespace basis
} // namespace tcc
#endif /* end of include guard: TCC_BASIS_BASIS_H */
