#pragma once
#ifndef MPQC_BASIS_BASIS_H
#define MPQC_BASIS_BASIS_H

#include "basis_fwd.h"
#include "../molecule/molecule_fwd.h"

#include <vector>
#include <iosfwd>

#include "../common/typedefs.h"


namespace mpqc {
namespace basis {

class Basis {
  public:
    Basis();
    ~Basis();
    Basis(Basis const &);
    Basis(Basis &&);
    Basis &operator=(Basis const &);
    Basis &operator=(Basis &&);

    Basis(std::vector<ShellVec> cs);

    std::vector<ShellVec> const & cluster_shells() const;

    TiledArray::TiledRange1 create_trange1() const;

    int64_t max_nprim() const;
    int64_t max_am() const;
    int64_t nfunctions() const;

  private:
    std::vector<ShellVec> shells_;
};

std::ostream & operator<<(std::ostream &, Basis const &);

} // namespace basis
} // namespace mpqc
#endif /* end of include guard: MPQC_BASIS_BASIS_H */
