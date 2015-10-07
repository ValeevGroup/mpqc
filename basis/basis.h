#pragma once
#ifndef MPQC_BASIS_BASIS_H
#define MPQC_BASIS_BASIS_H

#include "basis_fwd.h"
#include "../molecule/molecule_fwd.h"

#include <vector>
#include <iosfwd>
#include <memory>

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
    Basis(std::vector<ShellVec> cs, molecule::Molecule const&);
    Basis(std::vector<ShellVec> cs, std::shared_ptr<molecule::Molecule> mol_ptr);

    std::vector<ShellVec> const & cluster_shells() const;

    // Returns the molecule, potentially will throw if no Molecule is 
    // avalible. 
    molecule::Molecule const &molecule() const;

    TiledArray::TiledRange1 create_trange1() const;

    int64_t max_nprim() const;
    int64_t max_am() const;
    int64_t nfunctions() const;
    int64_t nshells() const;
    int64_t nclusters() const { return shells_.size(); };

  private:
    std::vector<ShellVec> shells_;
    std::shared_ptr<molecule::Molecule> mol_;
};

std::ostream & operator<<(std::ostream &, Basis const &);

} // namespace basis
} // namespace mpqc
#endif /* end of include guard: MPQC_BASIS_BASIS_H */
