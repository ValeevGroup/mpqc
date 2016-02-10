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

    Basis join(const Basis & basis);

    std::vector<ShellVec> const & cluster_shells() const;

    TiledArray::TiledRange1 create_trange1() const;

    int64_t max_nprim() const;
    int64_t max_am() const;
    int64_t nfunctions() const;
    int64_t nshells() const;
    int64_t nclusters() const { return shells_.size(); };
    std::vector<Shell> flattened_shells() const;

  private:
    std::vector<ShellVec> shells_;
};

std::ostream & operator<<(std::ostream &, Basis const &);

/*! \brief reblock allows for reblocking a basis
 *
 * \warning If reblocking a basis with the intent to use it with tensors 
 * computed with the old basis you must be careful not to reorder the shells. 
 *
 * \param op should be a function that takes a std::vector<Shell> and returns 
 * a std::vector<std::vector<Shell>> for use in initializing a Basis.
 */
template<typename Op, typename... Args>
Basis reblock(Basis const &basis, Op op, Args... args){
    return Basis(op(basis.flattened_shells(),args...));
}

} // namespace basis
} // namespace mpqc
#endif /* end of include guard: MPQC_BASIS_BASIS_H */
